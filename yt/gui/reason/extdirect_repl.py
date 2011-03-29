"""
A read-eval-print-loop that is served up through Bottle and accepts its
commands through ExtDirect calls

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import json
import os
import cStringIO
import logging
import uuid
from yt.utilities.logger import ytLogger, ufstring

from .bottle_mods import preroute, BottleDirectRouter, notify_route, \
                         PayloadHandler
from .bottle import response, request, route
from .basic_repl import ProgrammaticREPL

try:
    import pygments
    import pygments.lexers
    import pygments.formatters
    def _highlighter():
        pl = pygments.lexers.PythonLexer()
        hf = pygments.formatters.HtmlFormatter(
                linenos='table', linenospecial=2)
        def __highlighter(a):
            return pygments.highlight(a, pl, hf)
        return __highlighter, hf.get_style_defs()
    # We could add an additional '.highlight_pyg' in the call
    highlighter, highlighter_css = _highlighter()
except ImportError:
    highlighter = lambda a: a
    highlight_css = ''

local_dir = os.path.dirname(__file__)

class ExtDirectREPL(ProgrammaticREPL, BottleDirectRouter):
    _skip_expose = ('index')
    my_name = "ExtDirectREPL"

    def __init__(self, base_extjs_path, locals=None):
        # First we do the standard initialization
        self.extjs_path = os.path.join(base_extjs_path, "ext-resources")
        self.extjs_theme_path = os.path.join(base_extjs_path, "ext-theme")
        ProgrammaticREPL.__init__(self, locals)
        # Now, since we want to only preroute functions we know about, and
        # since they have different arguments, and most of all because we only
        # want to add them to the routing tables (which are a singleton for the
        # entire interpreter state) we apply all the pre-routing now, rather
        # than through metaclasses or other fancy decorating.
        preroute_table = dict(index = ("/", "GET"),
                              _myapi = ("/resources/ext-repl-api.js", "GET"),
                              _resources = ("/resources/:path#.+#", "GET"),
                              _js = ("/js/:path#.+#", "GET"),
                              _images = ("/images/:path#.+#", "GET"),
                              _theme = ("/theme/:path#.+#", "GET"),
                              _session_py = ("/session.py", "GET"),
                              _highlighter_css = ("/highlighter.css", "GET"),
                              )
        for v, args in preroute_table.items():
            preroute(args[0], method=args[1])(getattr(self, v))
        # This has to be routed to the root directory
        self.api_url = "repl"
        BottleDirectRouter.__init__(self)
        self.pflist = ExtDirectParameterFileList()
        self.executed_cell_texts = []
        self.payload_handler = PayloadHandler()
        # Now we load up all the yt.mods stuff, but only after we've finished
        # setting up.
        self.execute("from yt.mods import *")
        self.locals['load_script'] = ext_load_script
        self.locals['_widgets'] = {}
        self.locals['add_widget'] = self._add_widget
        self.locals['test_widget'] = self._test_widget
        self._setup_logging_handlers()

    def _setup_logging_handlers(self):
        handler = PayloadLoggingHandler()
        formatter = logging.Formatter(ufstring)
        handler.setFormatter(formatter)
        ytLogger.addHandler(handler)

    def index(self):
        """Return an HTTP-based Read-Eval-Print-Loop terminal."""
        # For now this doesn't work!  We will need to move to a better method
        # for this.
        vals = open(os.path.join(local_dir, "html/index.html")).read()
        return vals

    def _resources(self, path):
        pp = os.path.join(self.extjs_path, path)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

    def _theme(self, path):
        pp = os.path.join(self.extjs_theme_path, path)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

    def _js(self, path):
        pp = os.path.join(local_dir, "html", "js", path)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

    def _images(self, path):
        pp = os.path.join(local_dir, "html", "images", path)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

    def _highlighter_css(self):
        return highlighter_css

    def execute(self, code):
        self.executed_cell_texts.append(code)
        result = ProgrammaticREPL.execute(self, code)
        payloads = self.payload_handler.deliver_payloads()
        return_value = {'output': result,
                        'input': highlighter(code),
                        'payloads': payloads}
        return return_value

    def get_history(self):
        return self.executed_cell_texts[:]

    def save_session(self, filename):
        if not filename.startswith('/'):
            filename = os.path.join(os.path.expanduser('~/'), filename)
        if os.path.exists(filename):
            return {'status': 'FAIL', 'filename': filename,
                    'error': 'File exists!'}
        try:
            f = open(filename, 'w')
            f.write("\n######\n".join(self.executed_cell_texts))
            f.close()
        except IOError as (errno, strerror):
            return {'status': 'FAIL', 'filename': filename,
                    'error': strerror}
        except:
            return {'status': 'FAIL', 'filename': filename,
                    'error': 'Unexpected error.'}
        return {'status': 'SUCCESS', 'filename': filename}

    def _session_py(self):
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs.seek(0)
        response.headers["content-disposition"] = "attachment;"
        return cs

    def _add_widget(self, widget):
        # This should be sanitized
        uu = str(uuid.uuid1()).replace("-","_")
        varname = "%s_%s" % (widget._widget_name, uu)
        widget._ext_widget_id = varname
        self.locals[varname] = widget
        payload = {'type': 'widget',
                   'widget_type': widget._widget_name,
                   'varname': varname}
        print payload
        self.payload_handler.add_payload(payload)

    def _test_widget(self):
        class tt(object):
            _widget_name = "plot_window"
        mm = tt()
        return mm

class ExtDirectParameterFileList(BottleDirectRouter):
    my_name = "ExtDirectParameterFileList"
    api_url = "pflist"

    def get_list_of_pfs(self):
        from yt.data_objects.static_output import _cached_pfs
        rv = []
        for fn, pf in sorted(_cached_pfs.items()):
            objs = []
            for obj in pf.h.objects:
                try:
                    name = str(obj)
                except ReferenceError:
                    continue
                objs.append(dict(name=name, type=obj._type_name))
            rv.append( dict(name = str(pf), objects = objs) )
        return rv

def ext_load_script(filename):
    contents = open(filename).read()
    payload_handler = PayloadHandler()
    payload_handler.add_payload(
        {'type': 'cell_contents',
         'value': contents}
    )
    return

class PayloadLoggingHandler(logging.StreamHandler):
    def __init__(self, *args, **kwargs):
        logging.StreamHandler.__init__(self, *args, **kwargs)
        self.payload_handler = PayloadHandler()

    def emit(self, record):
        msg = self.format(record)
        self.payload_handler.add_payload(
            {'type':'log_entry',
             'log_entry':msg})

if os.path.exists(os.path.expanduser("~/.yt/favicon.ico")):
    ico = os.path.expanduser("~/.yt/favicon.ico")
else:
    ico = os.path.join(local_dir, "html", "images", "favicon.ico")
@route("/favicon.ico", method="GET")
def _favicon_ico():
    response.headers['Content-Type'] = "image/x-icon"
    return open(ico).read()


