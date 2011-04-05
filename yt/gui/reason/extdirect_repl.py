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
import numpy as na

from yt.funcs import *
from yt.utilities.logger import ytLogger, ufstring
from yt.utilities.definitions import inv_axis_names

from .bottle_mods import preroute, BottleDirectRouter, notify_route, \
                         PayloadHandler, append_payloads
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
                              _help_html = ("/help.html", "GET"),
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
        self.execute("from yt.data_objects.static_output import _cached_pfs")
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
        # for this.  It should use the package data command.
        vals = open(os.path.join(local_dir, "html/index.html")).read()
        return vals

    def heartbeep(self):
        return {'alive': True}

    def _help_html(self):
        vals = open(os.path.join(local_dir, "html/help.html")).read()
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

    @append_payloads
    def execute(self, code):
        self.executed_cell_texts.append(code)
        result = ProgrammaticREPL.execute(self, code)
        self.payload_handler.add_payload(
            {'type': 'cell_results',
             'output': result,
             'input': highlighter(code)})

    def get_history(self):
        return self.executed_cell_texts[:]

    def save_session(self, filename):
        if filename.startswith('~'):
            filename = os.path.expanduser(filename)
        elif not filename.startswith('/'):
            filename = os.path.join(os.getcwd(), filename)
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

    def paste_session(self):
        import xmlrpclib, cStringIO
        p = xmlrpclib.ServerProxy(
            "http://paste.enzotools.org/xmlrpc/",
            allow_none=True)
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs = cs.getvalue()
        ret = p.pastes.newPaste('pytb', cs, None, '', '', True)
        site = "http://paste.enzotools.org/show/%s" % ret
        return {'status': 'SUCCESS', 'site': site}

    def _session_py(self):
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs.seek(0)
        response.headers["content-disposition"] = "attachment;"
        return cs

    @append_payloads
    def _add_widget(self, widget_name):
        # This should be sanitized
        widget = self.locals[widget_name]
        uu = str(uuid.uuid1()).replace("-","_")
        varname = "%s_%s" % (widget._widget_name, uu)
        widget._ext_widget_id = varname
        # THIS BREAKS THE SCRIPT DOWNLOAD!
        # We need to make the variable be bound via an execution mechanism
        payload = {'type': 'widget',
                   'widget_type': widget._widget_name,
                   'varname': varname}
        self.payload_handler.add_payload(payload)
        self.execute("%s = %s\n" % (varname, widget_name))

    @append_payloads
    def create_proj(self, pfname, axis, field, weight):
        if weight == "None": weight = None
        else: weight = "'%s'" % (weight)
        funccall = """
        _tpf = %(pfname)s
        _taxis = %(axis)s
        _tfield = "%(field)s"
        _tweight = %(weight)s
        _tsl = _tpf.h.proj(_taxis,_tfield, weight_field=_tweight)
        _txax, _tyax = x_dict[_taxis], y_dict[_taxis]
        DLE, DRE = _tpf.domain_left_edge, _tpf.domain_right_edge
        from yt.visualization.plot_window import PWViewerExtJS
        _tpw = PWViewerExtJS(_tsl, (DLE[_txax], DRE[_txax], DLE[_tyax], DRE[_tyax]), setup = False)
        _tpw._current_field = _tfield
        _tpw.set_log(_tfield, True)
        add_widget('_tpw')
        """ % dict(pfname = pfname,
                   axis = inv_axis_names[axis],
                   weight = weight,
                   field=field)
        # There is a call to do this, but I have forgotten it ...
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall)

    @append_payloads
    def create_slice(self, pfname, center, axis, field):
        funccall = """
        _tpf = %(pfname)s
        _tcenter = na.array([%(c1)0.20f, %(c2)0.20f, %(c3)0.20f], dtype='float64')
        _taxis = %(axis)s
        _tfield = "%(field)s"
        _tcoord = _tcenter[_taxis]
        _tsl = _tpf.h.slice(_taxis, _tcoord, center = _tcenter)
        _txax, _tyax = x_dict[_taxis], y_dict[_taxis]
        DLE, DRE = _tpf.domain_left_edge, _tpf.domain_right_edge
        from yt.visualization.plot_window import PWViewerExtJS
        _tpw = PWViewerExtJS(_tsl, (DLE[_txax], DRE[_txax], DLE[_tyax], DRE[_tyax]), setup = False)
        _tpw._current_field = _tfield
        _tpw.set_log(_tfield, True)
        add_widget('_tpw')
        """ % dict(pfname = pfname,
                   c1 = float(center[0]),
                   c2 = float(center[1]),
                   c3 = float(center[2]),
                   axis = inv_axis_names[axis],
                   field=field)
        # There is a call to do this, but I have forgotten it ...
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall)

    def _test_widget(self):
        class tt(object):
            _widget_name = "plot_window"
        mm = tt()
        return mm

class ExtDirectParameterFileList(BottleDirectRouter):
    my_name = "ExtDirectParameterFileList"
    api_url = "pflist"

    def get_list_of_pfs(self):
        # Note that this instantiates the hierarchy.  This can be a costly
        # event.  However, we're going to assume that it's okay, if you have
        # decided to load up the parameter file.
        from yt.data_objects.static_output import _cached_pfs
        rv = []
        for fn, pf in sorted(_cached_pfs.items()):
            objs = []
            pf_varname = "_cached_pfs['%s']" % (fn)
            fields = set(pf.h.field_list + pf.h.derived_field_list)
            fields = list(fields)
            fields.sort()
            for i,obj in enumerate(pf.h.objects):
                try:
                    name = str(obj)
                except ReferenceError:
                    continue
                objs.append(dict(name=name, type=obj._type_name,
                                 varname = "%s.h.objects[%s]" % (pf_varname, i)))
            rv.append( dict(name = str(pf), objects = objs, filename=fn,
                            varname = pf_varname, field_list = fields) )
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


