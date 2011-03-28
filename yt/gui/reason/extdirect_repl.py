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

from .bottle_mods import preroute, BottleDirectRouter, notify_route, \
                         PayloadHandler
from .bottle import response, request
from .basic_repl import ProgrammaticREPL

try:
    import pygments
    import pygments.lexers
    import pygments.formatters
    def _highlighter():
        pl = pygments.lexers.PythonLexer()
        hf = pygments.formatters.HtmlFormatter(linenos='table')
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
    _skip_expose = ('index', 'resources')
    my_name = "ExtDirectREPL"

    def __init__(self, locals=None):
        # First we do the standard initialization
        ProgrammaticREPL.__init__(self, locals)
        # Now, since we want to only preroute functions we know about, and
        # since they have different arguments, and most of all because we only
        # want to add them to the routing tables (which are a singleton for the
        # entire interpreter state) we apply all the pre-routing now, rather
        # than through metaclasses or other fancy decorating.
        preroute_table = dict(index = ("/", "GET"),
                              _myapi = ("/resources/ext-repl-api.js", "GET"),
                              resources = ("/resources/:path#.+#", "GET"),
                              _session_py = ("/session.py", "GET"),
                              _ace = ("/ace/:path#.+#", "GET"),
                              _highlighter_css = ("/highlighter.css", "GET"),
                              )
        for v, args in preroute_table.items():
            preroute(args[0], method=args[1])(getattr(self, v))
        self.api_url = "repl"
        BottleDirectRouter.__init__(self)
        self.pflist = ExtDirectParameterFileList()
        self.executed_cell_texts = []
        self.payload_handler = PayloadHandler()
        # Now we load up all the yt.mods stuff, but only after we've finished
        # setting up.
        self.execute("from yt.mods import *")
        self.locals['load_script'] = ext_load_script

    def index(self):
        """Return an HTTP-based Read-Eval-Print-Loop terminal."""
        # For now this doesn't work!  We will need to move to a better method
        # for this.
        vals = open(os.path.join(local_dir, "html/index.html")).read()
        return vals

    def resources(self, path):
        # This will need to be changed.
        pp = os.path.join(local_dir, "../../../../misc/ext/ext-3.3.1/", path)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

    def _ace(self, path):
        # This will need to be changed.
        pp = os.path.join(local_dir, "ace", path)
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

    def _session_py(self):
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs.seek(0)
        response.headers["content-disposition"] = "attachment; filename=session.py"
        return cs

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
