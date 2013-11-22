"""
A read-eval-print-loop that is served up through Bottle and accepts its
commands through HTTP



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import json
import os

from .bottle_mods import preroute
from .basic_repl import ProgrammaticREPL

local_dir = os.path.dirname(__file__)

class HTTPREPL(ProgrammaticREPL):

    def __init__(self, locals=None):
        # First we do the standard initialization
        ProgrammaticREPL.__init__(self, locals)
        # Now, since we want to only preroute functions we know about, and
        # since they have different arguments, and most of all because we only
        # want to add them to the routing tables (which are a singleton for the
        # entire interpreter state) we apply all the pre-routing now, rather
        # than through metaclasses or other fancy decorating.
        preroute_table = dict(index = ("/", "GET"),
                              push = ("/push", "POST"),
                              dir = ("/dir", "GET"),
                              doc = ("/doc", "GET"),
                              resources = ("/resources/:val", "GET"))
        for v, args in preroute_table.items():
            preroute(args[0], method=args[1])(getattr(self, v))

    def index(self):
        """Return an HTTP-based Read-Eval-Print-Loop terminal."""
        # For now this doesn't work!  We will need to move to a better method
        # for this.
        vals = open(os.path.join(local_dir, "httprepl.html")).read()
        return vals
        
    def push(self):
        """Push 'line' and return exec results as a bare response."""
        line = request.POST['line']
        result = ProgrammaticREPL.push(self, line)
        new_values = self.locals.pop("new_values", "")
        if result is None:
            # More input lines needed.
            response.status = 204
        return json.dumps( dict(text = result, new_values = new_values ))

    def dir(self):
        """Push 'line' and return result of eval on the final expr."""
        line = request.GET['line']
        result = ProgrammaticREPL.dir(self, line)
        if not result:
            response.status = 204
            return
        return repr(result)

    def doc(self):
        """Push 'line' and return result of getargspec on the final expr."""
        line = request.GET['line']
        result = ProgrammaticREPL.doc(self, line)
        if not result:
            response.status = 204
        return result

    def resources(self, val):
        pp = os.path.join(local_dir, "resources", val)
        if not os.path.exists(pp):
            response.status = 404
            return
        return open(pp).read()

