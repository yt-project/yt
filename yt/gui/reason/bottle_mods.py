"""
Modifications and extensions to Bottle, to make it slightly more useful for
yt's purposes

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

from .bottle import server_names, debug, route, run, request
import uuid
from extdirect_router import DirectRouter, DirectProviderDefinition

route_functions = {}

def preroute(future_route, *args, **kwargs):
    def router(func):
        route_functions[future_route] = (args, kwargs, func)
        return func
    return router

class BottleDirectRouter(DirectRouter):
    # This class implements a mechanism for auto-routing an ExtDirect-callable
    # object through Bottle.  It should be used as a base class of an object,
    # and the __init__ function will need to include the keyword argument
    # 'route' for it to work.
    def __init__(self, *args, **kwargs):
        future_route = kwargs.pop("route")
        super(BottleDirectRouter, self).__init__(*args, **kwargs)
        self.__name__ = ""
        route_functions[future_route] = ((), {'method':"POST"}, self)

    def _myapi(self):
        dpd = DirectProviderDefinition(self, self.api_url)
        return dpd.render()

    def __call__(self):
        val = request.body.read()
        rv = super(BottleDirectRouter, self).__call__(val)
        return rv

def uuid_serve_functions(pre_routed = None, open_browser=False):
    if pre_routed == None: pre_routed = route_functions
    debug(mode=True)
    token = uuid.uuid1()
    for r in pre_routed:
        args, kwargs, f = pre_routed[r]
        if r[0] == "/": r = r[1:]
        rp = "/%s/%s" % (token, r)
        func_name = getattr(f, 'func_name', str(f))
        print "Routing from %s => %s" % (rp, func_name)
        route(rp, *args, **kwargs)(f)
    print "Greetings! Your private token is %s ." % token
    print
    print "Please direct your browser to:"
    print
    print "     http://localhost:8080/%s/" % token
    print
    print
    if open_browser:
        # We do some fancy footwork so that we can open the browser while the
        # server starts up.  I got this from some recipe whose URL escapes me.
        # Thank you, to whoever wrote it!
        def local_browse():
            """Start a browser after waiting for half a second."""
            import webbrowser, threading
            def _local_browse():
                webbrowser.open('http://localhost:%s/%s/' % (8080, token))
            thread = threading.Timer(0.5, _local_browse)
            thread.start()
        local_browse()
    # Right now we only really support the built-in wsgiref, but this may
    # change if we start using Rocket.
    server_type = server_names.get("wsgiref")
    server = server_type(host='localhost', port=8080)
    #repl.locals['server'] = server
    run(server=server)
