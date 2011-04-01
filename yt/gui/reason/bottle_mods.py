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
import json
import logging
from yt.utilities.logger import ytLogger as mylog
from yt.funcs import *

route_functions = {}
route_watchers = []
payloads = []

def preroute(future_route, *args, **kwargs):
    def router(func):
        route_functions[future_route] = (args, kwargs, func)
        return func
    return router

def notify_route(watcher):
    route_watchers.append(watcher)

class PayloadHandler(object):
    _shared_state = {}
    _hold = False

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        self.payloads = []

    def deliver_payloads(self):
        if self._hold: return []
        payloads = self.payloads
        self.payloads = []
        return payloads

    def add_payload(self, to_add):
        self.payloads.append(to_add)

_ph = PayloadHandler()

def append_payloads(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        reset = not _ph._hold
        _ph._hold = True
        func(self, *args, **kwargs)
        # Assume it returns a dict
        if not reset: return
        # In case it sets it manually
        _ph._hold = False
        payloads = _ph.deliver_payloads()
        return payloads
    return wrapper

class BottleDirectRouter(DirectRouter):
    # This class implements a mechanism for auto-routing an ExtDirect-callable
    # object through Bottle.  It should be used as a base class of an object,
    # and the __init__ function will need to include the keyword argument
    # 'route' for it to work.
    _route_prefix = None
    def __init__(self, *args, **kwargs):
        future_route = self.api_url
        super(BottleDirectRouter, self).__init__(*args, **kwargs)
        self.__name__ = str(self.my_name)
        route_functions[future_route] = ((), {'method':"POST"}, self)
        preroute("/resources/ext-%s-api.js" % self.api_url, method="GET")(self._myapi)
        notify_route(self)

    def _myapi(self):
        dpd = DirectProviderDefinition(self, self.api_url, ns="yt_rpc")
        source = "Ext.Direct.addProvider(%s);" % json.dumps(dpd._config())
        return source

    def __call__(self):
        #print "Hi there, I just got this request:",
        val = request.body.read()
        print val
        #import pdb;pdb.set_trace()
        rv = super(BottleDirectRouter, self).__call__(val)
        #print "With this response:", rv
        return rv

def uuid_serve_functions(pre_routed = None, open_browser=False, port=9099):
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
    for w in route_watchers:
        if not hasattr(w, "_route_prefix"):
            print "WARNING: %s has no _route_prefix attribute.  Not notifying."
            continue
            w._route_prefix = token
    print
    print
    print "============================================================================="
    print "============================================================================="
    print "Greetings, and welcome to Reason!"
    print "Your private token is %s ." % token
    print "DO NOT SHARE THIS TOKEN."
    print
    print "Please direct your browser to:"
    print
    print "     http://localhost:%s/%s/" % (port, token)
    print
    print "============================================================================="
    print
    print "If you are currently ssh'd into a remote machine, you should be able"
    print "to create a new SSH tunnel by typing or copy/pasting this text"
    print "verbatim, while waiting to see the 'ssh>' prompt after the first line."
    print
    print "~C"
    print "-L%s:localhost:%s" % (port, port)
    print
    print "and then pointing a web browser on your local machine to the above URL."
    print
    print "============================================================================="
    print "============================================================================="
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
                webbrowser.open('http://localhost:%s/%s/' % (port, token))
            thread = threading.Timer(0.5, _local_browse)
            thread.start()
        local_browse()
    try:
        import rocket
        server_name = "rocket"
        log = logging.getLogger('Rocket')
        log.setLevel(logging.INFO)
        kwargs = {'timeout': 600, 'max_threads': 1}
    except ImportError:
        server_name = "wsgiref"
        kwargs = {}
    server_type = server_names.get(server_name)
    server = server_type(host='localhost', port=port, **kwargs)
    #repl.locals['server'] = server
    mylog.info("Starting up the server.")
    run(server=server)
