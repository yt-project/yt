"""
Modifications and extensions to Bottle, to make it slightly more useful for
yt's purposes

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
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

from yt.utilities.bottle import \
    server_names, debug, route, run, request, ServerAdapter, response
import uuid
from extdirect_router import DirectRouter, DirectProviderDefinition
import json
import logging, threading
from yt.utilities.logger import ytLogger as mylog
from yt.funcs import *
import sys
import urllib, urllib2

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

class BinaryDelivery(object):
    delivered = False
    payload = ""
    def __init__(self, payload):
        self.payload = payload

    def get(self):
        # We set our 
        p = self.payload
        if p == "":
            response.status = 404
            return
        self.payload = ""
        return p

class PayloadHandler(object):
    _shared_state = {}
    payloads = None
    binary_payloads = None
    recorded_payloads = None
    multicast_ids = None
    multicast_payloads = None
    lock = None
    record = False
    event = None
    count = 0
    debug = False
    _prefix = ""

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        if self.payloads is None: self.payloads = []
        if self.lock is None: self.lock = threading.Lock()
        if self.recorded_payloads is None: self.recorded_payloads = []
        if self.event is None: self.event = threading.Event()
        if self.multicast_payloads is None: self.multicast_payloads = {}
        if self.multicast_ids is None: self.multicast_ids = {}
        if self.binary_payloads is None: self.binary_payloads = []

    def deliver_payloads(self):
        with self.lock:
            payloads = self.payloads
            if self.record:
                self.recorded_payloads += self.payloads
            if self.debug:
                sys.__stderr__.write("**** Delivering %s payloads\n" % (len(payloads)))
                for p in payloads:
                    sys.__stderr__.write("****    %s\n" % p['type'])
            self.payloads = []
            self.event.clear()
            try:
                self.deliver_multicast()
            except Exception as exc:
                sys.__stderr__.write("%s" % exc)
        return payloads

    def add_payload(self, to_add):
        with self.lock:
            if "binary" in to_add:
                self._add_binary_payload(to_add)
            self.payloads.append(to_add)
            # Does this next part need to be in the lock?
            if to_add.get("widget_id", None) in self.multicast_ids:
                self.multicast_payloads[to_add["widget_id"]] = to_add
            self.count += 1
            self.event.set()
            if self.debug:
                sys.__stderr__.write("**** Adding payload of type %s\n" % (to_add['type']))

    def _add_binary_payload(self, bp):  
        # This shouldn't be called by anybody other than add_payload.
        bdata = bp.pop(bp['binary']) # Get the binary data
        bpserver = BinaryDelivery(bdata)
        self.binary_payloads.append(bpserver)
        uu = uuid.uuid4().hex
        bp['binary'] = uu
        route("%s/%s" % (self._prefix, uu))(bpserver.get)
        sys.__stderr__.write("**** Adding binary payload to %s\n" % (uu))

    def replay_payloads(self):
        return self.recorded_payloads

    def widget_payload(self, widget, data):
        data['type'] = 'widget_payload'
        data['widget_id'] = widget._ext_widget_id
        self.add_payload(data)

    def deliver_multicast(self):
        for widget_id in self.multicast_payloads:
            if widget_id not in self.multicast_payloads: continue
            server_id, session_token = self.multicast_ids[widget_id]
            # Now we execute a post to the correct location
            data = urllib.urlencode({
                'payload_session_id': server_id,
                'payload_session_token': session_token,
                'payload_data': self.multicast_payloads[widget_id],
                'payload_metadata': {}
            })
            urllib2.urlopen("http://localhost:8080/UpdatePayload", data = data)

class YTRocketServer(ServerAdapter):
    server_info = {} # Hack to get back at instance vars
    def run(self, handler):
        from yt.utilities.rocket import Rocket
        server = Rocket((self.host, self.port), 'wsgi', { 'wsgi_app' : handler })
        self.server_info[id(self)] = server
        server.start()

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
        route_functions[future_route] = ((), {'method':("POST", "GET")}, self)
        preroute("/resources/ext-%s-api.js" % self.api_url, method="GET")(self._myapi)
        notify_route(self)

    def _myapi(self):
        dpd = DirectProviderDefinition(self, self.api_url, ns="yt_rpc")
        source = "Ext.Direct.addProvider(%s);" % json.dumps(dpd._config())
        response.headers['Content-Type'] = "text/javascript"
        return source

    def __call__(self):
        #print "Hi there, I just got this request:",
        val = request.body.read()
        #print val
        #import pdb;pdb.set_trace()
        rv = super(BottleDirectRouter, self).__call__(val)
        #print "With this response:", rv
        return rv

def uuid_serve_functions(pre_routed = None, open_browser=False, port=9099,
                         repl = None, token = None):
    if pre_routed == None: pre_routed = route_functions
    debug(mode=True)
    if token is None: token = uuid.uuid1()
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
    repl._global_token = token
    repl.activate()
    repl.execution_thread.wait()
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
        import yt.utilities.rocket as rocket
        server_type = YTRocketServer
        log = logging.getLogger('Rocket')
        log.setLevel(logging.WARNING)
        kwargs = {'timeout': 600, 'max_threads': 2}
        if repl is not None:
            repl.server = YTRocketServer.server_info
    except ImportError:
        server_type = server_names.get("wsgiref")
        kwargs = {}
    server = server_type(host='localhost', port=port, **kwargs)
    mylog.info("Starting up the server.")
    run(server=server)

class MethodLock(object):
    _shared_state = {}
    locks = None

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        if self.locks is None: self.locks = {}

    def __call__(self, func):
        if str(func) not in self.locks:
            self.locks[str(func)] = threading.Lock()
        @wraps(func)
        def locker(*args, **kwargs):
            print "Acquiring lock on %s" % (str(func))
            with self.locks[str(func)]:
                rv = func(*args, **kwargs)
            print "Regained lock on %s" % (str(func))
            return rv
        return locker

lockit = MethodLock()

