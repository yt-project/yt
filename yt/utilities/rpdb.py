"""
Some simple localhost-only remote pdb hacks



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import cmd
import pdb
import socket
import sys
from yt.extern.six.moves import StringIO
import traceback
import signal
from yt.extern.six.moves.xmlrpc_server import SimpleXMLRPCServer
from yt.extern.six.moves.xmlrpc_client import ServerProxy
from yt.config import ytcfg

class PdbXMLRPCServer(SimpleXMLRPCServer):
    """
    shutdown-enabled XMLRPCServer from
      http://code.activestate.com/recipes/114579/
    """

    finished=False

    def register_signal(self, signum):
        signal.signal(signum, self.signal_handler)

    def signal_handler(self, signum, frame):
        print("Caught signal", signum)
        self.shutdown()

    def shutdown(self):
        self.finished=True
        return 1

    def serve_forever(self):
        while not self.finished: self.handle_request()
        print("DONE SERVING")

def rpdb_excepthook(exc_type, exc, tb):
    traceback.print_exception(exc_type, exc, tb)
    task = ytcfg.getint("yt", "__global_parallel_rank")
    size = ytcfg.getint("yt", "__global_parallel_size")
    print("Starting RPDB server on task %s ; connect with 'yt rpdb -t %s'" \
            % (task,task))
    handler = pdb_handler(tb)
    server = PdbXMLRPCServer(("localhost", 8010+task))
    server.register_introspection_functions()
    server.register_instance(handler)
    server.register_function(server.shutdown)
    server.serve_forever()
    server.server_close()
    if size > 1:
        from mpi4py import MPI
        # This COMM_WORLD is okay.  We want to barrierize here, while waiting
        # for shutdown from the rest of the parallel group.  If you are running
        # with --rpdb it is assumed you know what you are doing and you won't
        # let this get out of hand.
        MPI.COMM_WORLD.Barrier()

class pdb_handler(object):
    def __init__(self, tb):
        self.cin = StringIO()
        sys.stdin = self.cin
        self.cout = StringIO()
        sys.stdout = self.cout
        sys.stderr = self.cout
        self.debugger = pdb.Pdb(stdin=self.cin, stdout=self.cout)
        self.debugger.reset()
        self.debugger.setup(tb.tb_frame, tb)

    def execute(self, line):
        tt = self.cout.tell()
        self.debugger.onecmd(line)
        self.cout.seek(tt)
        return self.cout.read()

class rpdb_cmd(cmd.Cmd):
    def __init__(self, proxy):
        self.proxy = proxy
        cmd.Cmd.__init__(self)
        print(self.proxy.execute("bt"))

    def default(self, line):
        print(self.proxy.execute(line))

    def do_shutdown(self, args):
        print(self.proxy.shutdown())
        return True

    def do_help(self, line):
        print(self.proxy.execute("help %s" % line))

    def postcmd(self, stop, line):
        return stop

    def postloop(self):
        try:
            self.proxy.shutdown()
        except:
            pass

__header = \
"""
You're in a remote PDB session with task %(task)s

You can run PDB commands, and when you're done, type 'shutdown' to quit.
"""

def run_rpdb(task = None):
    port = 8010
    if task is None:
        try:
            task + int(sys.argv[-1])
        except: pass
    port += task
    sp = ServerProxy("http://localhost:%s/" % port)
    try:
        pp = rpdb_cmd(sp)
    except socket.error:
        print("Connection refused.  Is the server running?")
        sys.exit(1)
    pp.cmdloop(__header % dict(task = port-8010))
