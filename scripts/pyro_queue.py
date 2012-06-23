import Pyro4
import uuid

from yt.mods import *
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    _get_comm
from yt.gui.reason.pyro_queue import \
    PyroQueueRoot, \
    PyroQueueNonRoot

comm = _get_comm(())
my_rank = comm.comm.rank

if my_rank == 0:
    my_q = PyroQueueRoot(comm)
    Pyro4.config.HMAC_KEY = uuid.uuid4().hex
    print "HMAC KEY", Pyro4.config.HMAC_KEY
    Pyro4.Daemon.serveSimple(
        {my_q: "yt.executor"},
        ns=False)
else:
    my_q = PyroQueueNonRoot(comm)
    my_q.run()
