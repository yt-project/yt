from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
import os
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
    key_file = 'reason.key'
    fd = os.open(key_file, os.O_CREAT, 0600)
    os.close(fd)
    out_file = file(key_file, 'w')
    out_file.write("HMAC KEY: %s\n" % Pyro4.config.HMAC_KEY)
    out_file.close()
    mylog.info('See %s for HMAC key.', key_file)
    Pyro4.Daemon.serveSimple(
        {my_q: "yt.executor"},
        ns=False, verbose=True)
else:
    my_q = PyroQueueNonRoot(comm)
    my_q.run()
