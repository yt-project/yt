from yt.mods import *

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
