from yt.mods import *

max_rho = []
max_pos = []
times = []
for i in range(30):
    pf = lagos.EnzoStaticOutput("my_output%04i" % (i))
    v, c = pf.h.find_max("Density")
    max_rho.append(v)
    max_pos.append(c)
    times.append(pf["CosmologyCurrentTime"] * pf["years"])
pylab.loglog(times, max_rho)
pylab.clf()
