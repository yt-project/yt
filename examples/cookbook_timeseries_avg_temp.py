from yt.mods import *

avg_T = []
times = []
for i in range(30):
    pf = lagos.EnzoStaticOutput("my_output%04i" % (i))
    v, c = pf.h.find_max("Density")
    sp = pf.h.sphere(c, 10.0/pf['kpc'])
    avg_T.append(sp.quantities["WeightedAverageQuantity"]\
        ("Temperature", "CellMassMsun"))
    times.append(pf["CosmologyCurrentTime"] * pf["years"])
pylab.loglog(times, avg_T)
