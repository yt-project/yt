from yt.mods import *

pf = EnzoStaticOutput("DD0003/sb_L2x2_0003")
pf.parameters["DomainLeftEdge"] = na.zeros(3, dtype='float64')
pf.parameters["DomainRightEdge"] = na.ones(3, dtype='float64')

pf.h.gridRightEdge[:,2] = 1.0
pf.h.gridDimensions[:,2] = 1.0

for g in pf.h.grids:
    g.dz = 1.0
    g.LeftEdge[2] = 0.0
    g.RightEdge[2] = 1.0
    g.ActiveDimensions[2] = 1

pc = PlotCollection(pf, center=[0.5,0.5,0.25])

pc.add_slice("Density", 2)
pc.save("hi")
