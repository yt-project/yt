
from yt.mods import *
pf = lagos.EnzoStaticOutput("/scratch/cbabbage/enzorun/DD0001/data0001") #Substitute your parameter file location.
pc = raven.PlotCollection(pf)
pc.add_slice("Density", 0)
pc.add_slice("Density", 1)
pc.add_slice("Density", 2)
print pc.save("AwesomeSimulation")
