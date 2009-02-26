
from yt.mods import *

#Substitute your parameter file location.
file_name = "/scratch/cbabbage/enzorun/DD0001/data0001"
pf = lagos.EnzoStaticOutput(file_name)
pc = raven.PlotCollection(pf)
pc.add_slice("Density", 0)
pc.add_slice("Density", 1)
pc.add_slice("Density", 2)
print pc.save("enzorun")
