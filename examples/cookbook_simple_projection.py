from yt.mods import *

pf = load("my_data") # Open "my_data"

pc = PlotCollection(pf)
pc.add_slice("Density", 0)
pc.add_slice("Density", 1)
pc.add_slice("Density", 2)

print pc.save("enzorun")
