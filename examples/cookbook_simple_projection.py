from yt.mods import *

pf = load("my_data") # Open "my_data"

pc = PlotCollection(pf)
pc.add_projection("Density", 0)
pc.add_projection("Density", 1)
pc.add_projection("Density", 2)

print pc.save("enzorun")
