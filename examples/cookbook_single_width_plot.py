from yt.mods import *
pf = load("my_data") # Open "my_data"

pc = PlotCollection(pf)
pc.add_slice("Density",0)
pc.add_slice("Density",1)
pc.add_slice("Density",2)
pc.set_width(100.0,'kpc')
pc.save("my_data0001_100kpc")
