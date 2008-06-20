from yt.mods import *
pf = get_pf()

fields = ["Density", "Temperature", "x-velocity"]
pc = raven.PlotCollection(pf)
pc.add_slice(fields[0],0)
pc.add_slice(fields[0],1)
pc.add_slice(fields[0],2)
pc.set_width(100.0,'kpc')
for field in fields:
   pc.switch_field(field)
   pc.save("my_data0001_100kpc")
