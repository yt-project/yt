from yt.mods import *
pf = get_pf()

widths = [1000.0, 100.0, 10.0, 1.0]
units = ['mpc','kpc','pc','au']
my_pairs = [ (w,u) for u in units for w in widths ]

pc = raven.PlotCollection(pf)
pc.add_slice("Density",0)
pc.add_slice("Density",1)
pc.add_slice("Density",2)
for w, u in my_pairs:
    pc.set_width(w,u)
    pc.save("my_data0001_%05i%s" % (w, u))

