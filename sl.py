from yt.mods import *
pf = load("tests/DD0010/moving7_0010")
#sl = pf.slice(0, 0.5)
#print sl["Density"]
sl = SlicePlot(pf, "x", "Density")
sl.save()
