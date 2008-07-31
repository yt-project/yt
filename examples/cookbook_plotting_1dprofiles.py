from yt.mods import *
pf = get_pf()

pc = PlotCollection(pf)

pc.add_profile_sphere(10.0, 'kpc', ["Density", "Temperature"])

my_sphere = pf.h.sphere([0.5,0.5,0.5], 100.0/pf['kpc'])
extracted = my_sphere.extract_region(my_sphere["x-velocity"] > 1e5)
pc.add_profile_object(extracted, ["Density", "Temperature"])

pc.save("diagrams")
