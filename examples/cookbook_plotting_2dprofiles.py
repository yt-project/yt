from yt.mods import *
pf = get_pf()

pc = PlotCollection(pf)

pc.add_phase_sphere(10.0, 'kpc', ["Density", "Temperature", "VelocityMagnitude"])

my_sphere = pf.h.sphere([0.5,0.5,0.5], 100.0/pf['kpc'])
extracted = my_sphere.extract_region(my_sphere["x-velocity"] > 1e5)
pc.add_phase_object(extracted, ["Density", "Temperature", "CellMassMsun"],
                    weight=None)

pc.save("diagrams")
