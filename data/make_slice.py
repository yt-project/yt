from yt.mods import *

pf = EnzoStaticOutput("DD0003/sb_L2x2_0003")

lagos.fieldInfo["VelocityMagnitude"].take_log = True

pc = PlotCollection(pf, center=[0.5,0.5,0.5])

p = pc.add_slice("Density", 2)
p.add_callback(ContourCallback("TotalEnergy", ncont=5, factor=1, take_log=True))
p.add_callback(QuiverCallback("x-velocity", "y-velocity", 32))
p.add_callback(GridBoundaryCallback())

pc.set_width(0.9, '1')

pc.add_phase_sphere(1.0, '1', ["Density", "TotalEnergy", "VelocityMagnitude"], weight=None)

pc.save("hi")
