from yt.mods import *

min_output_number = 0
max_output_number = 30
skip = 1

rho_min = 1e-30
rho_max = 1e-20

frame_template = "frames/frame_%04i.png"
basename_template = "galaxy%04i.dir/galaxy%04i"

for i in range(min_output_number, max_output_number+1, skip):
    basename = basename_template % (i,i)
    pf = lagos.EnzoStaticOutput(basename)
    pc = raven.PlotCollection(pf, center=[0.5,0.5,0.5])
    pc.add_projection("Density",0)
    pc.set_zlim(rho_min, rho_max)
    # Override the name
    pc.save(frame_template % (i), override=True)
