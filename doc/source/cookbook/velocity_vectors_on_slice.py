from yt.mods import *

# Load the dataset.
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

p = SlicePlot(pf, "x", "density")
# Draw a velocity vector every 16 pixels.
p.annotate_velocity(factor = 16)
p.save()
