import yt

# Load the dataset.
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

p = yt.SlicePlot(ds, "x", ("gas", "density"))

# Draw a velocity vector every 16 pixels.
p.annotate_velocity(factor=16)
p.save()
