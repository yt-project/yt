from yt.mods import *

# Load the dataset.
ds = load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the density-weighted mean density.

ProjectionPlot(ds, "x", "density", weight_field = "density").save()
ProjectionPlot(ds, "y", "density", weight_field = "density").save()
ProjectionPlot(ds, "z", "density", weight_field = "density").save()
