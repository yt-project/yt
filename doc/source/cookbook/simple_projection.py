from yt.mods import *

# Load the dataset.
pf = load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the density-weighted mean density.

ProjectionPlot(pf, "x", "density", weight_field = "density").save()
ProjectionPlot(pf, "y", "density", weight_field = "density").save()
ProjectionPlot(pf, "z", "density", weight_field = "density").save()
