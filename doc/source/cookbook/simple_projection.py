from yt.mods import *

# Load the dataset.
pf = load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the density-weighted mean density.

ProjectionPlot(pf, "x", "Density", weight_field = "Density").save()
ProjectionPlot(pf, "y", "Density", weight_field = "Density").save()
ProjectionPlot(pf, "z", "Density", weight_field = "Density").save()
