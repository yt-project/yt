import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the density-weighted mean density.

yt.ProjectionPlot(ds, "x", "density", weight_field = "density").save()
yt.ProjectionPlot(ds, "y", "density", weight_field = "density").save()
yt.ProjectionPlot(ds, "z", "density", weight_field = "density").save()

