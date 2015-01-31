import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the density (non-weighted line integrals).

yt.ProjectionPlot(ds, "x", "density").save()
yt.ProjectionPlot(ds, "y", "density").save()
yt.ProjectionPlot(ds, "z", "density").save()
