import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create density-weighted projections of temperature (weighted line integrals)

yt.ProjectionPlot(ds, "x", "temperature", weight_field="density").save()
yt.ProjectionPlot(ds, "y", "temperature", weight_field="density").save()
yt.ProjectionPlot(ds, "z", "temperature", weight_field="density").save()
