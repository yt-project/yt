import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of temperature (with different methods)


for method in ["integrate", "min", "max"]:
    proj = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), method=method)
    proj.save(f"projection_method_{method}.png")
