import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create projections of the gas density (non-weighted line integrals).

yt.ProjectionPlot(ds, "x", ("gas", "density")).save()
yt.ProjectionPlot(ds, "y", ("gas", "density")).save()
yt.ProjectionPlot(ds, "z", ("gas", "density")).save()
