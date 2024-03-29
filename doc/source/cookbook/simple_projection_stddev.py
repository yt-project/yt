import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create density-weighted projections of standard deviation of the velocity
# (weighted line integrals)

for normal in "xyz":
    yt.ProjectionPlot(
        ds,
        normal,
        ("gas", f"velocity_{normal}"),
        weight_field=("gas", "density"),
        moment=2,
    ).save()
