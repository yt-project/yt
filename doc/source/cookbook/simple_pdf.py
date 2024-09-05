import yt

# Load the dataset.
ds = yt.load("GalaxyClusterMerger/fiducial_1to3_b0.273d_hdf5_plt_cnt_0175")

# Create a data object that represents the whole box.
ad = ds.all_data()

# This is identical to the simple phase plot, except we supply
# the fractional=True keyword to divide the profile data by the sum.
plot = yt.PhasePlot(
    ad,
    ("gas", "density"),
    ("gas", "temperature"),
    ("gas", "mass"),
    weight_field=None,
    fractional=True,
)

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()
