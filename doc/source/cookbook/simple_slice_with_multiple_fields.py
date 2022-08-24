import yt

# Load the dataset
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create gas density slices of several fields along the x axis simultaneously
yt.SlicePlot(
    ds,
    "x",
    [("gas", "density"), ("gas", "temperature"), ("gas", "pressure")],
    width=(800.0, "kpc"),
).save()
