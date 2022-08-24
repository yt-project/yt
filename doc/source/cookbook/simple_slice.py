import yt

# Load the dataset.
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create gas density slices in all three axes.
yt.SlicePlot(ds, "x", ("gas", "density"), width=(800.0, "kpc")).save()
yt.SlicePlot(ds, "y", ("gas", "density"), width=(800.0, "kpc")).save()
yt.SlicePlot(ds, "z", ("gas", "density"), width=(800.0, "kpc")).save()
