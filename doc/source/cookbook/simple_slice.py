import yt

# Load the dataset.
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices in all three axes.
yt.SlicePlot(ds, 'x', "density", width = (800.0, 'kpc')).save()
yt.SlicePlot(ds, 'y', "density", width = (800.0, 'kpc')).save()
yt.SlicePlot(ds, 'z', "density", width = (800.0, 'kpc')).save()

