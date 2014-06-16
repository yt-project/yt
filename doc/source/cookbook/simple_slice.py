from yt.mods import *

# Load the dataset.
ds = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices in all three axes.
SlicePlot(ds, 'x', "density", width = (800.0, 'kpc')).save()
SlicePlot(ds, 'y', "density", width = (800.0, 'kpc')).save()
SlicePlot(ds, 'z', "density", width = (800.0, 'kpc')).save()
