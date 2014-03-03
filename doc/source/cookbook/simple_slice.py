from yt.mods import *

# Load the dataset.
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices in all three axes.
SlicePlot(pf, 'x', "density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'y', "density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'z', "density", width = (800.0, 'kpc')).save()
