from yt.mods import *

# Load the dataset.
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices in all three axes.
SlicePlot(pf, 'x', "Density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'y', "Density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'z', "Density", width = (800.0, 'kpc')).save()
