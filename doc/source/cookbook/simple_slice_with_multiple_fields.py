from yt.mods import *

# Load the dataset
pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices of several fields along the x axis
SlicePlot(pf, 'x', ['density','temperature','Pressure','VorticitySquared'], 
          width = (800.0, 'kpc')).save()
