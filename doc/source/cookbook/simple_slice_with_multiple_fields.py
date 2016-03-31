import yt

# Load the dataset
ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Create density slices of several fields along the x axis
yt.SlicePlot(ds, 'x', ['density','temperature','pressure'],
             width = (800.0, 'kpc')).save()
