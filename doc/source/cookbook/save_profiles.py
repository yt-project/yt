### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import yt
import matplotlib.pyplot as plt
import h5py as h5

ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Get a sphere

sp = ds.sphere(ds.domain_center, (500., "kpc"))

# Radial profile from the sphere

prof = yt.BinnedProfile1D(sp, 100, "Radiuskpc", 0.0, 500., log_space=False)
prof = yt.ProfilePlot(sp, 'radius', ['density', 'temperature'], weight_field="cell_mass")
prof.set_unit('radius', 'kpc')
prof.set_log('radius', False)
prof.set_xlim(0, 500)

# Write profiles to ASCII file

prof.write_out("%s_profile.dat" % ds, bin_style="center")

# Write profiles to HDF5 file

prof.write_out_h5("%s_profile.h5" % ds, bin_style="center")

# Now we will show how using NumPy, h5py, and Matplotlib the data in these
# files may be plotted.

# Plot density from ASCII file

# Open the text file using NumPy's "loadtxt" method. In order to get the 
# separate columns into separate NumPy arrays, it is essential to set unpack=True.

r, dens, std_dens, temp, std_temp = \
	np.loadtxt("sloshing_nomag2_hdf5_plt_cnt_0150_profile.dat", unpack=True)

fig1 = plt.figure()

ax = fig1.add_subplot(111)
ax.plot(r, dens)
ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{\rho\ (g\ cm^{-3})}$")
ax.set_title("Density vs. Radius")
fig1.savefig("%s_dens.png" % ds)

# Plot temperature from HDF5 file

# Get the file handle

f = h5py.File("%s_profile.h5" % ds, "r")

# Get the radius and temperature arrays from the file handle

r = f["/Radiuskpc-1d"].attrs["x-axis-Radiuskpc"][:]
temp = f["/Radiuskpc-1d/temperature"][:]

# Close the file handle

f.close()

fig2 = plt.figure()

ax = fig2.add_subplot(111)
ax.plot(r, temp)
ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{T\ (K)}$")
ax.set_title("temperature vs. Radius")
fig2.savefig("%s_temp.png" % ds)
