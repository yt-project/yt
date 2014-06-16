from yt.mods import *
import matplotlib.pyplot as plt

ds = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Get a sphere object

sphere = ds.sphere(ds.domain_center, (500., "kpc"))

# Bin up the data from the sphere into a radial profile

rad_profile = BinnedProfile1D(sphere, 100, "Radiuskpc", 0.0, 500., log_space=False)
rad_profile.add_fields("density","temperature")

# Make plots using matplotlib

fig = plt.figure()
ax = fig.add_subplot(111)

# Plot the density as a log-log plot using the default settings
dens_plot = ax.loglog(rad_profile["Radiuskpc"], rad_profile["density"])

# Here we set the labels of the plot axes

ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{\rho\ (g\ cm^{-3})}$")

# Save the default plot

fig.savefig("density_profile_default.png" % ds)

# The "dens_plot" object is a list of plot objects. In our case we only have one,
# so we index the list by '0' to get it. 

# Plot using dashed red lines

dens_plot[0].set_linestyle("--")
dens_plot[0].set_color("red")

fig.savefig("density_profile_dashed_red.png")

# Increase the line width and add points in the shape of x's

dens_plot[0].set_linewidth(5)
dens_plot[0].set_marker("x")
dens_plot[0].set_markersize(10)

fig.savefig("density_profile_thick_with_xs.png")

# Now get rid of the line on the axes plot

ax.lines = []

# Since the rad_profile object also includes the standard deviation in each bin,
# we'll use these as errorbars. We have to make a new plot for this:

dens_err_plot = ax.errorbar(rad_profile["Radiuskpc"], rad_profile["density"],
                            yerr=rad_profile["Density_std"])
                                                        
fig.savefig("density_profile_with_errorbars.png")
