from yt.mods import *
import matplotlib.pyplot as plt

ds = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Get the first sphere

sphere0 = ds.sphere(ds.domain_center, (500., "kpc"))

# Compute the bulk velocity from the cells in this sphere

bulk_vel = sphere0.quantities["BulkVelocity"]()

# Get the second sphere

sphere1 = ds.sphere(ds.domain_center, (500., "kpc"))

# Set the bulk velocity field parameter 
sphere1.set_field_parameter("bulk_velocity", bulk_vel)

# Radial profile without correction

rad_profile0 = BinnedProfile1D(sphere0, 100, "Radiuskpc", 0.0, 500., log_space=False)
rad_profile0.add_fields("RadialVelocity")

# Radial profile with correction for bulk velocity

rad_profile1 = BinnedProfile1D(sphere1, 100, "Radiuskpc", 0.0, 500., log_space=False)
rad_profile1.add_fields("RadialVelocity")

# Make a plot using matplotlib

fig = plt.figure()
ax = fig.add_subplot(111)

# Here we scale the velocities by 1.0e5 to get into km/s
ax.plot(rad_profile0["Radiuskpc"], rad_profile0["RadialVelocity"]/1.0e5,
		rad_profile1["Radiuskpc"], rad_profile1["RadialVelocity"]/1.0e5)

ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{v_r\ (km/s)}$")
ax.legend(["Without Correction", "With Correction"])

fig.savefig("%s_profiles.png" % ds)
