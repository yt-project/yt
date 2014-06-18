import yt
import matplotlib.pyplot as plt

ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Get the first sphere
sp0 = ds.sphere(ds.domain_center, (500., "kpc"))

# Compute the bulk velocity from the cells in this sphere
bulk_vel = sp0.quantities["BulkVelocity"]()


# Get the second sphere
sp1 = ds.sphere(ds.domain_center, (500., "kpc"))

# Set the bulk velocity field parameter 
sp1.set_field_parameter("bulk_velocity", bulk_vel)

# Radial profile without correction

rp0 = yt.ProfilePlot(sp0, 'radius', 'radial_velocity')
rp0.set_unit('radius', 'kpc')
rp0.set_log('radius', False)

# Radial profile with correction for bulk velocity

rp1 = yt.ProfilePlot(sp1, 'radius', 'radial_velocity')
rp1.set_unit('radius', 'kpc')
rp1.set_log('radius', False)

#rp0.save('radial_velocity_profile_uncorrected.png')
#rp1.save('radial_velocity_profile_corrected.png')

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
