from matplotlib import pyplot

from yt.mods import *

# Load the dataset.
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere of radius 1000 kpc centered on the max density.
sphere = pf.sphere("max", (1000, "kpc"))

# Calculate and store the bulk velocity for the sphere.
bulk_velocity = sphere.quantities['BulkVelocity']()
sphere.set_field_parameter('bulk_velocity', bulk_velocity)

# Create a 1D profile object for profiles over radius
# and add a velocity profile.
profile = BinnedProfile1D(sphere, 100, "Radiuskpc", 0.1, 1000.)
profile.add_fields('VelocityMagnitude')

# Plot the average velocity magnitude.
pyplot.loglog(profile['Radiuskpc'], profile['VelocityMagnitude'],
              label='mean')
# Plot the variance of the velocity madnitude.
pyplot.loglog(profile['Radiuskpc'], profile['VelocityMagnitude_std'],
              label='std')
pyplot.xlabel('r [kpc]')
pyplot.ylabel('v [cm/s]')
pyplot.legend()

pyplot.savefig('velocity_profiles.png')
