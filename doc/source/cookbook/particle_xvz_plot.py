import yt

# load the dataset
ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

# create our plot
p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_velocity_z', ['particle_mass'])

# pick some appropriate units
p.set_unit('particle_position_x', 'Mpc')
p.set_unit('particle_velocity_z', 'km/s')
p.set_unit('particle_mass', 'Msun')

# We want to plot position and velocity in linear scale
p.set_log('particle_position_x', False)
p.set_log('particle_velocity_z', False)

# save result
p.save()
