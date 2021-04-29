import yt

# load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# create our plot
p = yt.ParticlePlot(
    ds,
    ("all", "particle_position_x"),
    ("all", "particle_velocity_z"),
    [("all", "particle_mass")],
)

# pick some appropriate units
p.set_unit(("all", "particle_position_x"), "Mpc")
p.set_unit(("all", "particle_velocity_z"), "km/s")
p.set_unit(("all", "particle_mass"), "Msun")

# save result
p.save()
