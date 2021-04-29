import yt

# load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# create our plot
p = yt.ParticlePlot(
    ds,
    ("all", "particle_position_x"),
    ("all", "particle_position_y"),
    ("all", "particle_mass"),
    width=(0.5, 0.5),
)

# pick some appropriate units
p.set_axes_unit("kpc")
p.set_unit(("all", "particle_mass"), "Msun")

# save result
p.save()
