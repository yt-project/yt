import yt

# load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# create our plot
p = yt.ParticlePlot(
    ds, ("all", "particle_position_x"), ("all", "particle_position_y"), color="b"
)

# zoom in a little bit
p.set_width(500, "kpc")

# save result
p.save()
