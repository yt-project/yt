import yt

# Load the dataset
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a projection and save it with the default colormap ('algae')
p = yt.ProjectionPlot(ds, "z", ("gas", "density"), width=(100, "kpc"))
p.save()

# Change the colormap to 'dusk' and save again.  We must specify
# a different filename here or it will save it over the top of
# our first projection.
p.set_cmap(field=("gas", "density"), cmap="dusk")
p.save("proj_with_dusk_cmap.png")

# Change the colormap to 'hot' and save again.
p.set_cmap(field=("gas", "density"), cmap="hot")
p.save("proj_with_hot_cmap.png")
