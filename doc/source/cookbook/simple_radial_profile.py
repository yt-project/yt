import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere of radius 100 kpc in the center of the box.
my_sphere = ds.sphere("c", (100.0, "kpc"))

# Create a profile of the average density vs. radius.
plot = yt.ProfilePlot(
    my_sphere,
    ("index", "radius"),
    ("gas", "density"),
    weight_field=("gas", "mass"),
)

# Change the units of the radius into kpc (and not the default in cgs)
plot.set_unit(("index", "radius"), "kpc")

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()
