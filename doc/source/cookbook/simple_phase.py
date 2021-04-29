import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere of radius 100 kpc in the center of the domain.
my_sphere = ds.sphere("c", (100.0, "kpc"))

# Create a PhasePlot object.
# Setting weight to None will calculate a sum.
# Setting weight to a field will calculate an average
# weighted by that field.
plot = yt.PhasePlot(
    my_sphere,
    ("gas", "density"),
    ("gas", "temperature"),
    ("gas", "mass"),
    weight_field=None,
)

# Set the units of mass to be in solar masses (not the default in cgs)
plot.set_unit(("gas", "mass"), "Msun")

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()
