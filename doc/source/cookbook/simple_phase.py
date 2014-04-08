from yt.mods import *

# Load the dataset.
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere of radius 100 kpc in the center of the domain.
my_sphere = pf.sphere("c", (100.0, "kpc"))

# Create a PhasePlot object.
# Setting weight to None will calculate a sum.
# Setting weight to a field will calculate an average
# weighted by that field.
plot = PhasePlot(my_sphere, "density", "temperature", "cell_mass",
                 weight_field=None)

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()

