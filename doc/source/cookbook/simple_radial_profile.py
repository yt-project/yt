from yt.mods import *

# Load the dataset.
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a sphere of radius 100 kpc in the center of the box.
my_sphere = pf.sphere("c", (100.0, "kpc"))

# Create a profile of the average density vs. radius.
plot = ProfilePlot(my_sphere, "Radiuskpc", "density",
                   weight_field="cell_mass")

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()
