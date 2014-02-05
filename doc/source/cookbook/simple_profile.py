from yt.mods import *

# Load the dataset.
pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create a 1D profile within a sphere of radius 100 kpc
# of the average temperature and average x-velocity 
# vs. density, weighted by mass.
sphere = pf.h.sphere("c", (100., "kpc"))
plot = ProfilePlot(sphere, "Density", ["Temperature", "x-velocity"],
                   weight_field="CellMassMsun")

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
plot.save()
