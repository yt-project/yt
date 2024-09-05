import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Make a density projection centered on the 'm'aximum density location
# with a width of 10 Mpc..
p = yt.ProjectionPlot(ds, "y", ("gas", "density"), center="m", width=(10, "Mpc"))

# Modify the projection
# The argument specifies the region along the line of sight
# for which particles will be gathered.
# 1.0 signifies the entire domain in the line of sight
# p.annotate_particles(1.0)
# but in this case we only go 10 Mpc in depth
p.annotate_particles((10, "Mpc"))

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
p.save()
