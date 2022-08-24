import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Make a density projection.
p = yt.ProjectionPlot(ds, "y", ("gas", "density"))

# Modify the projection
# The argument specifies the region along the line of sight
# for which particles will be gathered.
# 1.0 signifies the entire domain in the line of sight.
p.annotate_grids()

# Save the image.
# Optionally, give a string as an argument
# to name files with a keyword.
p.save()
