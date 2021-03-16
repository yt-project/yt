import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0030/data0030")

# Make a projection that is the full width of the domain,
# but only 5 Mpc in depth.  This is done by creating a
# region object with this exact geometry and providing it
# as a data_source for the projection.

# Center on the domain center
center = ds.domain_center

# First make the left and right corner of the region based
# on the full domain.
left_corner = ds.domain_left_edge
right_corner = ds.domain_right_edge

# Now adjust the size of the region along the line of sight (x axis).
depth = ds.quan(5.0, "Mpc")
left_corner[0] = center[0] - 0.5 * depth
right_corner[0] = center[0] + 0.5 * depth

# Create the region
region = ds.box(left_corner, right_corner)

# Create a density projection and supply the region we have just created.
# Only cells within the region will be included in the projection.
# Try with another data container, like a sphere or disk.
plot = yt.ProjectionPlot(
    ds, "x", ("gas", "density"), weight_field=("gas", "density"), data_source=region
)

# Save the image with the keyword.
plot.save("Thin_Slice")
