import yt

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering, which will determine data bounds, use the first
# acceptable field in the field_list, and set up a default transfer function.

# This will save a file named 'data0043_Render_density.png' to disk.
im, sc = yt.volume_render(ds, field=("gas", "density"))
