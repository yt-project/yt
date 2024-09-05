import numpy as np

import yt

# Load the dataset
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering
sc = yt.create_scene(ds, field=("gas", "density"))

# Modify the transfer function

# First get the render source, in this case the entire domain,
# with field ('gas','density')
render_source = sc.get_source()

# Clear the transfer function
render_source.transfer_function.clear()

# Map a range of density values (in log space) to the Reds_r colormap
render_source.transfer_function.map_to_colormap(
    np.log10(ds.quan(5.0e-31, "g/cm**3")),
    np.log10(ds.quan(1.0e-29, "g/cm**3")),
    scale=30.0,
    colormap="RdBu_r",
)

sc.save("new_tf.png")
