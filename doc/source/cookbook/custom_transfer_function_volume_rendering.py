import yt
import numpy as np

# Load the dataset
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering
# NOTE: This should use yt.create_scene once that exists
im, sc = yt.volume_render(ds, field=('gas', 'density'))

# Modify the transfer function

# First get the render source, in this case the entire domain, with field ('gas','density')
render_source = sc.get_source(0)

# Clear the transfer function
render_source.transfer_function.clear()

# Map a range of density values (in log space) to the Reds_r colormap
render_source.transfer_function.map_to_colormap(
    np.log10(ds.quan(5.0e-31, 'g/cm**3')),
    np.log10(ds.quan(1.0e-29, 'g/cm**3')),
    scale=30.0, colormap='RdBu_r')

im = sc.render(fname='new_tf.png', clip_ratio=None)
