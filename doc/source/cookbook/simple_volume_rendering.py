import yt
import numpy as np

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a volume rendering, which will determine data bounds, use the first
# acceptable field in the field_list, and set up a default transfer function.
#im, sc = yt.volume_render(ds, fname="%s_volume_rendered.png" % ds, sigma_clip=8.0)

# You can easily specify a different field
im, sc = yt.volume_render(ds, field=('gas','density'), fname="%s_density_volume_rendered.png" % ds, sigma_clip=8.0)

# Now increase the resolution
sc.camera.resolution = (512, 512)
im = sc.render(fname='big.png', sigma_clip=8.0)

# Now modify the transfer function
# First get the render source, in this case the entire domain, with field ('gas','density')
render_source = sc.get_source(0)
# Clear the transfer function
render_source.transfer_function.clear()
# Map a range of density values (in log space) to the Reds_r colormap
render_source.transfer_function.map_to_colormap(
        np.log10(ds.quan(5.0e-31, 'g/cm**3')),
        np.log10(ds.quan(1.0e-29, 'g/cm**3')),
        scale=30.0, colormap='RdBu_r')
im = sc.render(fname='new_tf.png', sigma_clip=None)
