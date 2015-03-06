# Using AMRKDTree Homogenized Volumes to examine large datasets
# at lower resolution.

# In this example we will show how to use the AMRKDTree to take a simulation
# with 8 levels of refinement and only use levels 0-3 to render the dataset.

# We begin by loading up yt, and importing the AMRKDTree
import numpy as np

import yt
from yt.utilities.amr_kdtree.api import AMRKDTree

# Load up a dataset and define the kdtree
ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
im, sc = yt.volume_render(ds, 'density', fname='v0.png')
cam = sc.get_default_camera()
cam.set_width(ds.arr(1, 'Mpc'))
render_source = sc.get_source(0)
kd=render_source.volume

# Print out specifics of KD Tree
print "Total volume of all bricks = %i" % kd.count_volume()
print "Total number of cells = %i" % kd.count_cells()

kd_low_res = AMRKDTree(ds, max_level=3)
print kd_low_res.count_volume()
print kd_low_res.count_cells()

# Now we pass this in as the volume to our camera, and render the snapshot
# again.

render_source.set_volume(kd_low_res)
render_source.set_fields('density')
sc.render("v1.png")

# This operation was substantiall faster.  Now lets modify the low resolution
# rendering until we find something we like.

tf = render_source.transfer_function
tf.clear()
tf.add_layers(4, 0.01, col_bounds=[-27.5, -25.5],
              alpha=np.ones(4, dtype='float64'), colormap='RdBu_r')
sc.render("v2.png", clip_ratio=6.0)

# This looks better.  Now let's try turning on opacity.

tf.grey_opacity = True
sc.render("v3.png", clip_ratio=6.0)
#
## That seemed to pick out som interesting structures.  Now let's bump up the
## opacity.
#
tf.clear()
tf.add_layers(4, 0.01, col_bounds=[-27.5, -25.5],
              alpha=10.0 * np.ones(4, dtype='float64'), colormap='RdBu_r')
sc.render("v5.png", clip_ratio=6.0)
#
## This looks pretty good, now lets go back to the full resolution AMRKDTree
#
render_source.set_volume(kd)
sc.render("v6.png", clip_ratio=6.0)

# This looks great!
