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
kd = AMRKDTree(ds)

# Print out specifics of KD Tree
print "Total volume of all bricks = %i" % kd.count_volume()
print "Total number of cells = %i" % kd.count_cells()

# Define a camera and take an volume rendering.
tf = yt.ColorTransferFunction((-30, -22))
cam = ds.camera([0.5, 0.5, 0.5], [0.2, 0.3, 0.4], 0.10, 256,
                  tf, volume=kd)
tf.add_layers(4, 0.01, col_bounds=[-27.5, -25.5], colormap='RdBu_r')
cam.snapshot("v1.png", clip_ratio=6.0)

# This rendering is okay, but lets say I'd like to improve it, and I don't want
# to spend the time rendering the high resolution data.  What we can do is
# generate a low resolution version of the AMRKDTree and pass that in to the
# camera.  We do this by specifying a maximum refinement level of 6.

kd_low_res = AMRKDTree(ds, max_level=6)
print kd_low_res.count_volume()
print kd_low_res.count_cells()

# Now we pass this in as the volume to our camera, and render the snapshot
# again.

cam.volume = kd_low_res
cam.snapshot("v4.png", clip_ratio=6.0)

# This operation was substantiall faster.  Now lets modify the low resolution
# rendering until we find something we like.

tf.clear()
tf.add_layers(4, 0.01, col_bounds=[-27.5, -25.5],
              alpha=np.ones(4, dtype='float64'), colormap='RdBu_r')
cam.snapshot("v2.png", clip_ratio=6.0)

# This looks better.  Now let's try turning on opacity.

tf.grey_opacity = True
cam.snapshot("v4.png", clip_ratio=6.0)

# That seemed to pick out som interesting structures.  Now let's bump up the
# opacity.

tf.clear()
tf.add_layers(4, 0.01, col_bounds=[-27.5, -25.5],
              alpha=10.0 * np.ones(4, dtype='float64'), colormap='RdBu_r')
cam.snapshot("v3.png", clip_ratio=6.0)

# This looks pretty good, now lets go back to the full resolution AMRKDTree

cam.volume = kd
cam.snapshot("v4.png", clip_ratio=6.0)

# This looks great!
