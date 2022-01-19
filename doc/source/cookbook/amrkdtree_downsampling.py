# Using AMRKDTree Homogenized Volumes to examine large datasets
# at lower resolution.

# In this example we will show how to use the AMRKDTree to take a simulation
# with 8 levels of refinement and only use levels 0-3 to render the dataset.

# Currently this cookbook is flawed in that the data that is covered by the
# higher resolution data gets masked during the rendering.  This should be
# fixed by changing either the data source or the code in
# yt/utilities/amr_kdtree.py where data is being masked for the partitioned
# grid.  Right now the quick fix is to create a data_collection, but this
# will only work for patch based simulations that have ds.index.grids.

# We begin by loading up yt, and importing the AMRKDTree
import numpy as np

import yt
from yt.utilities.amr_kdtree.api import AMRKDTree

# Load up a dataset and define the kdtree
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
im, sc = yt.volume_render(ds, ("gas", "density"), fname="v0.png")
sc.camera.set_width(ds.arr(100, "kpc"))
render_source = sc.get_source()
kd = render_source.volume

# Print out specifics of KD Tree
print("Total volume of all bricks = %i" % kd.count_volume())
print("Total number of cells = %i" % kd.count_cells())

new_source = ds.all_data()
new_source.max_level = 3
kd_low_res = AMRKDTree(ds, data_source=new_source)
print(kd_low_res.count_volume())
print(kd_low_res.count_cells())

# Now we pass this in as the volume to our camera, and render the snapshot
# again.

render_source.set_volume(kd_low_res)
render_source.set_field(("gas", "density"))
sc.save("v1.png", sigma_clip=6.0)

# This operation was substantially faster.  Now lets modify the low resolution
# rendering until we find something we like.

tf = render_source.transfer_function
tf.clear()
tf.add_layers(
    4,
    0.01,
    col_bounds=[-27.5, -25.5],
    alpha=np.ones(4, dtype="float64"),
    colormap="RdBu_r",
)
sc.save("v2.png", sigma_clip=6.0)

# This looks better.  Now let's try turning on opacity.

tf.grey_opacity = True
sc.save("v3.png", sigma_clip=6.0)
#
## That seemed to pick out some interesting structures.  Now let's bump up the
## opacity.
#
tf.clear()
tf.add_layers(
    4,
    0.01,
    col_bounds=[-27.5, -25.5],
    alpha=10.0 * np.ones(4, dtype="float64"),
    colormap="RdBu_r",
)
tf.add_layers(
    4,
    0.01,
    col_bounds=[-27.5, -25.5],
    alpha=10.0 * np.ones(4, dtype="float64"),
    colormap="RdBu_r",
)
sc.save("v4.png", sigma_clip=6.0)
#
## This looks pretty good, now lets go back to the full resolution AMRKDTree
#
render_source.set_volume(kd)
sc.save("v5.png", sigma_clip=6.0)

# This looks great!
