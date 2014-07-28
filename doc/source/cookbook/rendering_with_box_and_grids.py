import yt
import numpy as np

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")

# Create a data container (like a sphere or region) that
# represents the entire domain.
ad = ds.all_data()

# Get the minimum and maximum densities.
mi, ma = ad.quantities.extrema("density")

# Create a transfer function to map field values to colors.
# We bump up our minimum to cut out some of the background fluid
tf = yt.ColorTransferFunction((np.log10(mi)+2.0, np.log10(ma)))

# Add three guassians, evenly spaced between the min and
# max specified above with widths of 0.02 and using the
# gist_stern colormap.
tf.add_layers(3, w=0.02, colormap="gist_stern")

# Choose a center for the render.
c = [0.5, 0.5, 0.5]

# Choose a vector representing the viewing direction.
L = [0.5, 0.2, 0.7]

# Set the width of the image.
# Decreasing or increasing this value
# results in a zoom in or out.
W = 1.0

# The number of pixels along one side of the image.
# The final image will have Npixel^2 pixels.
Npixels = 512

# Create a camera object.
# This object creates the images and
# can be moved and rotated.
cam = ds.camera(c, L, W, Npixels, tf)

# Create a snapshot.
# The return value of this function could also be accepted, modified (or saved
# for later manipulation) and then put written out using write_bitmap.
# clip_ratio applies a maximum to the function, which is set to that value
# times the .std() of the array.
im = cam.snapshot("%s_volume_rendered.png" % ds, clip_ratio=8.0)

# Add the domain edges, with an alpha blending of 0.3:
nim = cam.draw_domain(im, alpha=0.3)
nim.write_png('%s_vr_domain.png' % ds)

# Add the grids, colored by the grid level with the algae colormap
nim = cam.draw_grids(im, alpha=0.3, cmap='algae')
nim.write_png('%s_vr_grids.png' % ds)

# Here we can draw the coordinate vectors on top of the image by processing
# it through the camera. Then save it out.
cam.draw_coordinate_vectors(nim)
nim.write_png("%s_vr_vectors.png" % ds)
