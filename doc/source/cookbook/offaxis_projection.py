import numpy as np

import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Choose a center for the render.
c = [0.5, 0.5, 0.5]

# Our image plane will be normal to some vector.  For things like collapsing
# objects, you could set it the way you would a cutting plane -- but for this
# dataset, we'll just choose an off-axis value at random.  This gets normalized
# automatically.
L = [1.0, 0.0, 0.0]

# Our "width" is the width of the image plane as well as the depth.
# The first element is the left to right width, the second is the
# top-bottom width, and the last element is the back-to-front width
# (all in code units)
W = [0.04, 0.04, 0.4]

# The number of pixels along one side of the image.
# The final image will have Npixel^2 pixels.
Npixels = 512

# Create the off axis projection.
# Setting no_ghost to False speeds up the process, but makes a
# slightly lower quality image.
image = yt.off_axis_projection(ds, c, L, W, Npixels, ("gas", "density"), no_ghost=False)

# Write out the final image and give it a name
# relating to what our dataset is called.
# We save the log of the values so that the colors do not span
# many orders of magnitude.  Try it without and see what happens.
yt.write_image(np.log10(image), f"{ds}_offaxis_projection.png")
