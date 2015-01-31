# This shows how to save ImageArray objects, such as those returned from 
# volume renderings, to pngs with varying backgrounds.

# First we use the simple_volume_rendering.py recipe from above to generate
# a standard volume rendering.  The only difference is that we use 
# grey_opacity=True with our TransferFunction, as the colored background 
# functionality requires images with an opacity between 0 and 1. 

# We have removed all the comments from the volume rendering recipe for 
# brevity here, but consult the recipe for more details.

import yt
import numpy as np

ds = yt.load("Enzo_64/DD0043/data0043")
ad = ds.all_data()
mi, ma = ad.quantities.extrema("density")
tf = yt.ColorTransferFunction((np.log10(mi)+1, np.log10(ma)), grey_opacity=True)
tf.add_layers(5, w=0.02, colormap="spectral")
c = [0.5, 0.5, 0.5]
L = [0.5, 0.2, 0.7]
W = 1.0
Npixels = 512
cam = ds.camera(c, L, W, Npixels, tf)
im = cam.snapshot("original.png" % ds, clip_ratio=8.0)

# Our image array can now be transformed to include different background
# colors.  By default, the background color is black.  The following
# modifications can be used on any image array.

# write_png accepts a background keyword argument that defaults to 'black'.
# Other choices include:
# black (0.,0.,0.,1.)
# white (1.,1.,1.,1.)
# None  (0.,0.,0.,0.) <-- Transparent!
# any rgba list/array: [r,g,b,a], bounded by 0..1

# We include the clip_ratio=8 keyword here to bring out more contrast between
# the background and foreground, but it is entirely optional.

im.write_png('black_bg.png', background='black', clip_ratio=8.0)
im.write_png('white_bg.png', background='white', clip_ratio=8.0)
im.write_png('green_bg.png', background=[0.,1.,0.,1.], clip_ratio=8.0)
im.write_png('transparent_bg.png', background=None, clip_ratio=8.0)
