# This shows how to save ImageArray objects, such as those returned from 
# volume renderings, to pngs with varying backgrounds.

import yt
import numpy as np

# Lets make a fake "rendering" that has 4 channels and looks like a linear
# gradient from the bottom to top.

im = np.zeros([64,128,4])
for i in xrange(im.shape[0]):
    for k in xrange(im.shape[2]):
        im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])
im_arr = yt.ImageArray(im)

# in this case you would have gotten im_arr from something like:
# im_arr = cam.snapshot() 

# To save it with the default settings, we can just use write_png, where it 
# rescales the image and uses a black background.

im_arr.write_png('standard.png')
 
# write_png accepts a background keyword argument that defaults to 'black'.
# Other choices include:
# black (0.,0.,0.,1.)
# white (1.,1.,1.,1.)
# None  (0.,0.,0.,0.) <-- Transparent!
# any rgba list/array: [r,g,b,a], bounded by 0..1

im_arr.write_png('black_bg.png', background='black')
im_arr.write_png('white_bg.png', background='white')
im_arr.write_png('green_bg.png', background=[0.,1.,0.,1.])
im_arr.write_png('transparent_bg.png', background=None)
