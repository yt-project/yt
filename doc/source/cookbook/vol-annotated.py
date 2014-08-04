#!/usr/bin/env python

import numpy as np
import pylab

import yt
import yt.visualization.volume_rendering.api as vr

ds = yt.load("maestro_subCh_plt00248")

dd = ds.all_data()

# field in the dataset we will visualize
field = ('boxlib', 'radial_velocity')

# the values we wish to highlight in the rendering.  We'll put a Gaussian
# centered on these with width sigma        
vals = [-1.e7, -5.e6, -2.5e6, 2.5e6, 5.e6, 1.e7]
sigma = 2.e5
        
mi, ma = min(vals), max(vals)

# Instantiate the ColorTransferfunction.
tf =  vr.ColorTransferFunction((mi, ma))

for v in vals:
    tf.sample_colormap(v, sigma**2, colormap="coolwarm")


# volume rendering requires periodic boundaries.  This dataset has
# solid walls.  We need to hack it for now (this will be fixed in 
# a later yt)
ds.periodicity = (True, True, True)


# Set up the camera parameters: center, looking direction, width, resolution
c = np.array([0.0, 0.0, 0.0])
L = np.array([1.0, 1.0, 1.2])
W = 1.5*ds.domain_width
N = 720

# +z is "up" for our dataset
north=[0.0,0.0,1.0]

# Create a camera object
cam = vr.Camera(c, L, W, N, transfer_function=tf, ds=ds, 
                no_ghost=False, north_vector=north,
                fields = [field], log_fields = [False])

im = cam.snapshot()

# add an axes triad 
cam.draw_coordinate_vectors(im)

# add the domain box to the image
nim = cam.draw_domain(im)

# increase the contrast -- for some reason, the enhance default
# to save_annotated doesn't do the trick 
max_val = im[:,:,:3].std() * 4.0
nim[:,:,:3] /= max_val

# we want to write the simulation time on the figure, so create a 
# figure and annotate it
f = pylab.figure()

pylab.text(0.2, 0.85, "{:.3g} s".format(float(ds.current_time.d)),
           transform=f.transFigure, color="white")

# tell the camera to use our figure
cam._render_figure = f
    
# save annotated -- this added the transfer function values, 
# and the clear_fig=False ensures it writes onto our existing figure.
cam.save_annotated("vol_annotated.png", nim, dpi=145, clear_fig=False)
