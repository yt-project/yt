import yt
import numpy as np

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# We start by building a transfer function, and initializing a camera.

tf = yt.ColorTransferFunction((-30, -22))
cam = ds.camera([0.5, 0.5, 0.5], [0.2, 0.3, 0.4], 0.10, 256, tf)

# Now let's add some isocontours, and take a snapshot.

tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5], colormap = 'RdBu_r')
cam.snapshot("v1.png", clip_ratio=6.0)

# In this case, the default alphas used (np.logspace(-3,0,Nbins)) does not
# accentuate the outer regions of the galaxy. Let's start by bringing up the
# alpha values for each contour to go between 0.1 and 1.0

tf.clear()
tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=np.logspace(0,0,4), colormap = 'RdBu_r')
cam.snapshot("v2.png", clip_ratio=6.0)

# Now let's set the grey_opacity to True.  This should make the inner portions
# start to be obcured

tf.grey_opacity = True
cam.snapshot("v3.png", clip_ratio=6.0)

# That looks pretty good, but let's start bumping up the opacity.

tf.clear()
tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=10.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
cam.snapshot("v4.png", clip_ratio=6.0)

# Let's bump up again to see if we can obscure the inner contour.

tf.clear()
tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=30.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
cam.snapshot("v5.png", clip_ratio=6.0)

# Now we are losing sight of everything.  Let's see if we can obscure the next
# layer

tf.clear()
tf.add_layers(4, 0.01, col_bounds = [-27.5,-25.5],
        alpha=100.0*np.ones(4,dtype='float64'), colormap = 'RdBu_r')
cam.snapshot("v6.png", clip_ratio=6.0)

# That is very opaque!  Now lets go back and see what it would look like with
# grey_opacity = False

tf.grey_opacity=False
cam.snapshot("v7.png", clip_ratio=6.0)

# That looks pretty different, but the main thing is that you can see that the
# inner contours are somewhat visible again.  
