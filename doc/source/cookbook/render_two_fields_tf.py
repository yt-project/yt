import numpy as np

import yt
from yt.visualization.volume_rendering.api import Scene, create_volume_source

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# create a scene and add volume sources to it

sc = Scene()

# Add density

field = "density"

vol = create_volume_source(ds, field=field)
vol.use_ghost_zones = True

tf = yt.ColorTransferFunction([-28, -25])
tf.clear()
tf.add_layers(4, 0.02, alpha=np.logspace(-3, -1, 4), colormap="winter")

vol.set_transfer_function(tf)
sc.add_source(vol)

# Add temperature

field = "temperature"

vol2 = create_volume_source(ds, field=field)
vol2.use_ghost_zones = True

tf = yt.ColorTransferFunction([4.5, 7.5])
tf.clear()
tf.add_layers(4, 0.02, alpha=np.logspace(-0.2, 0, 4), colormap="autumn")

vol2.set_transfer_function(tf)
sc.add_source(vol2)

# setup the camera

cam = sc.add_camera(ds, lens_type="perspective")
cam.resolution = (1600, 900)
cam.zoom(20.0)

# Render the image.

sc.render()

sc.save_annotated(
    "render_two_fields_tf.png",
    sigma_clip=6.0,
    tf_rect=[0.88, 0.15, 0.03, 0.8],
    render=False,
)
