import yt
import numpy as np

# Follow the simple_volume_rendering cookbook for the first part of this.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")  # load data
ad = ds.all_data()
mi, ma = ad.quantities.extrema("density")

# Set up transfer function
tf = yt.ColorTransferFunction((np.log10(mi), np.log10(ma)))
tf.add_layers(6, w=0.05)

# Set up camera paramters
c = [0.5, 0.5, 0.5]  # Center
L = [1, 1, 1]  # Normal Vector
W = 1.0  # Width
Nvec = 512  # Pixels on a side

# Specify a north vector, which helps with rotations.
north_vector = [0., 0., 1.]

# Find the maximum density location, store it in max_c
v, max_c = ds.find_max('density')

# Initialize the Camera
cam = ds.camera(c, L, W, (Nvec, Nvec), tf, north_vector=north_vector)
frame = 0

# Do a rotation over 5 frames
for i, snapshot in enumerate(cam.rotation(np.pi, 5, clip_ratio=8.0)):
    snapshot.write_png('camera_movement_%04i.png' % frame)
    frame += 1

# Move to the maximum density location over 5 frames
for i, snapshot in enumerate(cam.move_to(max_c, 5, clip_ratio=8.0)):
    snapshot.write_png('camera_movement_%04i.png' % frame)
    frame += 1

# Zoom in by a factor of 10 over 5 frames
for i, snapshot in enumerate(cam.zoomin(10.0, 5, clip_ratio=8.0)):
    snapshot.write_png('camera_movement_%04i.png' % frame)
    frame += 1
