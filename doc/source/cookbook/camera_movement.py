import yt
import numpy as np

# Follow the simple_volume_rendering cookbook for the first part of this.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")  # load data
im, sc = yt.volume_render(ds)
cam = sc.camera
cam.resolution = (512, 512)
cam.set_width(ds.domain_width/20.0)

# Find the maximum density location, store it in max_c
v, max_c = ds.find_max('density')

frame = 0
# Move to the maximum density location over 5 frames
for _ in cam.move_to(max_c, 5):
    sc.render('camera_movement_%04i.png' % frame, clip_ratio=8.0)
    frame += 1

# Zoom in by a factor of 10 over 5 frames
for _ in cam.zoomin(10.0, 5):
    sc.render('camera_movement_%04i.png' % frame, clip_ratio=8.0)
    frame += 1

# Do a rotation over 5 frames
for _ in cam.rotation(np.pi, 5):
    sc.render('camera_movement_%04i.png' % frame, clip_ratio=8.0)
    frame += 1
