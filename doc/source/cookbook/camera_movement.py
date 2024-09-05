import numpy as np

import yt

ds = yt.load("MOOSE_sample_data/out.e-s010")
sc = yt.create_scene(ds)
cam = sc.camera

# save an image at the starting position
frame = 0
sc.save("camera_movement_%04i.png" % frame)
frame += 1

# Zoom out by a factor of 2 over 5 frames
for _ in cam.iter_zoom(0.5, 5):
    sc.save("camera_movement_%04i.png" % frame)
    frame += 1

# Move to the position [-10.0, 10.0, -10.0] over 5 frames
pos = ds.arr([-10.0, 10.0, -10.0], "code_length")
for _ in cam.iter_move(pos, 5):
    sc.save("camera_movement_%04i.png" % frame)
    frame += 1

# Rotate by 180 degrees over 5 frames
for _ in cam.iter_rotate(np.pi, 5):
    sc.save("camera_movement_%04i.png" % frame)
    frame += 1
