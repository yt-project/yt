"""
Run a simple volume rendering
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import yt
import numpy as np
from yt.testing import \
    fake_random_ds

ds = fake_random_ds(32)
im, sc = yt.volume_render(ds)

angle = 2*np.pi
frames = 10
for i in range(frames):
    sc.camera.yaw(angle/frames)
    sc.render('test_rot_%04i.png' % i, clip_ratio=6.0)
