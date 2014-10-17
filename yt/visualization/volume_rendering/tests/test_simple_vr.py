"""
Test Simple Volume Rendering Scene

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.mods import *
from yt.testing import \
    fake_random_ds
from yt.visualization.volume_rendering.api import volume_render

ds = fake_random_ds(32)

im, sc = volume_render(ds, fname='test.png', clip_ratio=4.0)
