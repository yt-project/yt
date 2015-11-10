"""
Test for off_axis_projection and write_projection



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.testing import \
    fake_random_ds
from yt.visualization.image_writer import write_projection
from yt.visualization.volume_rendering.api import off_axis_projection


test_ds = fake_random_ds(64)
c = test_ds.domain_center 
norm = [0.5, 0.5, 0.5]
W = test_ds.domain_width.max() 
N = 64
field = ('gas', "density")
oap_args = [test_ds, c, norm, W, N, field]

image, sc = off_axis_projection(*oap_args)
write_projection(image, 'oap.png')

