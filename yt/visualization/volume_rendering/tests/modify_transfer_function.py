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

volume_source = sc.get_source(0)
tf = volume_source.transfer_function
tf.clear()
tf.grey_opacity=True
tf.add_layers(3, colormap='RdBu')
sc.render("new_tf.png")

