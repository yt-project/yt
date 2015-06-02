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
volume_source.set_field(('gas','velocity_x'))
volume_source.build_default_transfer_function()
sc.render("render_x.png")

