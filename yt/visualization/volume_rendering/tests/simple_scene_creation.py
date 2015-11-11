"""
Create a simple scene object
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import yt
from yt.testing import \
    fake_random_ds

ds = fake_random_ds(32)
sc = yt.create_scene(ds)
