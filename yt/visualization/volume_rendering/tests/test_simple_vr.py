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
import yt
from yt.testing import fake_random_ds

def test_simple_vr():
    ds = fake_random_ds(32)
    im, sc = yt.volume_render(ds, fname='test.png', clip_ratio=4.0)
    print(sc)
    return im, sc

if __name__ == "__main__":
    im, sc = test_simple_vr()
