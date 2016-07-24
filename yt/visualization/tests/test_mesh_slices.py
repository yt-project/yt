"""
Tests for making unstructured mesh slices

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os
import os.path
import tempfile
import shutil
import numpy as np
import yt
from yt.testing import fake_tetrahedral_ds
from yt.testing import fake_hexahedral_ds


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_mesh_slices():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    np.random.seed(0x4d3d3d3)

    # tetrahedral ds
    ds = fake_tetrahedral_ds()

    for field in ds.field_list:
        for idir in [0, 1, 2]:
            sl = yt.SlicePlot(ds, idir, field)
            sl.annotate_mesh_lines()
            sl.save()

    # hexahedral ds
    ds = fake_hexahedral_ds()

    for field in ds.field_list:
        for idir in [0, 1, 2]:
            sl = yt.SlicePlot(ds, idir, field)
            sl.annotate_mesh_lines()
            sl.save()

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
