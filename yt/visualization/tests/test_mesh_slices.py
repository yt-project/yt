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
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    GenericImageTest


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, field, test_prefix, decimals=12):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, 'z', field)
        sl.set_log('all', False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = "slice_{}".format(test_prefix)
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    return test

tri2 = "SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e"

@requires_ds(tri2)
def test_tri2():
    ds = data_dir_load(tri2, kwargs={'step':-1})
    for field in ds.field_list:
        yield compare(ds, field, "answers_tri2_%s_%s" % (field[0], field[1]))

multi_region = "MultiRegion/two_region_example_out.e"

@requires_ds(multi_region)
def test_multi_region():
    ds = data_dir_load(multi_region, kwargs={'step':-1})
    for field in ds.field_list:
        yield compare(ds, field, "answers_multi_region_%s_%s" % (field[0], field[1]))

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
