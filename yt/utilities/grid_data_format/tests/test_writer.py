"""
Testsuite for writing yt data to GDF



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import tempfile
import shutil
import os
from yt.utilities.on_demand_imports import _h5py as h5
from yt.testing import \
    fake_random_ds, assert_equal
from yt.utilities.grid_data_format.writer import \
    write_to_gdf
from yt.frontends.gdf.data_structures import \
    GDFDataset
from yt.mods import \
    load

TEST_AUTHOR = "yt test runner"
TEST_COMMENT = "Testing write_to_gdf"


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_write_gdf():
    """Main test suite for write_gdf"""
    tmpdir = tempfile.mkdtemp()
    tmpfile = os.path.join(tmpdir, 'test_gdf.h5')

    try:
        test_ds = fake_random_ds(64)
        write_to_gdf(test_ds, tmpfile, data_author=TEST_AUTHOR,
                     data_comment=TEST_COMMENT)
        del test_ds
        assert isinstance(load(tmpfile), GDFDataset)

        h5f = h5.File(tmpfile, 'r')
        gdf = h5f['gridded_data_format'].attrs
        assert_equal(gdf['data_author'], TEST_AUTHOR)
        assert_equal(gdf['data_comment'], TEST_COMMENT)
        h5f.close()

    finally:
        shutil.rmtree(tmpdir)
