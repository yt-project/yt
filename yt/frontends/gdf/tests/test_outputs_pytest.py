"""
Title: test_gdf.py
Purpose: GDF frontend tests
Notes:
    Copyright (c) 2016, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.gdf.api import GDFDataset
from yt.testing import requires_file, units_override_check
from yt.utilities.answer_testing.answer_tests import small_patch_amr

# Test data
sedov = "sedov/sedov_tst_0004.h5"


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestGDF:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [sedov], indirect=True)
    def test_sedov_tunnel(self, a, d, w, f, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.parametrize("ds", [sedov], indirect=True)
    def test_GDFDataset(self, ds):
        assert isinstance(ds, GDFDataset)

    @requires_file(sedov)
    def test_units_override(self):
        units_override_check(sedov)
