import os
import tempfile
from collections import OrderedDict

import pytest

from yt.frontends.arepo.api import ArepoHDF5Dataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.answer_tests import sph_answer
from yt.utilities.answer_testing.utils import data_dir_load, requires_ds

bullet_h5 = "ArepoBullet/snapshot_150.hdf5"
tng59_h5 = "TNGHalo/halo_59.hdf5"


_tng59_bbox = [[45135.0, 51343.0], [51844.0, 56184.0], [60555.0, 63451.0]]


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestArepo:
    @pytest.mark.parametrize('ds', [bullet_h5], indirect=True)
    def test_arepo_hdf5_selection(self, ds):
        assert isinstance(ds, ArepoHDF5Dataset)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.usefixtures('hashing')
    @pytest.mark.parametrize('ds', [bullet_h5], indirect=True)
    def test_arepo_bullet(self, f, ds):
        sph = sph_answer(ds, 'snapshot_150', 26529600, f) 
        self.hashes.update({'sph_answer' : sph})

    @pytest.mark.usefixtures('hashing')
    @requires_ds(tng59_h5)
    def test_arepo_tng59(self, f):
        ds = data_dir_load(tng59_h5, kwargs = {'bounding_box': _tng59_bbox})
        sph = sph_answer(ds, 'halo_59', 10107142, f)
        self.hashes.update({'sph_answer' : sph})

    @requires_ds(tng59_h5)
    def test_index_override(self):
        # This tests that we can supply an index_filename, and that when we do, it
        # doesn't get written if our bounding_box is overwritten.
        tmpfd, tmpname = tempfile.mkstemp(suffix=".ewah")
        os.close(tmpfd)
        ds = data_dir_load(tng59_h5, kwargs = {'index_filename': tmpname,
                                               'bounding_box': _tng59_bbox})
        assert isinstance(ds, ArepoHDF5Dataset)
        ds.index
        assert len(open(tmpname, "r").read()) == 0

    @pytest.mark.parametrize('ds', [tng59_h5], indirect=True)
    def test_arepo_tng59_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()
