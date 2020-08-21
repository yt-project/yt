import os
import tempfile

import numpy as np
import pytest

from yt.frontends.arepo.api import ArepoHDF5Dataset
from yt.testing import ParticleSelectionComparison
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    pixelized_projection_values,
    sph_validation,
)
from yt.utilities.answer_testing.utils import data_dir_load, requires_ds

# Test data
bullet_h5 = "ArepoBullet/snapshot_150.hdf5"
tng59_h5 = "TNGHalo/halo_59.hdf5"

_tng59_bbox = [[45135.0, 51343.0], [51844.0, 56184.0], [60555.0, 63451.0]]
tng_kwargs = {"kwargs": {"bounding_box": _tng59_bbox}}

# Test parameters
val_params = [
    (bullet_h5, "snapshot_150", 26529600),
    ([tng59_h5, tng_kwargs], "halo_59", 10107142),
]

a_list = [0, 1, 2]

d_bullet = [None, ("sphere", ("c", (0.1, "unitary")))]
d_tng = [None, ("sphere", ("c", (0.5, "unitary")))]

f_bullet = [
    ("gas", "density"),
    ("gas", "temperature"),
    ("gas", "temperature"),
    ("gas", "velocity_magnitude"),
]
f_tng = [
    ("gas", "density"),
    ("gas", "temperature"),
    ("gas", "temperature"),
    ("gas", "H_number_density"),
    ("gas", "H_p0_number_density"),
    ("gas", "H_p1_number_density"),
    ("gas", "El_number_density"),
    ("gas", "C_number_density"),
    ("gas", "velocity_magnitude"),
    ("gas", "magnetic_field_strength"),
]

w_bullet = [
    None,
    None,
    ("gas", "density"),
    None,
]
w_tng = [
    None,
    None,
    ("gas", "density"),
    None,
    None,
    None,
    None,
    None,
    None,
    None,
]

pair_list = [
    [bullet_h5, f_bullet, w_bullet, d_bullet],
    [tng59_h5, f_tng, w_tng, d_tng],
]


ppv_pairs = [
    (i[0], f[0], f[1], d) for i in pair_list for f in zip(i[1], i[2]) for d in i[3]
]

fv_pairs = [(i[0], f, d) for i in pair_list for f in i[1] for d in i[3]]


@pytest.mark.answer_test
class TestArepo:
    answer_file = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, w, d", ppv_pairs, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    def test_arepo_ppv(self, a, d, w, f, ds):
        particle_type = f[0] in ds.particle_types
        if not particle_type:
            ppv = pixelized_projection_values(ds, a, f, w, d)
            self.hashes.update({"pixelized_projection_values": ppv})
        # So we have something to save for this test in the answer file
        else:
            self.hashes.update({"pixelized_projection_values": np.array(-1)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d", fv_pairs, indirect=True)
    def test_arepo_fv(self, d, f, ds):
        particle_type = f[0] in ds.particle_types
        fv = field_values(ds, f, d, particle_type=particle_type)
        self.hashes.update({"field_values": fv})

    @pytest.mark.parametrize("ds, ds_repr, N", val_params, indirect=True)
    def test_arepo_validation(self, ds, ds_repr, N):
        sph_validation(ds, ds_repr, N)

    @pytest.mark.parametrize("ds", [bullet_h5], indirect=True)
    def test_arepo_hdf5_selection(self, ds):
        assert isinstance(ds, ArepoHDF5Dataset)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @requires_ds(tng59_h5)
    def test_index_override(self):
        # This tests that we can supply an index_filename, and that when we do, it
        # doesn't get written if our bounding_box is overwritten.
        tmpfd, tmpname = tempfile.mkstemp(suffix=".ewah")
        os.close(tmpfd)
        ds = data_dir_load(
            tng59_h5, kwargs={"index_filename": tmpname, "bounding_box": _tng59_bbox}
        )
        assert isinstance(ds, ArepoHDF5Dataset)
        ds.index
        assert len(open(tmpname, "r").read()) == 0

    @pytest.mark.parametrize("ds", [tng59_h5], indirect=True)
    def test_arepo_tng59_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()
