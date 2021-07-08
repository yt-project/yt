import os
import tempfile
from collections import OrderedDict

from yt.frontends.arepo.api import ArepoHDF5Dataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds, sph_answer

bullet_h5 = "ArepoBullet/snapshot_150.hdf5"
tng59_h5 = "TNGHalo/halo_59.hdf5"
_tng59_bbox = [[45135.0, 51343.0], [51844.0, 56184.0], [60555.0, 63451.0]]


@requires_file(bullet_h5)
def test_arepo_hdf5_selection():
    ds = data_dir_load(bullet_h5)
    assert isinstance(ds, ArepoHDF5Dataset)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


bullet_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
    ]
)


@requires_ds(bullet_h5)
def test_arepo_bullet():
    ds = data_dir_load(bullet_h5)
    for test in sph_answer(ds, "snapshot_150", 26529600, bullet_fields):
        test_arepo_bullet.__name__ = test.description
        yield test


tng59_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "H_number_density"), None),
        (("gas", "H_p0_number_density"), None),
        (("gas", "H_p1_number_density"), None),
        (("gas", "El_number_density"), None),
        (("gas", "C_number_density"), None),
        (("gas", "velocity_magnitude"), None),
        (("gas", "magnetic_field_strength"), None),
    ]
)


@requires_ds(tng59_h5)
def test_arepo_tng59():
    ds = data_dir_load(tng59_h5, kwargs={"bounding_box": _tng59_bbox})
    for test in sph_answer(ds, "halo_59", 10107142, tng59_fields):
        test_arepo_tng59.__name__ = test.description
        yield test


@requires_ds(tng59_h5)
def test_index_override():
    # This tests that we can supply an index_filename, and that when we do, it
    # doesn't get written if our bounding_box is overwritten.
    tmpfd, tmpname = tempfile.mkstemp(suffix=".ewah")
    os.close(tmpfd)
    ds = data_dir_load(
        tng59_h5, kwargs={"index_filename": tmpname, "bounding_box": _tng59_bbox}
    )
    assert isinstance(ds, ArepoHDF5Dataset)
    ds.index
    assert len(open(tmpname).read()) == 0


@requires_file(tng59_h5)
def test_arepo_tng59_selection():
    ds = data_dir_load(tng59_h5, kwargs={"bounding_box": _tng59_bbox})
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()
