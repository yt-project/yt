from collections import OrderedDict
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.arepo.api import ArepoHDF5Dataset

bullet_h5 = "ArepoBullet/snapshot_150.hdf5"
tng59_h5 = "TNGHalo/halo_59.hdf5"


@requires_file(bullet_h5)
def test_arepo_hdf5():
    assert isinstance(data_dir_load(bullet_h5),
                      ArepoHDF5Dataset)


bullet_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None)
    ]
)


@requires_ds(bullet_h5)
def test_arepo_bullet():
    ds = data_dir_load(bullet_h5)
    for test in sph_answer(ds, 'snapshot_150', 40193369, 
                           bullet_fields):
        test_arepo_bullet.__name__ = test.description
        yield test


@requires_file(tng59_h5)
def test_tng_hdf5():
    assert isinstance(data_dir_load(tng59_h5),
                      ArepoHDF5Dataset)

tng59_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (('gas', 'magnetic_field_strength'), None)
    ]
)


@requires_ds(tng59_h5)
def test_arepo_tng59():
    ds = data_dir_load(tng59_h5)
    for test in sph_answer(ds, 'halo_59', 10107142,
                           tng59_fields):
        test_arepo_tng59.__name__ = test.description
        yield test

