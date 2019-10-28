"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
gc = "GaussianCloud/data.0077.3d.hdf5"
tb = "TurbBoxLowRes/data.0005.3d.hdf5"
iso = "IsothermalSphere/data.0000.3d.hdf5"
zp = "ZeldovichPancake/plt32.2d.hdf5"
kho = "KelvinHelmholtz/data.0004.hdf5"




def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ in ['test_gc', 'test_tb', 'test_iso', 'test_kho']:
        axes = [0, 1, 2]
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        fields = ("density", "velocity_magnitude", "magnetic_field_x")
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=['density', 'velocity_magnitude',
            'magnetic_field_x'])
    if metafunc.function.__name__ == 'test_zp':
        axes = [0, 1, 2]
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "rhs"]
        fields = ("rhs", "phi")
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'rhs'])
        metafunc.parametrize('f', fields, ids=['rhs', 'phi'])


@pytest.fixture(scope='class')
def ds_gc():
    ds = utils.data_dir_load(gc)
    assert str(ds) ==  "data.0077.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_tb():
    ds = utils.data_dir_load(tb)
    assert str(ds) == "data.0005.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_iso():
    ds = utils.data_dir_load(iso)
    assert str(ds) == "data.0000.3d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_zp():
    ds = utils.data_dir_load(zp)
    assert str(ds) == "plt32.2d.hdf5"
    return ds

@pytest.fixture(scope='class')
def ds_kho():
    ds = utils.data_dir_load(kho)
    assert str(ds) == "data.0004.hdf5"
    return ds
