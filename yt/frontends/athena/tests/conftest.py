"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
cloud = "ShockCloud/id0/Cloud.0050.vtk"
blast = "MHDBlast/id0/Blast.0100.vtk"
stripping = "RamPressureStripping/id0/rps.0062.vtk"

axes = [0, 1, 2]
center = "max"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "density"]


def pytest_generate_tests(metafunc):
    param_tests = ['test_cloud', 'test_blast', 'test_stripping']
    if metafunc.function.__name__ == 'test_cloud':
        fields = ("scalar[0]", "density", "total_energy")
        metafunc.parametrize('f', fields, ids=['scalar0', 'density', 'total_energy'])
    if metafunc.function.__name__ == 'test_blast':
        fields = ("temperature", "density", "velocity_magnitude")
        metafunc.parametrize('f', fields, ids=['temperature', 'density', 'velocity_magnitude'])
    if metafunc.function.__name__ == 'test_stripping':
        fields = ("temperature", "density", "specific_scalar[0]")
        metafunc.parametrize('f', fields, ids=['temperature', 'density', 'specific_scalar0'])
    if metafunc.function.__name__ in param_tests:
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])


@pytest.fixture(scope='class')
def ds_cloud():
    ds = utils.data_dir_load(cloud)
    assert_equal(str(ds), "Cloud.0050")
    return ds

@pytest.fixture(scope='class')
def ds_blast():
    ds = utils.data_dir_load(blast)
    assert_equal(str(ds), "Blast.0100")
    return ds

@pytest.fixture(scope='class')
def ds_stripping():
    uo_stripping = {"time_unit":3.086e14,
                    "length_unit":8.0236e22,
                    "mass_unit":9.999e-30*8.0236e22**3}
    ds = utils.data_dir_load(stripping, kwargs={"units_override":uo_stripping})
    assert_equal(str(ds), "rps.0062")
    return ds
