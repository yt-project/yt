"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.testing import assert_equal
from yt.utilities.answer_testing import utils


# Test data
radadvect = "RadAdvect/plt00000"
rt = "RadTube/plt00500"
star = "StarParticles/plrd01000"
LyA = "Nyx_LyA/plt00000"
RT_particles = "RT_particles/plt00050"
plasma = "PlasmaAcceleration/plt00030_v2"
langmuir = "LangmuirWave/plt00020_v2"
beam = "GaussianBeam/plt03008"
raw_fields = "Laser/plt00015"


orion_tests = ['test_radadvect', 'test_radtube', 'test_star']
warpx_tests = ['test_langmuir', 'test_plasma', 'test_beam']
axes = [0, 1, 2]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ in orion_tests:
        weights = [None, "density"]
        _orion_fields = ("temperature", "density", "velocity_magnitude")
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        metafunc.parametrize('f', _orion_fields, ids=['temperature', 'density',
            'velocity_magnitude'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
    if metafunc.function.__name__ == 'test_LyA':
        _nyx_fields = ("Ne", "Temp", "particle_mass_density")
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ne"]
        metafunc.parametrize('f', _nyx_fields, ids=['Ne', 'Temp', 'particle_mass_density'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'Ne'])
    if metafunc.function.__name__ == 'test_RT_particles':
        _castro_fields = ("Temp", "density", "particle_count")
        center = "max"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "density"]
        metafunc.parametrize('f', _castro_fields, ids=['Temp', 'density', 'count'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
    if metafunc.function.__name__ in warpx_tests:
        _warpx_fields = ("Ex", "By", "jz")
        center = "c"
        ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
        weights = [None, "Ex"]
        metafunc.parametrize('f', _warpx_fields, ids=['Ex', 'By', 'jz'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'Ex'])
    if metafunc.function.__name__ == 'test_raw_fields':
        _raw_fields = [('raw', 'Bx'), ('raw', 'Ey'), ('raw', 'jz')]
        metafunc.parametrize('f', _raw_fields, ids=['Bx', 'Ey', 'jz'])


@pytest.fixture(scope='class')
def ds_radadvect():
    ds = utils.data_dir_load(radadvect)
    assert (str(ds), "plt00000")
    return ds

@pytest.fixture(scope='class')
def ds_rt():
    ds = utils.data_dir_load(rt)
    assert str(ds) == "plt00500"
    return ds

@pytest.fixture(scope='class')
def ds_star():
    ds = utils.data_dir_load(star)
    assert str(ds) == "plrd01000"
    return ds

@pytest.fixture(scope='class')
def ds_LyA():
    ds = utils.data_dir_load(LyA)
    assert str(ds) == "plt00000"
    return ds

@pytest.fixture(scope='class')
def ds_RT_particles():
    ds = utils.data_dir_load(RT_particles)
    assert str(ds) == "plt00050"
    return ds

@pytest.fixture(scope='class')
def ds_plasma():
    ds = utils.data_dir_load(plasma)
    assert str(ds) == "plt00030_v2"
    return ds

@pytest.fixture(scope='class')
def ds_langmuir():
    ds = utils.data_dir_load(langmuir)
    assert_equal(str(ds), "plt00020_v2")
    return ds

@pytest.fixture(scope='class')
def ds_beam():
    ds = utils.data_dir_load(beam)
    assert (str(ds), "plt03008")
    return ds

@pytest.fixture(scope='class')
def ds_raw_fields():
    ds = utils.data_dir_load(raw_fields)
    return ds
