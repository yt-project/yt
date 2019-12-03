"""
Title: conftest.py
Purpose: Generates parameters and loads data for tests.
"""
import pytest

from yt.utilities.answer_testing import utils


# Data files
blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2D0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"


# Global test parameters
axes = [0, 1, 2]
center = "max"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "density"]


# Tests that use the hashing fixture
hash_tests = [
    'test_bw_polar_2d',
    'test_blastwave_cartesian_3D',
    'test_blaswave_spherical_2D',
    'test_blastwave_cylindrical_3D',
    'test_khi_cartesian_2D',
    'test_khi_cartesian_3D',
    'test_jet_cylindrical_25D',
    'test_riemann_cartesian_175D',
]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_bw_polar_2d':
        fields = utils.data_dir_load(blastwave_polar_2D).field_list
    if metafunc.function.__name__ == 'test_blastwave_cartesian_3D':
        fields = utils.data_dir_load(blastwave_cartesian_3D).field_list
    if metafunc.function.__name__ == 'test_blastwave_spherical_2D':
        fields = utils.data_dir_load(blastwave_spherical_2D).field_list
    if metafunc.function.__name__ == 'test_blastwave_cylindrical_3D':
        fields = utils.data_dir_load(blastwave_cylindrical_3D).field_list
    if metafunc.function.__name__ == 'test_khi_cartesian_2D':
        fields = utils.data_dir_load(khi_cartesian_2D).field_list
    if metafunc.function.__name__ == 'test_khi_cartesian_3D':
        fields = utils.data_Dir_load(khi_cartesian_3D).field_list
    if metafunc.function.__name__ == 'test_jet_cylindrical_25D':
        fields = utils.data_dir_load(jet_cylindrical_25D).field_list
    if metafunc.function.__name__ == 'test_riemann_cartesian_175D':
        fields = utils.data_dir_load(riemann_cartesian_175D).field_list
    if metafunc.function.__name__ in hash_tests:
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])
        metafunc.parametrize('f', fields, ids=[str(i) for i in range(len(fields))])


@pytest.fixture(scope='class')
def ds_bw_polar_2D():
    ds = utils.data_dir_load(blastwave_polar_2D)
    return ds

@pytest.fixture(scope='class')
def ds_blastwave_cartesian_3D():
    ds = utils.data_dir_load(blastwave_cartesian_3D)
    return ds

@pytest.fixture(scope='class')
def ds_blastwave_spherical_2D():
    ds = utils.data_dir_load(blastwave_spherical_2D)
    return ds

@pytest.fixture(scope='class')
def ds_blastwave_cylindrical_3D():
    ds = utils.data_dir_load(blastwave_cylindrical_3D)
    return ds

@pytest.fixture(scope='class')
def ds_khi_cartesian_2D():
    ds = utils.data_dir_load(khi_cartesian_2D)
    return ds

@pytest.fixture(scope='class')
def ds_khi_cartesian_3D():
    ds = utils.data_dir_load(khi_cartesian_3D)
    return ds

@pytest.fixture(scope='class')
def ds_jet_cylindrical_25D():
    ds = utils.data_dir_load(jet_cylindrical_25D)
    return ds

@pytest.fixture(scope='class')
def ds_riemann_cartesian_175D():
    ds = utils.data_dir_load(riemann_cartesian_175D)
    return ds
