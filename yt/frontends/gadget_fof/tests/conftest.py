"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"
g298 = "gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5"
g56 = "gadget_halos/data/groups_056/fof_subhalo_tab_056.0.hdf5"
g76 = "gadget_halos/data/groups_076/fof_subhalo_tab_076.0.hdf5"


fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_velocity_x",
           "particle_velocity_y", "particle_velocity_z",
           "particle_mass", "particle_identifier")


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ in ['test_fields_g5', 'test_fields_g42']:
        metafunc.parametrize('field', fields, ids=['x', 'y', 'z', 'vx', 'vy',
            'vz', 'mass', 'identifier'])


@pytest.fixture(scope='class')
def ds_g5():
    ds = utils.data_dir_load(g5)
    return ds

@pytest.fixture(scope='class')
def ds_g42():
    ds = utils.data_dir_load(g42)
    return ds

@pytest.fixture(scope='class')
def ds_g298():
    ds = utils.data_dir_load(g298)
    return ds

@pytest.fixture(scope='class')
def ds_g56():
    ds = utils.data_dir_load(g56)
    return ds

@pytest.fixture(scope='class')
def ds_g76():
    ds = utils.data_dir_load(g76)
    return ds
