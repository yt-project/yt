"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"




def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_datacontainer_data':
        fields = [('grid', 'density'), ('all', 'particle_mass')]
        metafunc.parametrize('field', fields, ids=['density', 'mass'])
    if metafunc.function.__name__ == 'test_covering_grid_datacontainer_data':
        fields = [('grid', 'density'), ('all', 'particle_mass')]
        metafunc.parametrize('field', fields, ids=['density', 'mass'])
    if metafunc.function.__name__ == 'test_arbitrary_grid_datacontainer_data':
        fields = [('grid', 'density'), ('all', 'particle_mass')]
        metafunc.parametrize('field', fields, ids=['density', 'mass'])
    if metafunc.function.__name__ == 'test_frb_datacontainer_data':
        fields = ['density']
        metafunc.parametrize('field', fields, ids=['density'])
    if metafunc.function.__name__ == 'test_spatial_data':
        fields = [('grid', 'density')]
        metafunc.parametrize('field', fields, ids=['density'])
    if metafunc.function.__name__ == 'test_profile_data1':
        fields = ['temperature', 'x', 'density']
        metafunc.parametrize('field', fields, ids=['temperature', 'x', 'density'])
    if metafunc.function.__name__ == 'test_profile_data2':
        fields = ['density', 'x', 'temperature', 'y', 'cell_mass']
        metafunc.parametrize('field', fields, ids=['density', 'x', 'temperature',
            'y', 'cell_mass'])
    if metafunc.function.__name__ == 'test_nonspatial_data1':
        fields = ['region_density', 'sphere_density']
        metafunc.parametrize('field', fields, ids=['region_density', 'sphere_density'])
    if metafunc.function.__name__ == 'test_nonspatial_data2':
        fields = ['density']
        metafunc.parametrize('field', fields, ids=['density''])



@pytest.fixture(scope='class')
def ds_enzotiny():
    ds = utils.data_dir_load(enzotiny)
    return ds
