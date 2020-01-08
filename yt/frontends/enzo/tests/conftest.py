"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
toro1d = "ToroShockTube/DD0001/data0001"
kh2d = "EnzoKelvinHelmholtz/DD0011/DD0011"
m7 = "DD0010/moving7_0010"
g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
ecp = "enzo_cosmology_plus/DD0046/DD0046"
two_sphere_test = 'ActiveParticleTwoSphere/DD0011/DD0011'
active_particle_cosmology = 'ActiveParticleCosmology/DD0046/DD0046'
mhdctot = "MHDCTOrszagTang/DD0004/data0004"
dnz = "DeeplyNestedZoom/DD0025/data0025"
p3mini = "PopIII_mini/DD0034/DD0034"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_toro1d' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [utils.data_dir_load(toro1d).field_list,
            [str(i) for i in range(len(utils.data_dir_load(toro1d).field_list))]]  
    },
    'test_kh2d' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [utils.data_dir_load(kh2d).field_list,
            [str(i) for i in range(len(utils.data_dir_load(kh2d).field_list))]]  
    },
    'test_moving7' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('temperature', 'density', 'velocity_magnitude', 'velocity_divergence'),
            ('temperature', 'density', 'velocity_magnitude', 'velocity_divergence')]
    },
    'test_galaxy0030' : {
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'density')],
        'f' : [('temperature', 'density', 'velocity_magnitude', 'velocity_divergence'),
            ('temperature', 'density', 'velocity_magnitude', 'velocity_divergence')]
    },
    'test_simulated_halo_mass_function' : {
        'finder' : [('fof', 'hop'), ('fof', 'hop')]
    },
    'test_analytic_halo_mass_function' : {
        'fit' : [[i for i in range(1,6)], [str(i) for i in range(1,6)]]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])

@pytest.fixture(scope='class')
def ds_toro1d():
    ds = utils.data_dir_load(toro1d)
    assert str(ds) == 'data0001'
    return ds

@pytest.fixture(scope='class')
def ds_kh2d():
    ds = utils.data_dir_load(kh2d)
    assert str(ds) == 'DD0011'
    return ds

@pytest.fixture(scope='class')
def ds_m7():
    ds = utils.data_dir_load(m7)
    assert str(ds) == 'moving7_0010'
    return ds

@pytest.fixture(scope='class')
def ds_g30():
    ds = utils.data_dir_load(g30)
    assert str(ds) == 'galaxy0030'
    return ds

@pytest.fixture(scope='class')
def ds_enzotiny():
    ds = utils.data_dir_load(enzotiny)
    return ds

@pytest.fixture(scope='class')
def ds_ecp():
    ds = utils.data_dir_load(ecp)
    return ds

@pytest.fixture(scope='class')
def ds_two_sphere_test():
    ds = utils.data_dir_load(two_sphere_test)
    return ds

@pytest.fixture(scope='class')
def ds_active_particle_cosmology():
    ds = utils.data_dir_load(active_particle_cosmology)
    return ds

@pytest.fixture(scope='class')
def ds_mhdctot():
    ds = utils.data_dir_load(mhdctot)
    return ds

@pytest.fixture(scope='class')
def ds_dnz():
    ds = utils.data_dir_load(dnz)
    return ds

@pytest.fixture(scope='class')
def ds_p3mini():
    ds = utils.data_dir_load(p3mini)
    return ds
