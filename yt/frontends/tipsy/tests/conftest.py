"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.frontends.tipsy.api import TipsyDataset
from yt.utilities.answer_testing import utils


# Test data
pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = 'TipsyGalaxy/galaxy.00300'

_fields = (("deposit", "all_density"),
           ("deposit", "all_count"),
           ("deposit", "DarkMatter_density"),
)
tg_fields = [
        [('gas', 'density'), None],
        [('gas', 'temperature'), None],
        [('gas', 'temperature'), ('gas', 'density')],
        [('gas', 'velocity_magnitude'), None],
        [('gas', 'Fe_fraction'), None],
        [('Stars', 'Metals'), None]
    ]


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'tes_pkdgrav' : {
        'f' : [_fields, ('all_density', 'all_count', 'DM_density')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.3, 'unitary')))), ('None', 'sphere')],
        'w' : [(None,), ('None',)]
    },
    'tes_gasoline_dmonly' : {
        'f' : [_fields, ('all_density', 'all_count', 'DM_density')],
        'a' : [(0, 1, 2), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('c', (0.3, 'unitary')))), ('None', 'sphere')],
        'w' : [(None,), ('None',)]
    },
    'test_tipsy_galaxy' : {
        'f, w' : [((pair[0], pair[1]) for pair in tg_fields), ('dens', 'temp-None',
                    'temp_dens', 'velocity_magnitude', 'Fe_fraction', 'Metals')],
        'd' : [(None, ('sphere', ('c', (0.1, 'unitary')))), ('None', 'sphere')],
        'a' : [(0, 1, 2), ('0', '1', '2')]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


@pytest.fixture(scope='function')
def ds_pkdgrav():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(field_dtypes = {"Coordinates": "d"},
                  cosmology_parameters = cosmology_parameters,
                  unit_base = {'length': (60.0, "Mpccm/h")},
                  n_ref = 64)
    ds = utils.data_dir_load(pkdgrav, TipsyDataset, (), kwargs)
    assert str(ds) == "halo1e11_run1.00400"
    return ds

@pytest.fixture(scope='class')
def ds_gasoline_dmonly():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(cosmology_parameters = cosmology_parameters,
                  unit_base = {'length': (60.0, "Mpccm/h")},
                  n_ref = 64)
    ds = utils.data_dir_load(gasoline_dmonly, TipsyDataset, (), kwargs)
    assert str(ds) == "agora_1e11.00400"
    return ds

@pytest.fixture(scope='class')
def ds_tipsy_gal():
    ds = utils.data_dir_load(tipsy_gal)
    return ds
