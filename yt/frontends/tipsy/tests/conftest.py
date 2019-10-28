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


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ in ['test_pkdgrav', 'test_gasoline_dm_only']:
        dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
        axes = [0, 1, 2]
        weights = [None]
        metafunc.parametrize('f', _fields, ids=['all_density', 'all_count', 'DM_density'])
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None'])
    if metafunc.function.__name__ == 'test_tipsy_galaxy':
        dso = [None, ('sphere', ('c', (0.1, 'unitary')))]
        metafunc.parametrize('f,w', [(pair[0], pair[1]) for pair in tg_fields],
            ids=['dens', 'temp-None', 'temp-dens', 'velocity_magnitude',
            'Fe_fraction', 'Metals'])
        metafunc.parametrize('d', dso, ids=['None', 'sphere'])
        metafunc.parametrize('a', [0, 1, 2], ids=['0', '1', '2'])


@pytest.fixture(scope='class')
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
