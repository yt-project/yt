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
