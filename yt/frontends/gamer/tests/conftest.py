"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
jet         = "InteractingJets/jet_000002"
psiDM       = "WaveDarkMatter/psiDM_000020"


@pytest.fixture(scope='class')
def ds_jet():
    jet_units   = {"length_unit":(1.0,"kpc"),
                   "time_unit"  :(3.08567758096e+13,"s"),
                   "mass_unit"  :(1.4690033e+36,"g")}
    ds = utils.data_dir_load(jet, kwargs={"units_override":jet_units})
    assert str(ds) == "jet_000002"
    return ds

@pytest.fixture(scope='class')
def ds_psiDM():
    ds = utils.data_dir_load(psiDM)
    assert str(ds) == "psiDM_000020"
    return ds
