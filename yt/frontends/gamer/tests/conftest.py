"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
jet         = "InteractingJets/jet_000002"
psiDM       = "WaveDarkMatter/psiDM_000020"
plummer     = "Plummer/plummer_000000"
mhd_vortex   = "MHDOrszagTangVortex/Data_000018"


axes = [0, 1, 2]
center = "max"
ds_objs = [None, ("sphere", (center, (0.1, 'unitary')))]
weights = [None, "density"]
fields_mhdvortex = (("gamer","CCMagX"), ("gamer","CCMagY"), ("gas","magnetic_energy"))


def pytest_generate_tests(metafunc):
    param_tests = ['test_jet', 'test_psiDM', 'test_plummer']
    if metafunc.function.__name__ == 'test_jet':
        fields = ("temperature", "density", "velocity_magnitude")
        ids = ['temperature', 'density', 'velocity_magnitude']
    if metafunc.function.__name__ == 'test_psiDM':
        fields = ("Dens", "Real", "Imag")
        ids = ['dens', 'real', 'imag']
    if metafunc.function.__name__ == 'test_plummer':
        fields = ( ("gamer","ParDens"), ("deposit","io_cic") )
        ids = ['pardens', 'io_cic']
    if metafunc.function.__name__ == 'test_mhdvortex':
        fields = fields_mhdvortex
        ids = ['CCMagX', 'CCMagY', 'magnetic_energy']
    if metafunc.function.__name__ in param_tests:
        metafunc.parametrize('f', fields, ids=ids)
        metafunc.parametrize('a', axes, ids=['0', '1', '2'])
        metafunc.parametrize('d', ds_objs, ids=['None', 'sphere'])
        metafunc.parametrize('w', weights, ids=['None', 'density'])


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

@pytest.fixture(scope='class')
def ds_plummer():
    ds = utils.data_dir_load(plummer)
    assert str(ds) == "plummer_000000"
    return ds

@pytest.fixture(scope='class')
def ds_mhd_vortex():
    ds = utils.data_dir_load(mhd_vortex)
    assert str(ds) == "Data_000018"
    return ds
