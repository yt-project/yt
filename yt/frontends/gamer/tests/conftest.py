"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing.utils import data_dir_load

# Test data
jet = "InteractingJets/jet_000002"


test_params = {
    "test_jet": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [
            ("temperature", "density", "velocity_magnitude"),
            ("temperature", "density", "velocity_magnitude"),
        ],
    },
    "test_psiDM": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [("Dens", "Real", "Imag"), ("Dens", "Real", "Imag")],
    },
    "test_plummer": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [(("gamer", "ParDens"), ("deposit", "io_cic")), ("parDens", "iocic")],
    },
    "test_mhdvortex": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": [
            (("gamer", "CCMagX"), ("gamer", "CCMagY"), ("gas", "magnetic_energy")),
            ("CCMagX", "CCMagY", "magnetic_energy"),
        ],
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


@pytest.fixture(scope="class")
def ds_jet():
    jet_units = {
        "length_unit": (1.0, "kpc"),
        "time_unit": (3.08567758096e13, "s"),
        "mass_unit": (1.4690033e36, "g"),
    }
    ds = data_dir_load(jet, kwargs={"units_override": jet_units})
    assert str(ds) == "jet_000002"
    return ds
