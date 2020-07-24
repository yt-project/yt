"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.frontends.fits.data_structures import (
    FITSDataset,
    SpectralCubeFITSDataset,
    SkyDataFITSDataset,
    EventsFITSDataset,
)
from yt.testing import assert_equal
from yt.utilities.answer_testing.utils import data_dir_load


# Test data
grs = "radio_fits/grs-50-cube.fits"
vf = "UnigridData/velocity_field_20.fits"
acis = "xray_fits/acisf05356N003_evt2.fits"
A2052 = "xray_fits/A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"


_fields_grs = ("temperature",)
_fields_vels = ("velocity_x", "velocity_y", "velocity_z")
_fields_acis = ("counts_0.1-2.0", "counts_2.0-5.0")
_fields_A2052 = ("flux",)


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    "test_grs": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "ones"), ("None", "ones")],
        "f": _fields_grs,
    },
    "test_velocity_field": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "ones"), ("None", "ones")],
        "f": _fields_vels,
    },
    "test_acis": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "ones"), ("None", "ones")],
        "f": _fields_acis,
    },
    "test_A2052": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("c", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "ones"), ("None", "ones")],
        "f": _fields_A2052,
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


@pytest.fixture(scope="class")
def ds_grs():
    ds = data_dir_load(grs, cls=SpectralCubeFITSDataset, kwargs={"nan_mask": 0.0})
    assert str(ds) == "grs-50-cube.fits"
    return ds


@pytest.fixture(scope="class")
def ds_vf():
    ds = data_dir_load(vf, cls=FITSDataset)
    assert str(ds) == "velocity_field_20.fits"
    return ds


@pytest.fixture(scope="class")
def ds_acis():
    ds = data_dir_load(acis, cls=EventsFITSDataset)
    from yt.frontends.fits.misc import setup_counts_fields

    ebounds = [(0.1, 2.0), (2.0, 5.0)]
    setup_counts_fields(ds, ebounds)
    assert_equal(str(ds), "acisf05356N003_evt2.fits")
    return ds


@pytest.fixture(scope="class")
def ds_A2052():
    ds = data_dir_load(A2052, cls=SkyDataFITSDataset)
    assert str(ds) == "A2052_merged_0.3-2_match-core_tmap_bgecorr.fits"
    return ds
