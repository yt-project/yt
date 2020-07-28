"""
Title: conftest.py
Purpose: Generates parameters and loads data for tests.
"""
import pytest

from yt.utilities.answer_testing.utils import data_dir_load
from yt.utilities.exceptions import YTOutputNotIdentified


def _get_fields_to_check(fname):
    # This function is called during test collection. If this frontend
    # is not being run, and therefore the data isn't present, this try
    # except block prevents pytest from failing needlessly
    try:
        ds = data_dir_load(fname)
        fields = [("gas", "density"), ("gas", "velocity_magnitude")]
        field_ids = ["density", "velocity_magnitude"]
        raw_fields_labels = [fname for ftype, fname in ds.field_list]
        if "b1" in raw_fields_labels:
            fields.append(("gas", "magnetic_energy_density"))
            field_ids.append("magnetic_energy_density")
        if "e" in raw_fields_labels:
            fields.append(("gas", "energy_density"))
            field_ids.append("energy_density")
        if "rhod1" in raw_fields_labels:
            fields.append(("gas", "total_dust_density"))
            field_ids.append("total_dust_density")
            # note : not hitting dust velocity fields
        return [fields, field_ids]
    except YTOutputNotIdentified:
        return [[pytest.param(None, marks=pytest.mark.skip),], ['Data not found',]]


bw_polar_2d = _get_fields_to_check("amrvac/bw_polar_2D0000.dat")
bw_cart_3d = _get_fields_to_check("amrvac/bw_3d0000.dat")
bw_sph_2d = _get_fields_to_check("amrvac/bw_2d0000.dat")
bw_cyl_3d = _get_fields_to_check("amrvac/bw_cylindrical_3D0000.dat")
khi_cart_2d = _get_fields_to_check("amrvac/kh_2d0000.dat")
khi_cart_3d = _get_fields_to_check("amrvac/kh_3D0000.dat")
jet_cyl_25d = _get_fields_to_check("amrvac/Jet0003.dat")
rie_cart_175d = _get_fields_to_check("amrvac/R_1d0005.dat")
rmi_cart_dust_2d = _get_fields_to_check(
    "amrvac/Richtmyer_Meshkov_dust_2D/RM2D_dust_Kwok0000.dat"
)

test_params = {
    "test_bw_polar_2d": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": bw_polar_2d,
    },
    "test_blastwave_cartesian_3D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": bw_cart_3d,
    },
    "test_blastwave_spherical_2D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": bw_sph_2d,
    },
    "test_blastwave_cylindrical_3D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": bw_cyl_3d,
    },
    "test_khi_cartesian_2D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": khi_cart_2d,
    },
    "test_khi_cartesian_3D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": khi_cart_3d,
    },
    "test_jet_cylindrical_25D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": jet_cyl_25d,
    },
    "test_riemann_cartesian_175D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": rie_cart_175d,
    },
    "test_rmi_cartesian_dust_2D": {
        "a": [(0, 1, 2), ("0", "1", "2")],
        "d": [(None, ("sphere", ("max", (0.1, "unitary")))), ("None", "sphere")],
        "w": [(None, "density"), ("None", "density")],
        "f": rmi_cart_dust_2d,
    },
}


def pytest_generate_tests(metafunc):
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])
