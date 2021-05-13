import numpy as np

from yt.loaders import load_uniform_grid
from yt.testing import assert_almost_equal
from yt.utilities.physical_constants import mu_0


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_magnetic_fields():

    ddims = (16, 16, 16)
    data1 = {
        "magnetic_field_x": (np.random.random(size=ddims), "T"),
        "magnetic_field_y": (np.random.random(size=ddims), "T"),
        "magnetic_field_z": (np.random.random(size=ddims), "T"),
    }
    data2 = {}
    for field in data1:
        data2[field] = (data1[field][0] * 1.0e4, "gauss")

    ds1 = load_uniform_grid(data1, ddims, unit_system="cgs")
    ds2 = load_uniform_grid(data2, ddims, unit_system="mks")
    # For this test dataset, code units are cgs units
    ds3 = load_uniform_grid(data2, ddims, unit_system="code")

    ds1.index
    ds2.index
    ds3.index

    dd1 = ds1.all_data()
    dd2 = ds2.all_data()
    dd3 = ds3.all_data()

    assert ds1.fields.gas.magnetic_field_strength.units == "G"
    assert ds1.fields.gas.magnetic_field_poloidal.units == "G"
    assert ds1.fields.gas.magnetic_field_toroidal.units == "G"
    assert ds2.fields.gas.magnetic_field_strength.units == "T"
    assert ds2.fields.gas.magnetic_field_poloidal.units == "T"
    assert ds2.fields.gas.magnetic_field_toroidal.units == "T"
    assert ds3.fields.gas.magnetic_field_strength.units == "code_magnetic"
    assert ds3.fields.gas.magnetic_field_poloidal.units == "code_magnetic"
    assert ds3.fields.gas.magnetic_field_toroidal.units == "code_magnetic"

    emag1 = (
        dd1[("gas", "magnetic_field_x")] ** 2
        + dd1[("gas", "magnetic_field_y")] ** 2
        + dd1[("gas", "magnetic_field_z")] ** 2
    ) / (8.0 * np.pi)
    emag1.convert_to_units("dyne/cm**2")

    emag2 = (
        dd2[("gas", "magnetic_field_x")] ** 2
        + dd2[("gas", "magnetic_field_y")] ** 2
        + dd2[("gas", "magnetic_field_z")] ** 2
    ) / (2.0 * mu_0)
    emag2.convert_to_units("Pa")

    emag3 = (
        dd3[("gas", "magnetic_field_x")] ** 2
        + dd3[("gas", "magnetic_field_y")] ** 2
        + dd3[("gas", "magnetic_field_z")] ** 2
    ) / (2.0 * mu_0)
    emag3.convert_to_units("code_pressure")

    assert_almost_equal(emag1, dd1[("gas", "magnetic_energy_density")])
    assert_almost_equal(emag2, dd2[("gas", "magnetic_energy_density")])
    assert_almost_equal(emag3, dd3[("gas", "magnetic_energy_density")])

    assert str(emag1.units) == str(dd1[("gas", "magnetic_energy_density")].units)
    assert str(emag2.units) == str(dd2[("gas", "magnetic_energy_density")].units)
    assert str(emag3.units) == str(dd3[("gas", "magnetic_energy_density")].units)

    assert_almost_equal(emag1.in_cgs(), emag2.in_cgs())
    assert_almost_equal(emag1.in_cgs(), emag3.in_cgs())
