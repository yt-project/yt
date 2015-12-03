import numpy as np
from yt.frontends.stream.api import load_uniform_grid
from yt.utilities.physical_constants import mu_0
from yt.testing import assert_array_almost_equal_nulp

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_magnetic_fields():

    ddims = (32,32,32)

    data_cgs = {"magnetic_field_x":(np.random.random(size=ddims), "gauss"),
                "magnetic_field_y":(np.random.random(size=ddims), "gauss"),
                "magnetic_field_z":(np.random.random(size=ddims), "gauss")}
    data_mks = {"magnetic_field_x":(np.random.random(size=ddims), "T"),
                "magnetic_field_y":(np.random.random(size=ddims), "T"),
                "magnetic_field_z":(np.random.random(size=ddims), "T")}

    ds_cgs = load_uniform_grid(data_cgs, ddims)
    ds_mks = load_uniform_grid(data_mks, ddims, unit_system="mks")

    ds_cgs.index
    ds_mks.index

    dd_cgs = ds_cgs.all_data()
    dd_mks = ds_mks.all_data()

    assert ds_cgs.field_info["gas","magnetic_field_strength"].units == "gauss"
    assert ds_cgs.field_info["gas","magnetic_field_poloidal"].units == "gauss"
    assert ds_cgs.field_info["gas","magnetic_field_toroidal"].units == "gauss"
    assert ds_mks.field_info["gas","magnetic_field_strength"].units == "T"
    assert ds_mks.field_info["gas","magnetic_field_poloidal"].units == "T"
    assert ds_mks.field_info["gas","magnetic_field_toroidal"].units == "T"

    emag_cgs = (dd_cgs["magnetic_field_x"]**2 +
                dd_cgs["magnetic_field_y"]**2 +
                dd_cgs["magnetic_field_z"]**2)/(8.0*np.pi)

    emag_mks = (dd_mks["magnetic_field_x"]**2 +
                dd_mks["magnetic_field_y"]**2 +
                dd_mks["magnetic_field_z"]**2)/(2.0*mu_0)

    yield assert_array_almost_equal_nulp, emag_cgs, dd_cgs["magnetic_energy"], 2
    yield assert_array_almost_equal_nulp, emag_mks, dd_mks["magnetic_energy"], 2
