import numpy as np
from yt.utilities.physical_constants import mu_0
from yt.testing import assert_almost_equal, fake_random_ds

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_magnetic_fields():

    fields = ("magnetic_field_x", "magnetic_field_y", "magnetic_field_z")
    cgs_units = ("gauss",)*3
    mks_units = ("T",)*3

    ds_cgs = fake_random_ds(16, fields=fields, units=cgs_units, nprocs=16)
    ds_mks = fake_random_ds(16, fields=fields, units=mks_units, nprocs=16,
                            unit_system="mks")

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
    emag_cgs.convert_to_units("dyne/cm**2")

    emag_mks = (dd_mks["magnetic_field_x"]**2 +
                dd_mks["magnetic_field_y"]**2 +
                dd_mks["magnetic_field_z"]**2)/(2.0*mu_0)
    emag_mks.convert_to_units("Pa")

    yield assert_almost_equal, emag_cgs, dd_cgs["magnetic_energy"]
    yield assert_almost_equal, emag_mks, dd_mks["magnetic_energy"]

    assert str(emag_cgs.units) == str(dd_cgs["magnetic_energy"].units)
    assert str(emag_mks.units) == str(dd_mks["magnetic_energy"].units)

