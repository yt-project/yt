"""
Test unit systems. 

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.units.unit_object import Unit, unit_system_registry
from yt.units.unit_systems import UnitSystem
from yt.units import dimensions
from yt.convenience import load
from yt.testing import assert_almost_equal, assert_allclose, requires_file
from yt.config import ytcfg

def test_unit_systems():
    goofy_unit_system = UnitSystem("goofy", "ly", "lbm", "hr", temperature_unit="R",
                                   angle_unit="arcsec", current_mks_unit="mA")
    assert goofy_unit_system["temperature"] == Unit("R")
    assert goofy_unit_system[dimensions.solid_angle] == Unit("arcsec**2")
    assert goofy_unit_system["energy"] == Unit("lbm*ly**2/hr**2")
    goofy_unit_system["energy"] = "eV"
    assert goofy_unit_system["energy"] == Unit("eV")
    assert goofy_unit_system["magnetic_field_mks"] == Unit("lbm/(hr**2*mA)")
    assert "goofy" in unit_system_registry

test_units = {}
test_units["mks"] = {"density": "kg/m**3",
                     "kinetic_energy": "Pa",
                     "velocity_magnitude": "m/s",
                     "velocity_divergence": "1/s",
                     "density_gradient_x": "kg/m**4"}
test_units["imperial"] = {"density": "lbm/ft**3",
                          "kinetic_energy": "lbf/ft**2",
                          "velocity_magnitude": "ft/s",
                          "velocity_divergence": "1/s",
                          "density_gradient_x": "lbm/ft**4"}
test_units["galactic"] = {"density": "Msun/kpc**3",
                          "kinetic_energy": "Msun/(Myr**2*kpc)",
                          "velocity_magnitude": "kpc/Myr",
                          "velocity_divergence": "1/Myr",
                          "density_gradient_x": "Msun/kpc**4"}
test_units["code"] = {"density": "code_mass/code_length**3",
                      "kinetic_energy": "code_pressure",
                      "velocity_magnitude": "code_velocity",
                      "velocity_divergence": "code_velocity/code_length",
                      "density_gradient_x": "code_mass/code_length**4"}

test_fields = ["density",
               "kinetic_energy",
               "velocity_divergence",
               "density_gradient_x",
               "velocity_magnitude"]

gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
@requires_file(gslr)
def test_fields_diff_systems_sloshing():
    ytcfg["yt","skip_dataset_cache"] = "True"

    ds_cgs = load(gslr)
    dd_cgs = ds_cgs.sphere("c", (15., "kpc"))

    for us in test_units:
        ds = load(gslr, unit_system=us)
        dd = ds.sphere("c", (15.,"kpc"))
        for field in test_fields:
            if us == "code":
                # For this dataset code units are cgs
                v1 = dd_cgs[field]
            else:
                v1 = dd_cgs[field].in_base(us)
            v2 = dd[field]
            assert_almost_equal(v1.v, v2.v)
            assert str(v2.units) == test_units[us][field]

etc = "enzo_tiny_cosmology/DD0046/DD0046"
@requires_file(etc)
def test_fields_diff_systems_etc():
    ytcfg["yt","skip_dataset_cache"] = "True"

    ds_cgs = load(etc)
    dd_cgs = ds_cgs.sphere("max", (500., "kpc"))

    for us in test_units:
        ds = load(etc, unit_system=us)
        dd = ds.sphere("max", (500.,"kpc"))
        for field in test_fields:
            if us == "code":
                v1 = dd_cgs[field].in_units(test_units["code"][field])
            else:
                v1 = dd_cgs[field].in_base(us)
            v2 = dd[field]
            assert_almost_equal(v1.v, v2.v)
            assert str(v2.units) == test_units[us][field]

wdm = 'WDMerger_hdf5_chk_1000/WDMerger_hdf5_chk_1000.hdf5'
@requires_file(wdm)
def test_tesla_magnetic_unit():
    ytcfg["yt", "skip_dataset_cache"] = "True"
    for us in ['cgs', 'mks', 'code']:
        ds = load(wdm, unit_system=us,
                  units_override={'magnetic_unit': (1.0, 'T')})
        ad = ds.all_data()
        dens = ad['density']
        magx = ad['magx']
        magnetic_field_x = ad['magnetic_field_x']

        if us == 'cgs':
            assert str(dens.units) == 'g/cm**3'
            assert str(magx.units) == 'code_magnetic'
            assert magx.uq == ds.quan(1e4, 'G')
            assert str(magnetic_field_x.units) == 'gauss'
            assert_allclose(magx.value, magnetic_field_x.value/1e4)
            assert_allclose(
                magnetic_field_x.to_equivalent('T', 'SI').value,
                magnetic_field_x.value/1e4)

        if us == 'mks':
            assert str(dens.units) == 'kg/m**3'
            assert str(magx.units) == 'code_magnetic'
            assert magx.uq == ds.quan(1, 'T')
            assert str(magnetic_field_x.units) == 'T'
            assert_allclose(magx.value, magnetic_field_x.value)
            assert_allclose(magnetic_field_x.to_equivalent('G', 'CGS').value,
                            magnetic_field_x.value*1e4)

        if us == 'code':
            assert str(dens.units) == 'code_mass/code_length**3'
            assert str(magx.units) == 'code_magnetic'
            assert magx.uq == ds.quan(1, 'T')
            assert str(magnetic_field_x.units) == 'code_magnetic'
            assert_allclose(magx.value, magnetic_field_x.value)
            assert_allclose(magnetic_field_x.to_equivalent('G', 'CGS').value,
                            magnetic_field_x.value*1e4)
