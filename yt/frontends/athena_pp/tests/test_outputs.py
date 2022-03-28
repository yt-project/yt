import numpy as np

from yt.frontends.athena_pp.api import AthenaPPDataset
from yt.loaders import load
from yt.testing import (
    assert_allclose_units,
    assert_equal,
    disable_dataset_cache,
    requires_file,
    units_override_check,
)
from yt.units import dimensions
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

# Deactivating this problematic test until the dataset type can be
# handled properly, see https://github.com/yt-project/yt/issues/3619
"""
_fields_disk = ("density", "velocity_r")

disk = "KeplerianDisk/disk.out1.00000.athdf"


@requires_ds(disk)
def test_disk():
    ds = data_dir_load(disk)
    assert_equal(str(ds), "disk.out1.00000")
    dd = ds.all_data()
    vol = (ds.domain_right_edge[0] ** 3 - ds.domain_left_edge[0] ** 3) / 3.0
    vol *= np.cos(ds.domain_left_edge[1]) - np.cos(ds.domain_right_edge[1])
    vol *= ds.domain_right_edge[2].v - ds.domain_left_edge[2].v
    assert_allclose(dd.quantities.total_quantity(("gas", "cell_volume")), vol)
    for field in _fields_disk:

        def field_func(name):
            return dd[field]

        yield GenericArrayTest(ds, field_func, args=[field])
"""

_fields_AM06 = ("temperature", "density", "velocity_magnitude", "magnetic_field_x")

AM06 = "AM06/AM06.out1.00400.athdf"


@requires_ds(AM06)
def test_AM06():
    ds = data_dir_load(AM06)
    assert_equal(str(ds), "AM06.out1.00400")
    for test in small_patch_amr(ds, _fields_AM06):
        test_AM06.__name__ = test.description
        yield test


uo_AM06 = {
    "length_unit": (1.0, "kpc"),
    "mass_unit": (1.0, "Msun"),
    "time_unit": (1.0, "Myr"),
}


@requires_file(AM06)
def test_AM06_override():
    # verify that overriding units causes derived unit values to be updated.
    # see issue #1259
    ds = load(AM06, units_override=uo_AM06)
    assert_equal(float(ds.magnetic_unit.in_units("gauss")), 9.01735778342523e-08)


@requires_file(AM06)
def test_units_override():
    units_override_check(AM06)


@requires_file(AM06)
def test_AthenaPPDataset():
    assert isinstance(data_dir_load(AM06), AthenaPPDataset)


@requires_file(AM06)
def test_magnetic_units():
    ds = load(AM06, unit_system="code")
    assert ds.magnetic_unit.units.dimensions == dimensions.magnetic_field_cgs
    assert (ds.magnetic_unit**2).units.dimensions == dimensions.pressure


@requires_file(AM06)
@disable_dataset_cache
def test_mag_factor():
    ds1 = load(AM06, units_override=uo_AM06, magnetic_normalization="gaussian")
    assert ds1.magnetic_unit == np.sqrt(
        4.0 * np.pi * ds1.mass_unit / (ds1.time_unit**2 * ds1.length_unit)
    )
    sp1 = ds1.sphere("c", (100.0, "kpc"))
    pB1a = (
        sp1["athena_pp", "Bcc1"] ** 2
        + sp1["athena_pp", "Bcc2"] ** 2
        + sp1["athena_pp", "Bcc3"] ** 2
    ) / (8.0 * np.pi)
    pB1b = (
        sp1["gas", "magnetic_field_x"] ** 2
        + sp1["gas", "magnetic_field_y"] ** 2
        + sp1["gas", "magnetic_field_z"] ** 2
    ) / (8.0 * np.pi)
    pB1a.convert_to_units("dyn/cm**2")
    pB1b.convert_to_units("dyn/cm**2")
    assert_allclose_units(pB1a, pB1b)
    assert_allclose_units(pB1a, sp1["magnetic_pressure"])
    ds2 = load(AM06, units_override=uo_AM06, magnetic_normalization="lorentz_heaviside")
    assert ds2.magnetic_unit == np.sqrt(
        ds2.mass_unit / (ds2.time_unit**2 * ds2.length_unit)
    )
    sp2 = ds2.sphere("c", (100.0, "kpc"))
    pB2a = (
        sp2["athena_pp", "Bcc1"] ** 2
        + sp2["athena_pp", "Bcc2"] ** 2
        + sp2["athena_pp", "Bcc3"] ** 2
    ) / 2.0
    pB2b = (
        sp2["gas", "magnetic_field_x"] ** 2
        + sp2["gas", "magnetic_field_y"] ** 2
        + sp2["gas", "magnetic_field_z"] ** 2
    ) / 2.0
    pB2a.convert_to_units("dyn/cm**2")
    pB2b.convert_to_units("dyn/cm**2")
    assert_allclose_units(pB2a, pB2b)
    assert_allclose_units(pB2a, sp2["magnetic_pressure"])
    assert_allclose_units(pB1a, pB2a)
