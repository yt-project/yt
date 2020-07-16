import os

import numpy as np

import yt
from yt.config import ytcfg
from yt.frontends.ramses.api import RAMSESDataset
from yt.frontends.ramses.field_handlers import DETECTED_FIELDS, HydroFieldFileHandler
from yt.testing import (
    assert_equal,
    requires_file,
    requires_module,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    PixelizedProjectionValuesTest,
    create_obj,
    data_dir_load,
    requires_ds,
)
from yt.utilities.on_demand_imports import _f90nml as f90nml

_fields = ("temperature", "density", "velocity_magnitude")

output_00080 = "output_00080/info_00080.txt"


@requires_ds(output_00080)
def test_output_00080():
    ds = data_dir_load(output_00080)
    assert_equal(str(ds), "info_00080")
    assert_equal(ds.particle_type_counts, {"io": 1090895, "nbody": 0})
    dso = [None, ("sphere", ("max", (0.1, "unitary")))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        output_00080, axis, field, weight_field, dobj_name
                    )
            yield FieldValuesTest(output_00080, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)


@requires_file(output_00080)
def test_RAMSESDataset():
    assert isinstance(data_dir_load(output_00080), RAMSESDataset)


@requires_file(output_00080)
def test_units_override():
    units_override_check(output_00080)


ramsesNonCosmo = "DICEGalaxyDisk_nonCosmological/output_00002/info_00002.txt"


@requires_file(ramsesNonCosmo)
def test_non_cosmo_detection():
    ds = yt.load(ramsesNonCosmo, cosmological=False)
    assert_equal(ds.cosmological_simulation, 0)

    ds = yt.load(ramsesNonCosmo, cosmological=None)
    assert_equal(ds.cosmological_simulation, 0)

    ds = yt.load(ramsesNonCosmo)
    assert_equal(ds.cosmological_simulation, 0)


@requires_file(ramsesNonCosmo)
def test_unit_non_cosmo():
    for force_cosmo in [False, None]:
        ds = yt.load(ramsesNonCosmo, cosmological=force_cosmo)

        expected_raw_time = 0.0299468077820411  # in ramses unit
        assert_equal(ds.current_time.value, expected_raw_time)

        expected_time = 14087886140997.336  # in seconds
        assert_equal(ds.current_time.in_units("s").value, expected_time)


ramsesCosmo = "output_00080/info_00080.txt"


@requires_file(ramsesCosmo)
def test_cosmo_detection():
    ds = yt.load(ramsesCosmo, cosmological=True)
    assert_equal(ds.cosmological_simulation, 1)

    ds = yt.load(ramsesCosmo, cosmological=None)
    assert_equal(ds.cosmological_simulation, 1)

    ds = yt.load(ramsesCosmo)
    assert_equal(ds.cosmological_simulation, 1)


@requires_file(ramsesCosmo)
def test_unit_cosmo():
    for force_cosmo in [True, None]:
        ds = yt.load(ramsesCosmo, cosmological=force_cosmo)

        expected_raw_time = 1.119216564055017  # in ramses unit
        assert_equal(ds.current_time.value, expected_raw_time)

        expected_time = 3.756241729312462e17  # in seconds
        assert_equal(ds.current_time.in_units("s").value, expected_time)


ramsesExtraFieldsSmall = "ramses_extra_fields_small/output_00001"


@requires_file(ramsesExtraFieldsSmall)
def test_extra_fields():
    extra_fields = [("particle_family", "I"), ("particle_pointer", "I")]
    ds = yt.load(
        os.path.join(ramsesExtraFieldsSmall, "info_00001.txt"),
        extra_particle_fields=extra_fields,
    )

    # the dataset should contain the fields
    for field, _ in extra_fields:
        assert ("all", field) in ds.field_list

    # Check the family (they should equal 100, for tracer particles)
    dd = ds.all_data()
    families = dd[("all", "particle_family")]
    assert all(families == 100)


@requires_file(ramsesExtraFieldsSmall)
def test_extra_fields_2():
    extra_fields = ["particle_extra_field_%s" % (i + 1) for i in range(2)]
    ds = yt.load(os.path.join(ramsesExtraFieldsSmall, "info_00001.txt"))

    # the dataset should contain the fields
    for field in extra_fields:
        assert ("io", field) in ds.field_list

    # In the dataset, the fields are integers, so we cannot test
    # that they are accessed correctly.


ramses_rt = "ramses_rt_00088/output_00088/info_00088.txt"


@requires_file(ramses_rt)
def test_ramses_rt():
    ds = yt.load(ramses_rt)
    ad = ds.all_data()

    expected_fields = [
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pres_IR",
        "Pressure",
        "Metallicity",
        "HII",
        "HeII",
        "HeIII",
    ]

    for field in expected_fields:
        assert ("ramses", field) in ds.field_list

        # test that field access works
        ad["ramses", field]

    # test that special derived fields for RT datasets work
    special_fields = [("gas", "temp_IR")]
    species = ["H_p1", "He_p1", "He_p2"]
    for specie in species:
        special_fields.extend(
            [
                ("gas", specie + "_fraction"),
                ("gas", specie + "_density"),
                ("gas", specie + "_mass"),
            ]
        )

    for field in special_fields:
        assert field in ds.derived_field_list
        ad[field]


ramses_sink = "ramses_sink_00016/output_00016/info_00016.txt"


@requires_file(ramses_sink)
def test_ramses_sink():
    expected_fields = [
        "BH_bondi_accretion",
        "BH_eddington_accretion",
        "BH_efficiency",
        "BH_esave",
        "BH_real_accretion",
        "BH_spin",
        "BH_spin_x",
        "BH_spin_y",
        "BH_spin_z",
        "gas_spin_x",
        "gas_spin_y",
        "gas_spin_z",
        "particle_birth_time",
        "particle_identifier",
        "particle_mass",
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_prop_0_0",
        "particle_prop_0_1",
        "particle_prop_0_2",
        "particle_prop_0_3",
        "particle_prop_1_0",
        "particle_prop_1_1",
        "particle_prop_1_2",
        "particle_velocity_x",
        "particle_velocity_y",
        "particle_velocity_z",
    ]

    # Check that sinks are autodetected
    ds = yt.load(ramses_sink)
    ad = ds.all_data()

    for field in expected_fields:
        assert ("sink", field) in ds.field_list

        # test that field access works
        ad["sink", field]

    # Checking that sinks are autodetected
    ds = yt.load(ramsesCosmo)
    ad = ds.all_data()

    for field in expected_fields:
        assert ("sink", "field") not in ds.field_list


ramses_new_format = "ramses_new_format/output_00002/info_00002.txt"


@requires_file(ramses_new_format)
def test_new_format():
    expected_particle_fields = [
        ("star", "particle_identity"),
        ("star", "particle_level"),
        ("star", "particle_mass"),
        ("star", "particle_metallicity"),
        ("star", "particle_position_x"),
        ("star", "particle_position_y"),
        ("star", "particle_position_z"),
        ("star", "particle_tag"),
        ("star", "particle_velocity_x"),
        ("star", "particle_velocity_y"),
        ("star", "particle_velocity_z"),
    ]

    ds = yt.load(ramses_new_format)
    ad = ds.all_data()

    # Check all the expected fields exist and can be accessed
    for f in expected_particle_fields:
        assert f in ds.derived_field_list
        ad[f]

    # Check there is only stars with tag 0 (it should be right)
    assert all(ad["star", "particle_family"] == 2)
    assert all(ad["star", "particle_tag"] == 0)
    assert len(ad["star", "particle_tag"]) == 600


@requires_file(ramses_sink)
def test_ramses_part_count():
    ds = yt.load(ramses_sink)
    pcount = ds.particle_type_counts

    assert_equal(pcount["io"], 17132, err_msg="Got wrong number of io particle")
    assert_equal(pcount["sink"], 8, err_msg="Got wrong number of sink particle")


@requires_file(ramsesCosmo)
def test_custom_particle_def():
    ytcfg.add_section("ramses-particles")
    ytcfg[
        "ramses-particles", "fields"
    ] = """particle_position_x, d
         particle_position_y, d
         particle_position_z, d
         particle_velocity_x, d
         particle_velocity_y, d
         particle_velocity_z, d
         particle_mass, d
         particle_identifier, i
         particle_refinement_level, I
         particle_birth_time, d
         particle_foobar, d
    """
    ds = yt.load(ramsesCosmo)

    def check_unit(array, unit):
        assert str(array.in_cgs().units) == unit

    try:
        assert ("io", "particle_birth_time") in ds.derived_field_list
        assert ("io", "particle_foobar") in ds.derived_field_list

        check_unit(ds.r["io", "particle_birth_time"], "s")
        check_unit(ds.r["io", "particle_foobar"], "dimensionless")
    finally:
        ytcfg.remove_section("ramses-particles")


@requires_file(ramsesCosmo)
def test_custom_hydro_def():
    ytcfg.add_section("ramses-hydro")
    ytcfg[
        "ramses-hydro", "fields"
    ] = """
    Density
    x-velocity
    y-velocity
    z-velocity
    Foo
    Bar
    """
    ds = yt.load(ramsesCosmo)

    def check_unit(array, unit):
        assert str(array.in_cgs().units) == unit

    try:
        assert ("ramses", "Foo") in ds.derived_field_list
        assert ("ramses", "Bar") in ds.derived_field_list

        check_unit(ds.r["ramses", "Foo"], "dimensionless")
        check_unit(ds.r["ramses", "Bar"], "dimensionless")
    finally:
        ytcfg.remove_section("ramses-hydro")


@requires_file(output_00080)
def test_grav_detection():
    ds = yt.load(output_00080)

    # Test detection
    for k in "xyz":
        assert ("gravity", "%s-acceleration" % k) in ds.field_list
        assert ("gas", "acceleration_%s" % k) in ds.derived_field_list

    # Test access
    for k in "xyz":
        ds.r["gas", "acceleration_%s" % k]


@requires_file(ramses_sink)
@requires_file(output_00080)
def test_ramses_field_detection():
    ds1 = yt.load(ramses_rt)

    assert "ramses" not in DETECTED_FIELDS

    # Now they are detected !
    ds1.index
    P1 = HydroFieldFileHandler.parameters
    assert DETECTED_FIELDS[ds1.unique_identifier]["ramses"] is not None
    fields_1 = set(DETECTED_FIELDS[ds1.unique_identifier]["ramses"])

    # Check the right number of variables has been loaded
    assert P1["nvar"] == 10
    assert len(fields_1) == P1["nvar"]

    # Now load another dataset
    ds2 = yt.load(output_00080)
    ds2.index
    P2 = HydroFieldFileHandler.parameters
    fields_2 = set(DETECTED_FIELDS[ds2.unique_identifier]["ramses"])

    # Check the right number of variables has been loaded
    assert P2["nvar"] == 6
    assert len(fields_2) == P2["nvar"]


@requires_file(ramses_new_format)
@requires_file(ramsesCosmo)
@requires_file(ramsesNonCosmo)
def test_formation_time():
    extra_particle_fields = [
        ("particle_birth_time", "d"),
        ("particle_metallicity", "d"),
    ]

    # test semantics for cosmological dataset
    ds = yt.load(ramsesCosmo, extra_particle_fields=extra_particle_fields)
    assert ("io", "particle_birth_time") in ds.field_list
    assert ("io", "conformal_birth_time") in ds.field_list
    assert ("io", "particle_metallicity") in ds.field_list

    ad = ds.all_data()
    whstars = ad["conformal_birth_time"] != 0
    assert np.all(ad["star_age"][whstars] > 0)

    # test semantics for non-cosmological new-style output format
    ds = yt.load(ramses_new_format)
    ad = ds.all_data()
    assert ("io", "particle_birth_time") in ds.field_list
    assert np.all(ad["particle_birth_time"] > 0)

    # test semantics for non-cosmological old-style output format
    ds = yt.load(ramsesNonCosmo, extra_particle_fields=extra_particle_fields)
    ad = ds.all_data()
    assert ("io", "particle_birth_time") in ds.field_list
    # the dataset only includes particles with arbitrarily old ages
    # and particles that formed in the very first timestep
    assert np.all(ad["particle_birth_time"] <= 0)
    whdynstars = ad["particle_birth_time"] == 0
    assert np.all(ad["star_age"][whdynstars] == ds.current_time)


@requires_file(ramses_new_format)
def test_cooling_fields():

    # Test the field is being loaded correctly
    ds = yt.load(ramses_new_format)

    # Derived cooling fields
    assert ("gas", "cooling_net") in ds.derived_field_list
    assert ("gas", "cooling_total") in ds.derived_field_list
    assert ("gas", "heating_total") in ds.derived_field_list
    assert ("gas", "number_density") in ds.derived_field_list

    # Original cooling fields
    assert ("gas", "cooling_primordial") in ds.derived_field_list
    assert ("gas", "cooling_compton") in ds.derived_field_list
    assert ("gas", "cooling_metal") in ds.derived_field_list
    assert ("gas", "heating_primordial") in ds.derived_field_list
    assert ("gas", "heating_compton") in ds.derived_field_list
    assert ("gas", "cooling_primordial_prime") in ds.derived_field_list
    assert ("gas", "cooling_compton_prime") in ds.derived_field_list
    assert ("gas", "cooling_metal_prime") in ds.derived_field_list
    assert ("gas", "heating_primordial_prime") in ds.derived_field_list
    assert ("gas", "heating_compton_prime") in ds.derived_field_list
    assert ("gas", "mu") in ds.derived_field_list

    # Abundances
    assert ("gas", "Electron_number_density") in ds.derived_field_list
    assert ("gas", "HI_number_density") in ds.derived_field_list
    assert ("gas", "HII_number_density") in ds.derived_field_list
    assert ("gas", "HeI_number_density") in ds.derived_field_list
    assert ("gas", "HeII_number_density") in ds.derived_field_list
    assert ("gas", "HeIII_number_density") in ds.derived_field_list

    def check_unit(array, unit):
        assert str(array.in_cgs().units) == unit

    check_unit(ds.r[("gas", "cooling_total")], "cm**5*g/s**3")
    check_unit(ds.r[("gas", "cooling_primordial_prime")], "cm**5*g/(K*s**3)")
    check_unit(ds.r[("gas", "number_density")], "cm**(-3)")
    check_unit(ds.r[("gas", "mu")], "dimensionless")
    check_unit(ds.r[("gas", "Electron_number_density")], "cm**(-3)")


ramses_rt = "ramses_rt_00088/output_00088/info_00088.txt"


@requires_file(ramses_rt)
def test_ramses_mixed_files():
    # Test that one can use derived fields that depend on different
    # files (here hydro and rt files)
    ds = yt.load(ramses_rt)

    def _mixed_field(field, data):
        return data["rt", "photon_density_1"] / data["gas", "H_nuclei_density"]

    ds.add_field(("gas", "mixed_files"), function=_mixed_field, sampling_type="cell")

    # Access the field
    ds.r[("gas", "mixed_files")]


ramses_empty_record = "ramses_empty_record/output_00003/info_00003.txt"


@requires_file(ramses_empty_record)
def test_ramses_empty_record():
    # Test that yt can load datasets with empty records
    ds = yt.load(ramses_empty_record)

    # This should not fail
    ds.index

    # Access some field
    ds.r[("gas", "density")]


@requires_file(ramses_new_format)
@requires_module("f90nml")
def test_namelist_reading():
    ds = data_dir_load(ramses_new_format)
    namelist_fname = os.path.join(ds.directory, "namelist.txt")
    with open(namelist_fname, "r") as f:
        ref = f90nml.read(f)

    nml = ds.parameters["namelist"]

    assert nml == ref


@requires_file(ramses_empty_record)
@requires_file(output_00080)
@requires_module("f90nml")
def test_namelist_reading_should_not_fail():

    for ds_name in (ramses_empty_record, output_00080):
        # Test that the reading does not fail for malformed namelist.txt files
        ds = data_dir_load(ds_name)
        ds.index  # should work


ramses_mhd_128 = "ramses_mhd_128/output_00027/info_00027.txt"


@requires_file(ramses_mhd_128)
def test_magnetic_field_aliasing():
    # Test if RAMSES magnetic fields are correctly aliased to yt magnetic fields and if
    # derived magnetic quantities are calculated
    ds = data_dir_load(ramses_mhd_128)
    ad = ds.all_data()

    for field in [
        "magnetic_field_x",
        "magnetic_field_magnitude",
        "alfven_speed",
        "magnetic_field_divergence",
    ]:
        assert ("gas", field) in ds.derived_field_list
        ad[("gas", field)]
