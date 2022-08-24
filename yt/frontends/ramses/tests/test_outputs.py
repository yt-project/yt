import os

import numpy as np

import yt
from yt.config import ytcfg
from yt.fields.field_detector import FieldDetector
from yt.frontends.ramses.api import RAMSESDataset
from yt.frontends.ramses.field_handlers import DETECTED_FIELDS, HydroFieldFileHandler
from yt.testing import (
    assert_equal,
    assert_raises,
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

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
    ("gas", "pressure_gradient_magnitude"),
)

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
                for weight_field in [None, ("gas", "density")]:
                    yield PixelizedProjectionValuesTest(
                        output_00080, axis, field, weight_field, dobj_name
                    )
            yield FieldValuesTest(output_00080, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj[("index", "ones")].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)


@requires_file(output_00080)
def test_RAMSESDataset():
    assert isinstance(data_dir_load(output_00080), RAMSESDataset)


@requires_file(output_00080)
def test_RAMSES_alternative_load():
    # Test that we can load a RAMSES dataset by giving it the folder name,
    # the info file name or an amr file name
    base_dir, info_file_fname = os.path.split(output_00080)
    for fname in (base_dir, os.path.join(base_dir, "amr_00080.out00001")):
        assert isinstance(data_dir_load(fname), RAMSESDataset)


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

        # NOTE: these are the old test values, which used 3.08e24 as
        # the Mpc to cm conversion factor
        # expected_raw_time = 1.119216564055017 # in ramses unit
        # expected_time = 3.756241729312462e17 # in seconds

        expected_raw_time = 1.121279694787743  # in ramses unit
        assert_equal(ds.current_time.value, expected_raw_time)

        expected_time = 3.7631658742904595e17  # in seconds
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
    extra_fields = [f"particle_extra_field_{i + 1}" for i in range(2)]
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
                ("gas", f"{specie}_fraction"),
                ("gas", f"{specie}_density"),
                ("gas", f"{specie}_mass"),
            ]
        )

    for field in special_fields:
        assert field in ds.derived_field_list
        ad[field]

    # Test the creation of rt fields
    rt_fields = [("rt", "photon_density_1")] + [
        ("rt", f"photon_flux_{d}_1") for d in "xyz"
    ]
    for field in rt_fields:
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
        "particle_prop_0_4",
        "particle_prop_0_5",
        "particle_prop_0_6",
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
        assert ("sink", field) not in ds.field_list


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


@requires_file(output_00080)
@requires_file(ramses_sink)
def test_grav_detection():
    for path, has_potential in ((output_00080, False), (ramses_sink, True)):
        ds = yt.load(path)

        # Test detection
        for k in "xyz":
            assert ("gravity", f"{k}-acceleration") in ds.field_list
            assert ("gas", f"acceleration_{k}") in ds.derived_field_list

        if has_potential:
            assert ("gravity", "Potential") in ds.field_list
            assert ("gas", "potential") in ds.derived_field_list
            assert ("gas", "potential_energy") in ds.derived_field_list

        # Test access
        for k in "xyz":
            ds.r["gas", f"acceleration_{k}"].to("m/s**2")

        if has_potential:
            ds.r["gas", "potential"].to("m**2/s**2")
            ds.r["gas", "potential_energy"].to("kg*m**2/s**2")


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
    whstars = ad[("io", "conformal_birth_time")] != 0
    assert np.all(ad[("io", "star_age")][whstars] > 0)

    # test semantics for non-cosmological new-style output format
    ds = yt.load(ramses_new_format)
    ad = ds.all_data()
    assert ("io", "particle_birth_time") in ds.field_list
    assert np.all(ad[("io", "particle_birth_time")] > 0)

    # test semantics for non-cosmological old-style output format
    ds = yt.load(ramsesNonCosmo, extra_particle_fields=extra_particle_fields)
    ad = ds.all_data()
    assert ("io", "particle_birth_time") in ds.field_list
    # the dataset only includes particles with arbitrarily old ages
    # and particles that formed in the very first timestep
    assert np.all(ad[("io", "particle_birth_time")] <= 0)
    whdynstars = ad[("io", "particle_birth_time")] == 0
    assert np.all(ad[("io", "star_age")][whdynstars] == ds.current_time)


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
    with open(namelist_fname) as f:
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


@requires_file(output_00080)
def test_field_accession():
    ds = yt.load(output_00080)
    fields = [
        ("gas", "density"),  # basic ones
        ("gas", "pressure"),
        ("gas", "pressure_gradient_magnitude"),  # requires ghost zones
    ]
    # Check accessing gradient works for a variety of spatial domains
    for reg in (
        ds.all_data(),
        ds.sphere([0.1] * 3, 0.01),
        ds.sphere([0.5] * 3, 0.05),
        ds.box([0.1] * 3, [0.2] * 3),
    ):
        for field in fields:
            reg[field]


@requires_file(output_00080)
def test_ghost_zones():
    ds = yt.load(output_00080)

    def gen_dummy(ngz):
        def dummy(field, data):
            if not isinstance(data, FieldDetector):
                shape = data["gas", "mach_number"].shape[:3]
                np.testing.assert_equal(shape, (2 + 2 * ngz, 2 + 2 * ngz, 2 + 2 * ngz))
            return data["gas", "mach_number"]

        return dummy

    fields = []
    for ngz in (1, 2, 3):
        fname = ("gas", f"density_ghost_zone_{ngz}")
        ds.add_field(
            fname,
            gen_dummy(ngz),
            sampling_type="cell",
            units="",
            validators=[yt.ValidateSpatial(ghost_zones=ngz)],
        )
        fields.append(fname)

    box = ds.box([0, 0, 0], [0.1, 0.1, 0.1])
    for f in fields:
        print("Getting ", f)
        box[f]


@requires_file(output_00080)
def test_max_level():
    ds = yt.load(output_00080)

    assert any(ds.r["index", "grid_level"] > 2)

    # Should work
    for ds in (
        yt.load(output_00080, max_level=2, max_level_convention="yt"),
        yt.load(output_00080, max_level=8, max_level_convention="ramses"),
    ):
        assert all(ds.r["index", "grid_level"] <= 2)
        assert any(ds.r["index", "grid_level"] == 2)


@requires_file(output_00080)
def test_invalid_max_level():
    invalid_value_args = (
        (1, None),
        (1, "foo"),
        (1, "bar"),
        (-1, "yt"),
    )
    for lvl, convention in invalid_value_args:
        with assert_raises(ValueError):
            yt.load(output_00080, max_level=lvl, max_level_convention=convention)

    invalid_type_args = (
        (1.0, "yt"),  # not an int
        ("invalid", "yt"),
    )
    # Should fail with value errors
    for lvl, convention in invalid_type_args:
        with assert_raises(TypeError):
            yt.load(output_00080, max_level=lvl, max_level_convention=convention)


@requires_file(ramses_new_format)
def test_print_stats():
    ds = yt.load(ramses_new_format)

    # Should work
    ds.print_stats()

    # FIXME #3197: use `capsys` with pytest to make sure the print_stats function works as intended
