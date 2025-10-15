import numpy as np
from numpy.testing import assert_allclose, assert_equal

from yt.frontends.amrex.api import (
    AMReXDataset,
    CastroDataset,
    MaestroDataset,
    NyxDataset,
    OrionDataset,
    QuokkaDataset,
    WarpXDataset,
)
from yt.loaders import load
from yt.testing import (
    disable_dataset_cache,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    GridValuesTest,
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

# We don't do anything needing ghost zone generation right now, because these
# are non-periodic datasets.
_orion_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
)
_nyx_fields = (
    ("boxlib", "Ne"),
    ("boxlib", "Temp"),
    ("boxlib", "particle_mass_density"),
)
_warpx_fields = (("mesh", "Ex"), ("mesh", "By"), ("mesh", "jz"))
_castro_fields = (
    ("boxlib", "Temp"),
    ("gas", "density"),
    ("boxlib", "particle_count"),
)

radadvect = "RadAdvect/plt00000"


@requires_ds(radadvect)
def test_radadvect():
    ds = data_dir_load(radadvect)
    assert_equal(str(ds), "plt00000")
    for test in small_patch_amr(ds, _orion_fields):
        test_radadvect.__name__ = test.description
        yield test


rt = "RadTube/plt00500"


@requires_ds(rt)
def test_radtube():
    ds = data_dir_load(rt)
    assert_equal(str(ds), "plt00500")
    for test in small_patch_amr(ds, _orion_fields):
        test_radtube.__name__ = test.description
        yield test


star = "StarParticles/plrd01000"


@requires_ds(star)
def test_star():
    ds = data_dir_load(star)
    assert_equal(str(ds), "plrd01000")
    for test in small_patch_amr(ds, _orion_fields):
        test_star.__name__ = test.description
        yield test


LyA = "Nyx_LyA/plt00000"


@requires_ds(LyA)
def test_LyA():
    ds = data_dir_load(LyA)
    assert_equal(str(ds), "plt00000")
    for test in small_patch_amr(
        ds, _nyx_fields, input_center="c", input_weight=("boxlib", "Ne")
    ):
        test_LyA.__name__ = test.description
        yield test


@requires_file(LyA)
def test_nyx_particle_io():
    ds = data_dir_load(LyA)

    grid = ds.index.grids[0]
    npart_grid_0 = 7908  # read directly from the header
    assert_equal(grid["all", "particle_position_x"].size, npart_grid_0)
    assert_equal(grid["DM", "particle_position_y"].size, npart_grid_0)
    assert_equal(grid["all", "particle_position_z"].size, npart_grid_0)

    ad = ds.all_data()
    npart = 32768  # read directly from the header
    assert_equal(ad["all", "particle_velocity_x"].size, npart)
    assert_equal(ad["DM", "particle_velocity_y"].size, npart)
    assert_equal(ad["all", "particle_velocity_z"].size, npart)

    assert np.all(ad["all", "particle_mass"] == ad["all", "particle_mass"][0])

    left_edge = ds.arr([0.0, 0.0, 0.0], "code_length")
    right_edge = ds.arr([4.0, 4.0, 4.0], "code_length")
    center = 0.5 * (left_edge + right_edge)

    reg = ds.region(center, left_edge, right_edge)

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_x"] <= right_edge[0],
            reg["all", "particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_y"] <= right_edge[1],
            reg["all", "particle_position_y"] >= left_edge[1],
        )
    )

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_z"] <= right_edge[2],
            reg["all", "particle_position_z"] >= left_edge[2],
        )
    )


RT_particles = "RT_particles/plt00050"


@requires_ds(RT_particles)
def test_RT_particles():
    ds = data_dir_load(RT_particles)
    assert_equal(str(ds), "plt00050")
    for test in small_patch_amr(ds, _castro_fields):
        test_RT_particles.__name__ = test.description
        yield test


@requires_file(RT_particles)
def test_castro_particle_io():
    ds = data_dir_load(RT_particles)

    grid = ds.index.grids[2]
    npart_grid_2 = 49  # read directly from the header
    assert_equal(grid["all", "particle_position_x"].size, npart_grid_2)
    assert_equal(grid["Tracer", "particle_position_y"].size, npart_grid_2)
    assert_equal(grid["all", "particle_position_y"].size, npart_grid_2)

    ad = ds.all_data()
    npart = 49  # read directly from the header
    assert_equal(ad["all", "particle_velocity_x"].size, npart)
    assert_equal(ad["Tracer", "particle_velocity_y"].size, npart)
    assert_equal(ad["all", "particle_velocity_y"].size, npart)

    left_edge = ds.arr([0.0, 0.0, 0.0], "code_length")
    right_edge = ds.arr([0.25, 1.0, 1.0], "code_length")
    center = 0.5 * (left_edge + right_edge)

    reg = ds.region(center, left_edge, right_edge)

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_x"] <= right_edge[0],
            reg["all", "particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_y"] <= right_edge[1],
            reg["all", "particle_position_y"] >= left_edge[1],
        )
    )


langmuir = "LangmuirWave/plt00020_v2"


@requires_ds(langmuir)
def test_langmuir():
    ds = data_dir_load(langmuir)
    assert_equal(str(ds), "plt00020_v2")
    for test in small_patch_amr(
        ds, _warpx_fields, input_center="c", input_weight=("mesh", "Ex")
    ):
        test_langmuir.__name__ = test.description
        yield test


plasma = "PlasmaAcceleration/plt00030_v2"


@requires_ds(plasma)
def test_plasma():
    ds = data_dir_load(plasma)
    assert_equal(str(ds), "plt00030_v2")
    for test in small_patch_amr(
        ds, _warpx_fields, input_center="c", input_weight=("mesh", "Ex")
    ):
        test_plasma.__name__ = test.description
        yield test


beam = "GaussianBeam/plt03008"


@requires_ds(beam)
def test_beam():
    ds = data_dir_load(beam)
    assert_equal(str(ds), "plt03008")
    for param in ("number of boxes", "maximum zones"):
        # PR 2807
        # these parameters are only populated if the config file attached to this
        # dataset is read correctly
        assert param in ds.parameters
    for test in small_patch_amr(
        ds, _warpx_fields, input_center="c", input_weight=("mesh", "Ex")
    ):
        test_beam.__name__ = test.description
        yield test


@requires_file(plasma)
def test_warpx_particle_io():
    ds = data_dir_load(plasma)
    grid = ds.index.grids[0]

    # read directly from the header
    npart0_grid_0 = 344
    npart1_grid_0 = 69632

    assert_equal(grid["particle0", "particle_position_x"].size, npart0_grid_0)
    assert_equal(grid["particle1", "particle_position_y"].size, npart1_grid_0)
    assert_equal(grid["all", "particle_position_z"].size, npart0_grid_0 + npart1_grid_0)

    # read directly from the header
    npart0 = 1360
    npart1 = 802816
    ad = ds.all_data()
    assert_equal(ad["particle0", "particle_velocity_x"].size, npart0)
    assert_equal(ad["particle1", "particle_velocity_y"].size, npart1)
    assert_equal(ad["all", "particle_velocity_z"].size, npart0 + npart1)

    np.all(ad["particle1", "particle_mass"] == ad["particle1", "particle_mass"][0])
    np.all(ad["particle0", "particle_mass"] == ad["particle0", "particle_mass"][0])

    left_edge = ds.arr([-7.5e-5, -7.5e-5, -7.5e-5], "code_length")
    right_edge = ds.arr([2.5e-5, 2.5e-5, 2.5e-5], "code_length")
    center = 0.5 * (left_edge + right_edge)

    reg = ds.region(center, left_edge, right_edge)

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_x"] <= right_edge[0],
            reg["all", "particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_y"] <= right_edge[1],
            reg["all", "particle_position_y"] >= left_edge[1],
        )
    )

    assert np.all(
        np.logical_and(
            reg["all", "particle_position_z"] <= right_edge[2],
            reg["all", "particle_position_z"] >= left_edge[2],
        )
    )


_raw_fields = [("raw", "Bx"), ("raw", "Ey"), ("raw", "jz")]

laser = "Laser/plt00015"


@requires_ds(laser)
def test_raw_fields():
    for field in _raw_fields:
        yield GridValuesTest(laser, field)


@requires_file(rt)
def test_OrionDataset():
    assert isinstance(data_dir_load(rt), OrionDataset)


@requires_file(LyA)
def test_NyxDataset():
    assert isinstance(data_dir_load(LyA), NyxDataset)


@requires_file("nyx_small/nyx_small_00000")
def test_NyxDataset_2():
    assert isinstance(data_dir_load("nyx_small/nyx_small_00000"), NyxDataset)


@requires_file(RT_particles)
def test_CastroDataset():
    assert isinstance(data_dir_load(RT_particles), CastroDataset)


@requires_file("castro_sod_x_plt00036")
def test_CastroDataset_2():
    assert isinstance(data_dir_load("castro_sod_x_plt00036"), CastroDataset)


@requires_file("castro_sedov_1d_cyl_plt00150")
def test_CastroDataset_3():
    assert isinstance(data_dir_load("castro_sedov_1d_cyl_plt00150"), CastroDataset)


@requires_file(plasma)
def test_WarpXDataset():
    assert isinstance(data_dir_load(plasma), WarpXDataset)


@disable_dataset_cache
@requires_file(plasma)
def test_magnetic_units():
    ds1 = load(plasma)
    assert_allclose(ds1.magnetic_unit.value, 1.0)
    assert str(ds1.magnetic_unit.units) == "T"
    mag_unit1 = ds1.magnetic_unit.to("code_magnetic")
    assert_allclose(mag_unit1.value, 1.0)
    assert str(mag_unit1.units) == "code_magnetic"
    ds2 = load(plasma, unit_system="cgs")
    assert_allclose(ds2.magnetic_unit.value, 1.0e4)
    assert str(ds2.magnetic_unit.units) == "G"
    mag_unit2 = ds2.magnetic_unit.to("code_magnetic")
    assert_allclose(mag_unit2.value, 1.0)
    assert str(mag_unit2.units) == "code_magnetic"


@requires_ds(laser)
def test_WarpXDataset_2():
    assert isinstance(data_dir_load(laser), WarpXDataset)


@requires_file("plt.Cavity00010")
def test_AMReXDataset():
    ds = data_dir_load("plt.Cavity00010", kwargs={"cparam_filename": "inputs"})
    assert isinstance(ds, AMReXDataset)


@requires_file(rt)
def test_units_override():
    units_override_check(rt)


nyx_no_particles = "nyx_sedov_plt00086"


@requires_file(nyx_no_particles)
def test_nyx_no_part():
    assert isinstance(data_dir_load(nyx_no_particles), NyxDataset)

    fields = sorted(
        [
            ("boxlib", "H"),
            ("boxlib", "He"),
            ("boxlib", "MachNumber"),
            ("boxlib", "Ne"),
            ("boxlib", "Rank"),
            ("boxlib", "StateErr"),
            ("boxlib", "Temp"),
            ("boxlib", "X(H)"),
            ("boxlib", "X(He)"),
            ("boxlib", "density"),
            ("boxlib", "divu"),
            ("boxlib", "eint_E"),
            ("boxlib", "eint_e"),
            ("boxlib", "entropy"),
            ("boxlib", "forcex"),
            ("boxlib", "forcey"),
            ("boxlib", "forcez"),
            ("boxlib", "kineng"),
            ("boxlib", "logden"),
            ("boxlib", "magmom"),
            ("boxlib", "magvel"),
            ("boxlib", "magvort"),
            ("boxlib", "pressure"),
            ("boxlib", "rho_E"),
            ("boxlib", "rho_H"),
            ("boxlib", "rho_He"),
            ("boxlib", "rho_e"),
            ("boxlib", "soundspeed"),
            ("boxlib", "x_velocity"),
            ("boxlib", "xmom"),
            ("boxlib", "y_velocity"),
            ("boxlib", "ymom"),
            ("boxlib", "z_velocity"),
            ("boxlib", "zmom"),
        ]
    )

    ds = data_dir_load(nyx_no_particles)
    assert_equal(sorted(ds.field_list), fields)


msubch = "maestro_subCh_plt00248"


@requires_file(msubch)
def test_maestro_parameters():
    assert isinstance(data_dir_load(msubch), MaestroDataset)
    ds = data_dir_load(msubch)

    # Check a string parameter
    assert ds.parameters["plot_base_name"] == "subCh_hot_baserun_plt"
    assert type(ds.parameters["plot_base_name"]) is str  # noqa: E721

    # Check boolean parameters: T or F
    assert not ds.parameters["use_thermal_diffusion"]
    assert type(ds.parameters["use_thermal_diffusion"]) is bool  # noqa: E721

    assert ds.parameters["do_burning"]
    assert type(ds.parameters["do_burning"]) is bool  # noqa: E721

    # Check a float parameter with a decimal point
    assert ds.parameters["sponge_kappa"] == float("10.00000000")
    assert type(ds.parameters["sponge_kappa"]) is float  # noqa: E721

    # Check a float parameter with E exponent notation
    assert ds.parameters["small_dt"] == float("0.1000000000E-09")

    # Check an int parameter
    assert ds.parameters["s0_interp_type"] == 3
    assert type(ds.parameters["s0_interp_type"]) is int  # noqa: E721


quokka = "quokka_RadiatingParticles_plt026"


@requires_file(quokka)
def test_quokka():
    # Test QuokkaDataset instance type
    ds = data_dir_load(quokka)
    assert isinstance(ds, QuokkaDataset)

    # Test basic parameters
    assert ds.parameters["HydroMethod"] == "Quokka"
    assert ds.parameters["quokka_version"] == 25.03
    assert ds.parameters["dimensionality"] == 2
    assert_equal(
        ds.parameters["fields"],
        [
            "gasDensity",
            "x-GasMomentum",
            "y-GasMomentum",
            "z-GasMomentum",
            "gasEnergy",
            "gasInternalEnergy",
            "radEnergy-Group0",
            "x-RadFlux-Group0",
            "y-RadFlux-Group0",
            "z-RadFlux-Group0",
        ],
    )

    # Check domain dimensions and boundaries
    assert_equal(ds.domain_dimensions, [64, 64, 1])
    assert_equal(ds.domain_left_edge, [0.0, 0.0, 0.0])
    assert_equal(ds.domain_right_edge, [1.0, 1.0, 1.0])

    # Test derived fields and unit conversions
    ad = ds.all_data()
    assert ("gas", "density") in ds.derived_field_list
    density = ad["gas", "density"]
    assert str(density.units) == "g/cm**3"

    # Check momentum/velocity relationship which confirms velocity is correctly derived
    momentum_x = ad["boxlib", "x-GasMomentum"]
    gas_density = ad["boxlib", "gasDensity"]
    velocity_x = ad["gas", "velocity_x"]
    assert_allclose(momentum_x / gas_density, velocity_x)

    # Test radiation fields specifically
    assert ds.parameters["radiation_field_groups"] == 1
    assert "rad" in ds.fluid_types
    assert ("rad", "energy_density_0") in ds.derived_field_list

    # Test particle IO with specific expectations
    assert ds.parameters["particle_types"] == ("Rad_particles",)

    # Verify particle fields and units
    particle_info = ds.parameters["particle_info"]["Rad_particles"]
    assert particle_info["num_particles"] == 2
    assert particle_info["num_fields"] == 3
    assert set(particle_info["fields"]) == {"luminosity", "birth_time", "end_time"}
    assert particle_info["units"]["luminosity"] == "M^1 L^2 T^-3"
    assert particle_info["units"]["birth_time"] == "T^1"
    assert particle_info["units"]["end_time"] == "T^1"


quokka_face_centered = "quokka_HydroWave_plt00004"


@requires_file(quokka_face_centered)
def test_quokka_face_centered():
    ds = data_dir_load(quokka_face_centered)

    # Check if face-centered datasets were loaded
    assert hasattr(ds, "ds_fc_x")
    assert hasattr(ds, "ds_fc_y")
    assert hasattr(ds, "ds_fc_z")

    # Check that the face-centered datasets have the correct attributes
    assert ds.ds_fc_x.fc_direction == "x"
    assert ds.ds_fc_y.fc_direction == "y"
    assert ds.ds_fc_z.fc_direction == "z"

    # Check that the face-centered datasets have a reference to the parent dataset
    assert ds.ds_fc_x.parent_ds is ds
    assert ds.ds_fc_y.parent_ds is ds
    assert ds.ds_fc_z.parent_ds is ds

    assert isinstance(ds, QuokkaDataset)
    assert isinstance(ds.ds_fc_x, QuokkaDataset)


# test loading non-Cartesian coordinate systems in different dimensionalities


def check_coordsys_data(ds):
    # check that level_dds is consistent with domain_width
    assert_allclose(
        ds.index.level_dds[0] * ds.domain_dimensions,
        ds.domain_width.to_value("code_length"),
        rtol=1e-12,
        atol=0.0,
    )

    # check that we get the expected number of data points when selecting the
    # entire domain
    expected_size = sum(np.count_nonzero(g.child_mask) for g in ds.index.grids)
    ad = ds.all_data()
    assert ad["boxlib", "Temp"].size == expected_size


cyl_1d = "castro_sedov_1d_cyl_plt00150"
cyl_2d = "castro_sedov_2d_sph_in_cyl_plt00130"
sph_1d = "sedov_1d_sph_plt00120"
sph_2d = "xrb_spherical_smallplt00010"


@requires_file(cyl_1d)
def test_coordsys_1d_cylindrical():
    ds = data_dir_load(cyl_1d)
    assert ds.geometry == "cylindrical"
    assert ds.dimensionality == 1
    check_coordsys_data(ds)


@requires_file(cyl_2d)
def test_coordsys_2d_cylindrical():
    ds = data_dir_load(cyl_2d)
    assert ds.geometry == "cylindrical"
    assert ds.dimensionality == 2
    check_coordsys_data(ds)


@requires_file(sph_1d)
def test_coordsys_1d_spherical():
    ds = data_dir_load(sph_1d)
    assert ds.geometry == "spherical"
    assert ds.dimensionality == 1
    check_coordsys_data(ds)


@requires_file(sph_2d)
def test_coordsys_2d_spherical():
    ds = data_dir_load(sph_2d)
    assert ds.geometry == "spherical"
    assert ds.dimensionality == 2
    check_coordsys_data(ds)
