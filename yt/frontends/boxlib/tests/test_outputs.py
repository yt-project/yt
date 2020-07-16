import numpy as np

from yt.frontends.boxlib.api import (
    CastroDataset,
    MaestroDataset,
    NyxDataset,
    OrionDataset,
    WarpXDataset,
)
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    GridValuesTest,
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

# We don't do anything needing ghost zone generation right now, because these
# are non-periodic datasets.
_orion_fields = ("temperature", "density", "velocity_magnitude")
_nyx_fields = ("Ne", "Temp", "particle_mass_density")
_warpx_fields = ("Ex", "By", "jz")
_castro_fields = ("Temp", "density", "particle_count")

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
    for test in small_patch_amr(ds, _nyx_fields, input_center="c", input_weight="Ne"):
        test_LyA.__name__ = test.description
        yield test


@requires_file(LyA)
def test_nyx_particle_io():
    ds = data_dir_load(LyA)

    grid = ds.index.grids[0]
    npart_grid_0 = 7908  # read directly from the header
    assert_equal(grid["particle_position_x"].size, npart_grid_0)
    assert_equal(grid["DM", "particle_position_y"].size, npart_grid_0)
    assert_equal(grid["all", "particle_position_z"].size, npart_grid_0)

    ad = ds.all_data()
    npart = 32768  # read directly from the header
    assert_equal(ad["particle_velocity_x"].size, npart)
    assert_equal(ad["DM", "particle_velocity_y"].size, npart)
    assert_equal(ad["all", "particle_velocity_z"].size, npart)

    assert np.all(ad["particle_mass"] == ad["particle_mass"][0])

    left_edge = ds.arr([0.0, 0.0, 0.0], "code_length")
    right_edge = ds.arr([4.0, 4.0, 4.0], "code_length")
    center = 0.5 * (left_edge + right_edge)

    reg = ds.region(center, left_edge, right_edge)

    assert np.all(
        np.logical_and(
            reg["particle_position_x"] <= right_edge[0],
            reg["particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["particle_position_y"] <= right_edge[1],
            reg["particle_position_y"] >= left_edge[1],
        )
    )

    assert np.all(
        np.logical_and(
            reg["particle_position_z"] <= right_edge[2],
            reg["particle_position_z"] >= left_edge[2],
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
    assert_equal(grid["particle_position_x"].size, npart_grid_2)
    assert_equal(grid["Tracer", "particle_position_y"].size, npart_grid_2)
    assert_equal(grid["all", "particle_position_y"].size, npart_grid_2)

    ad = ds.all_data()
    npart = 49  # read directly from the header
    assert_equal(ad["particle_velocity_x"].size, npart)
    assert_equal(ad["Tracer", "particle_velocity_y"].size, npart)
    assert_equal(ad["all", "particle_velocity_y"].size, npart)

    left_edge = ds.arr([0.0, 0.0, 0.0], "code_length")
    right_edge = ds.arr([0.25, 1.0, 1.0], "code_length")
    center = 0.5 * (left_edge + right_edge)

    reg = ds.region(center, left_edge, right_edge)

    assert np.all(
        np.logical_and(
            reg["particle_position_x"] <= right_edge[0],
            reg["particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["particle_position_y"] <= right_edge[1],
            reg["particle_position_y"] >= left_edge[1],
        )
    )


langmuir = "LangmuirWave/plt00020_v2"


@requires_ds(langmuir)
def test_langmuir():
    ds = data_dir_load(langmuir)
    assert_equal(str(ds), "plt00020_v2")
    for test in small_patch_amr(ds, _warpx_fields, input_center="c", input_weight="Ex"):
        test_langmuir.__name__ = test.description
        yield test


plasma = "PlasmaAcceleration/plt00030_v2"


@requires_ds(plasma)
def test_plasma():
    ds = data_dir_load(plasma)
    assert_equal(str(ds), "plt00030_v2")
    for test in small_patch_amr(ds, _warpx_fields, input_center="c", input_weight="Ex"):
        test_plasma.__name__ = test.description
        yield test


beam = "GaussianBeam/plt03008"


@requires_ds(beam)
def test_beam():
    ds = data_dir_load(beam)
    assert_equal(str(ds), "plt03008")
    for test in small_patch_amr(ds, _warpx_fields, input_center="c", input_weight="Ex"):
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
            reg["particle_position_x"] <= right_edge[0],
            reg["particle_position_x"] >= left_edge[0],
        )
    )

    assert np.all(
        np.logical_and(
            reg["particle_position_y"] <= right_edge[1],
            reg["particle_position_y"] >= left_edge[1],
        )
    )

    assert np.all(
        np.logical_and(
            reg["particle_position_z"] <= right_edge[2],
            reg["particle_position_z"] >= left_edge[2],
        )
    )


_raw_fields = [("raw", "Bx"), ("raw", "Ey"), ("raw", "jz")]

raw_fields = "Laser/plt00015"


@requires_ds(raw_fields)
def test_raw_fields():
    ds_fn = raw_fields
    for field in _raw_fields:
        yield GridValuesTest(ds_fn, field)


@requires_file(rt)
def test_OrionDataset():
    assert isinstance(data_dir_load(rt), OrionDataset)


@requires_file(LyA)
def test_NyxDataset():
    assert isinstance(data_dir_load(LyA), NyxDataset)


@requires_file(RT_particles)
def test_CastroDataset():
    assert isinstance(data_dir_load(RT_particles), CastroDataset)


@requires_file(LyA)
def test_WarpXDataset():
    assert isinstance(data_dir_load(plasma), WarpXDataset)


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
    assert type(ds.parameters["plot_base_name"]) is str

    # Check boolean parameters: T or F
    assert not ds.parameters["use_thermal_diffusion"]
    assert type(ds.parameters["use_thermal_diffusion"]) is bool

    assert ds.parameters["do_burning"]
    assert type(ds.parameters["do_burning"]) is bool

    # Check a float parameter with a decimal point
    assert ds.parameters["sponge_kappa"] == float("10.00000000")
    assert type(ds.parameters["sponge_kappa"]) is float

    # Check a float parameter with E exponent notation
    assert ds.parameters["small_dt"] == float("0.1000000000E-09")

    # Check an int parameter
    assert ds.parameters["s0_interp_type"] == 3
    assert type(ds.parameters["s0_interp_type"]) is int
