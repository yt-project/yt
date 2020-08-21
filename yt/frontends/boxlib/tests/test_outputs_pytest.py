"""
Title: test_boxlib.py
Purpose: Boxlib frontend tests
Notes:
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.frontends.boxlib.api import (
    CastroDataset,
    MaestroDataset,
    NyxDataset,
    OrionDataset,
    WarpXDataset,
)
from yt.testing import assert_equal, units_override_check
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    grid_hierarchy,
    grid_values,
    parentage_relationships,
    projection_values,
)

# Test data
radadvect = "RadAdvect/plt00000"
rt = "RadTube/plt00500"
star = "StarParticles/plrd01000"
LyA = "Nyx_LyA/plt00000"
RT_particles = "RT_particles/plt00050"
langmuir = "LangmuirWave/plt00020_v2"
plasma = "PlasmaAcceleration/plt00030_v2"
beam = "GaussianBeam/plt03008"
raw_fields = "Laser/plt00015"
nyx_no_particles = "nyx_sedov_plt00086"
msubch = "maestro_subCh_plt00248"


ds_list = [
    radadvect,
    rt,
    star,
    LyA,
    RT_particles,
    langmuir,
    plasma,
    beam,
]
a_list = [0, 1, 2]

orion_fields = ["temperature", "density", "velocity_magnitude"]
nyx_fields = ["Ne", "Temp", "particle_mass_density"]
castro_fields = ["Temp", "density", "particle_count"]
warpx_fields = ["Ex", "By", "jz"]
raw_fields = [("raw", "Bx"), ("raw", "Ey"), ("raw", "jz")]

w_radadvect = [None, "density"]
w_rt = [None, "density"]
w_star = [None, "density"]
w_lya = [None, "Ne"]
w_RTP = [None, "density"]
w_warp = [None, "Ex"]

d_radadvect = [None, ("sphere", ("max", (0.1, "unitary")))]
d_rt = [None, ("sphere", ("max", (0.1, "unitary")))]
d_star = [None, ("sphere", ("max", (0.1, "unitary")))]
d_lya = [None, ("sphere", ("c", (0.1, "unitary")))]
d_RTP = [None, ("sphere", ("max", (0.1, "unitary")))]
d_warp = [None, ("sphere", ("c", (0.1, "unitary")))]

pair_list = [
    [radadvect, orion_fields, w_radadvect, d_radadvect],
    [rt, orion_fields, w_rt, d_rt],
    [star, orion_fields, w_star, d_star],
    [LyA, nyx_fields, w_lya, d_lya],
    [RT_particles, castro_fields, w_RTP, d_RTP],
    [langmuir, warpx_fields, w_warp, d_warp],
    [plasma, warpx_fields, w_warp, d_warp],
    [beam, warpx_fields, w_warp, d_warp],
]

gv_pairs = [(i[0], f) for i in pair_list for f in i[1]]
fv_pairs = [(i[0], f, d) for i in pair_list for f in i[1] for d in i[3]]
pv_pairs = [
    (i[0], f, d, w) for i in pair_list for f in i[1] for d in i[3] for w in i[2]
]


@pytest.mark.answer_test
class TestBoxLib:
    answer_file = None

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", ds_list, indirect=True)
    def test_gh_pr(self, ds):
        self.hashes.update({"grid_hierarchy": grid_hierarchy(ds)})
        self.hashes.update({"parentage_relationships": parentage_relationships(ds)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f", gv_pairs, indirect=True)
    def test_gv(self, f, ds):
        self.hashes.update({"grid_values": grid_values(ds, f)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d", fv_pairs, indirect=True)
    def test_fv(self, d, f, ds):
        self.hashes.update({"field_values": field_values(ds, f, d)})

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds, f, d, w", pv_pairs, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    def test_pv(self, a, d, w, f, ds):
        self.hashes.update({"projection_values": projection_values(ds, a, f, w, d)})

    @pytest.mark.parametrize("ds", [LyA], indirect=True)
    def test_nyx_particle_io(self, ds):
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

    @pytest.mark.parametrize("ds", [RT_particles], indirect=True)
    def test_castro_particle_io(self, ds):
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

    @pytest.mark.parametrize("ds", [plasma], indirect=True)
    def test_warpx_particle_io(self, ds):
        grid = ds.index.grids[0]
        # read directly from the header
        npart0_grid_0 = 344
        npart1_grid_0 = 69632
        assert_equal(grid["particle0", "particle_position_x"].size, npart0_grid_0)
        assert_equal(grid["particle1", "particle_position_y"].size, npart1_grid_0)
        assert_equal(
            grid["all", "particle_position_z"].size, npart0_grid_0 + npart1_grid_0
        )
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

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [raw_fields], indirect=True)
    @pytest.mark.parametrize("f", raw_fields, indirect=True)
    def test_raw_fields(self, f, ds):
        gv = grid_values(ds, f)
        self.hashes.update({"grid_values": gv})

    @pytest.mark.parametrize("ds", [rt], indirect=True)
    def test_OrionDataset(self, ds):
        assert isinstance(ds, OrionDataset)

    @pytest.mark.parametrize("ds", [LyA], indirect=True)
    def test_NyxDataset(self, ds):
        assert isinstance(ds, NyxDataset)

    @pytest.mark.parametrize("ds", [RT_particles], indirect=True)
    def test_CastroDataset(self, ds):
        assert isinstance(ds, CastroDataset)

    @pytest.mark.parametrize("ds", [plasma], indirect=True)
    def test_WarpXDataset(self, ds):
        assert isinstance(ds, WarpXDataset)

    def test_units_override(self):
        units_override_check(rt)

    @pytest.mark.parametrize("ds", [nyx_no_particles], indirect=True)
    def test_nyx_no_part(self, ds):
        assert isinstance(ds, NyxDataset)
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
        assert_equal(sorted(ds.field_list), fields)

    @pytest.mark.parametrize("ds", [msubch], indirect=True)
    def test_maestro_parameters(self, ds):
        assert isinstance(ds, MaestroDataset)
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
