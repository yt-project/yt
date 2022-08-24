import os
import shutil
import tempfile
import unittest

import numpy as np

import yt
from yt.data_objects.particle_filters import add_particle_filter
from yt.data_objects.profiles import Profile1D, Profile2D, Profile3D, create_profile
from yt.testing import (
    assert_equal,
    assert_raises,
    assert_rel_equal,
    fake_random_ds,
    fake_sph_orientation_ds,
    requires_module,
)
from yt.utilities.exceptions import YTIllDefinedProfile, YTProfileDataShape
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot

_fields = ("density", "temperature", "dinosaurs", "tribbles")
_units = ("g/cm**3", "K", "dyne", "erg")


def test_profiles():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, units=_units)
    nv = ds.domain_dimensions.prod()
    dd = ds.all_data()
    rt, tt, dt = dd.quantities["TotalQuantity"](
        [("gas", "density"), ("gas", "temperature"), ("stream", "dinosaurs")]
    )

    e1, e2 = 0.9, 1.1
    for nb in [8, 16, 32, 64]:
        for input_units in ["mks", "cgs"]:
            (rmi, rma), (tmi, tma), (dmi, dma) = (
                getattr(ex, f"in_{input_units}")()
                for ex in dd.quantities["Extrema"](
                    [
                        ("gas", "density"),
                        ("gas", "temperature"),
                        ("stream", "dinosaurs"),
                    ]
                )
            )
            # We log all the fields or don't log 'em all.  No need to do them
            # individually.
            for lf in [True, False]:
                direct_profile = Profile1D(
                    dd,
                    ("gas", "density"),
                    nb,
                    rmi * e1,
                    rma * e2,
                    lf,
                    weight_field=None,
                )
                direct_profile.add_fields([("index", "ones"), ("gas", "temperature")])

                indirect_profile_s = create_profile(
                    dd,
                    ("gas", "density"),
                    [("index", "ones"), ("gas", "temperature")],
                    n_bins=nb,
                    extrema={("gas", "density"): (rmi * e1, rma * e2)},
                    logs={("gas", "density"): lf},
                    weight_field=None,
                )

                indirect_profile_t = create_profile(
                    dd,
                    ("gas", "density"),
                    [("index", "ones"), ("gas", "temperature")],
                    n_bins=nb,
                    extrema={("gas", "density"): (rmi * e1, rma * e2)},
                    logs={("gas", "density"): lf},
                    weight_field=None,
                )

                for p1d in [direct_profile, indirect_profile_s, indirect_profile_t]:
                    assert_equal(p1d["index", "ones"].sum(), nv)
                    assert_rel_equal(tt, p1d["gas", "temperature"].sum(), 7)

                p2d = Profile2D(
                    dd,
                    ("gas", "density"),
                    nb,
                    rmi * e1,
                    rma * e2,
                    lf,
                    ("gas", "temperature"),
                    nb,
                    tmi * e1,
                    tma * e2,
                    lf,
                    weight_field=None,
                )
                p2d.add_fields([("index", "ones"), ("gas", "temperature")])
                assert_equal(p2d["index", "ones"].sum(), nv)
                assert_rel_equal(tt, p2d["gas", "temperature"].sum(), 7)

                p3d = Profile3D(
                    dd,
                    ("gas", "density"),
                    nb,
                    rmi * e1,
                    rma * e2,
                    lf,
                    ("gas", "temperature"),
                    nb,
                    tmi * e1,
                    tma * e2,
                    lf,
                    ("stream", "dinosaurs"),
                    nb,
                    dmi * e1,
                    dma * e2,
                    lf,
                    weight_field=None,
                )
                p3d.add_fields([("index", "ones"), ("gas", "temperature")])
                assert_equal(p3d["index", "ones"].sum(), nv)
                assert_rel_equal(tt, p3d["gas", "temperature"].sum(), 7)

        p1d = Profile1D(dd, ("index", "x"), nb, 0.0, 1.0, False, weight_field=None)
        p1d.add_fields(("index", "ones"))
        av = nv / nb
        assert_equal(p1d["index", "ones"], np.ones(nb) * av)

        # We re-bin ones with a weight now
        p1d = Profile1D(
            dd, ("index", "x"), nb, 0.0, 1.0, False, weight_field=("gas", "temperature")
        )
        p1d.add_fields([("index", "ones")])
        assert_equal(p1d["index", "ones"], np.ones(nb))

        # Verify we can access "ones" after adding a new field
        # See issue 988
        p1d.add_fields([("gas", "density")])
        assert_equal(p1d["index", "ones"], np.ones(nb))

        p2d = Profile2D(
            dd,
            ("index", "x"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "y"),
            nb,
            0.0,
            1.0,
            False,
            weight_field=None,
        )
        p2d.add_fields(("index", "ones"))
        av = nv / nb**2
        assert_equal(p2d["index", "ones"], np.ones((nb, nb)) * av)

        # We re-bin ones with a weight now
        p2d = Profile2D(
            dd,
            ("index", "x"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "y"),
            nb,
            0.0,
            1.0,
            False,
            weight_field=("gas", "temperature"),
        )
        p2d.add_fields([("index", "ones")])
        assert_equal(p2d["index", "ones"], np.ones((nb, nb)))

        p3d = Profile3D(
            dd,
            ("index", "x"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "y"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "z"),
            nb,
            0.0,
            1.0,
            False,
            weight_field=None,
        )
        p3d.add_fields(("index", "ones"))
        av = nv / nb**3
        assert_equal(p3d["index", "ones"], np.ones((nb, nb, nb)) * av)

        # We re-bin ones with a weight now
        p3d = Profile3D(
            dd,
            ("index", "x"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "y"),
            nb,
            0.0,
            1.0,
            False,
            ("index", "z"),
            nb,
            0.0,
            1.0,
            False,
            weight_field=("gas", "temperature"),
        )
        p3d.add_fields([("index", "ones")])
        assert_equal(p3d["index", "ones"], np.ones((nb, nb, nb)))

        p2d = create_profile(
            dd,
            ("gas", "density"),
            ("gas", "temperature"),
            weight_field=("gas", "cell_mass"),
            extrema={("gas", "density"): (None, rma * e2)},
        )
        assert_equal(p2d.x_bins[0], rmi - np.spacing(rmi))
        assert_equal(p2d.x_bins[-1], rma * e2)
        assert str(ds.field_info["gas", "cell_mass"].units) == str(p2d.weight.units)

        p2d = create_profile(
            dd,
            ("gas", "density"),
            ("gas", "temperature"),
            weight_field=("gas", "cell_mass"),
            extrema={("gas", "density"): (rmi * e2, None)},
        )
        assert_equal(p2d.x_bins[0], rmi * e2)
        assert_equal(p2d.x_bins[-1], rma + np.spacing(rma))


extrema_s = {("all", "particle_position_x"): (0, 1)}
logs_s = {("all", "particle_position_x"): False}

extrema_t = {("all", "particle_position_x"): (0, 1)}
logs_t = {("all", "particle_position_x"): False}


def test_particle_profiles():
    for nproc in [1, 2, 4, 8]:
        ds = fake_random_ds(32, nprocs=nproc, particles=32**3)
        dd = ds.all_data()

        p1d = Profile1D(
            dd, ("all", "particle_position_x"), 128, 0.0, 1.0, False, weight_field=None
        )
        p1d.add_fields([("all", "particle_ones")])
        assert_equal(p1d[("all", "particle_ones")].sum(), 32**3)

        p1d = create_profile(
            dd,
            [("all", "particle_position_x")],
            [("all", "particle_ones")],
            weight_field=None,
            n_bins=128,
            extrema=extrema_s,
            logs=logs_s,
        )
        assert_equal(p1d[("all", "particle_ones")].sum(), 32**3)

        p1d = create_profile(
            dd,
            [("all", "particle_position_x")],
            [("all", "particle_ones")],
            weight_field=None,
            n_bins=128,
            extrema=extrema_t,
            logs=logs_t,
        )
        assert_equal(p1d[("all", "particle_ones")].sum(), 32**3)

        p2d = Profile2D(
            dd,
            ("all", "particle_position_x"),
            128,
            0.0,
            1.0,
            False,
            ("all", "particle_position_y"),
            128,
            0.0,
            1.0,
            False,
            weight_field=None,
        )
        p2d.add_fields([("all", "particle_ones")])
        assert_equal(p2d[("all", "particle_ones")].sum(), 32**3)

        p3d = Profile3D(
            dd,
            ("all", "particle_position_x"),
            128,
            0.0,
            1.0,
            False,
            ("all", "particle_position_y"),
            128,
            0.0,
            1.0,
            False,
            ("all", "particle_position_z"),
            128,
            0.0,
            1.0,
            False,
            weight_field=None,
        )
        p3d.add_fields([("all", "particle_ones")])
        assert_equal(p3d[("all", "particle_ones")].sum(), 32**3)


def test_mixed_particle_mesh_profiles():
    ds = fake_random_ds(32, particles=10)
    ad = ds.all_data()
    assert_raises(
        YTIllDefinedProfile,
        ProfilePlot,
        ad,
        ("index", "radius"),
        ("all", "particle_mass"),
    )
    assert_raises(
        YTIllDefinedProfile,
        ProfilePlot,
        ad,
        "radius",
        [("all", "particle_mass"), ("all", "particle_ones")],
    )
    assert_raises(
        YTIllDefinedProfile,
        ProfilePlot,
        ad,
        ("index", "radius"),
        [("all", "particle_mass"), ("index", "ones")],
    )
    assert_raises(
        YTIllDefinedProfile,
        ProfilePlot,
        ad,
        ("all", "particle_radius"),
        ("all", "particle_mass"),
        ("gas", "cell_mass"),
    )
    assert_raises(
        YTIllDefinedProfile,
        ProfilePlot,
        ad,
        ("index", "radius"),
        ("gas", "cell_mass"),
        ("all", "particle_ones"),
    )

    assert_raises(
        YTIllDefinedProfile,
        PhasePlot,
        ad,
        ("index", "radius"),
        ("all", "particle_mass"),
        ("gas", "velocity_x"),
    )
    assert_raises(
        YTIllDefinedProfile,
        PhasePlot,
        ad,
        ("all", "particle_radius"),
        ("all", "particle_mass"),
        ("gas", "cell_mass"),
    )
    assert_raises(
        YTIllDefinedProfile,
        PhasePlot,
        ad,
        ("index", "radius"),
        ("gas", "cell_mass"),
        ("all", "particle_ones"),
    )
    assert_raises(
        YTIllDefinedProfile,
        PhasePlot,
        ad,
        ("all", "particle_radius"),
        ("all", "particle_mass"),
        ("all", "particle_ones"),
    )


def test_particle_profile_negative_field():
    # see Issue #1340
    n_particles = int(1e4)

    ppx, ppy, ppz = np.random.normal(size=[3, n_particles])
    pvx, pvy, pvz = -np.ones((3, n_particles))

    data = {
        "particle_position_x": ppx,
        "particle_position_y": ppy,
        "particle_position_z": ppz,
        "particle_velocity_x": pvx,
        "particle_velocity_y": pvy,
        "particle_velocity_z": pvz,
    }

    bbox = 1.1 * np.array(
        [[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]]
    )
    ds = yt.load_particles(data, bbox=bbox)
    ad = ds.all_data()

    profile = yt.create_profile(
        ad,
        [("all", "particle_position_x"), ("all", "particle_position_y")],
        ("all", "particle_velocity_x"),
        logs={
            ("all", "particle_position_x"): True,
            ("all", "particle_position_y"): True,
            ("all", "particle_position_z"): True,
        },
        weight_field=None,
    )
    assert profile[("all", "particle_velocity_x")].min() < 0
    assert profile.x_bins.min() > 0
    assert profile.y_bins.min() > 0

    profile = yt.create_profile(
        ad,
        [("all", "particle_position_x"), ("all", "particle_position_y")],
        ("all", "particle_velocity_x"),
        weight_field=None,
    )
    assert profile[("all", "particle_velocity_x")].min() < 0
    assert profile.x_bins.min() < 0
    assert profile.y_bins.min() < 0

    # can't use CIC deposition with log-scaled bin fields
    with assert_raises(RuntimeError):
        yt.create_profile(
            ad,
            [("all", "particle_position_x"), ("all", "particle_position_y")],
            ("all", "particle_velocity_x"),
            logs={
                ("all", "particle_position_x"): True,
                ("all", "particle_position_y"): False,
                ("all", "particle_position_z"): False,
            },
            weight_field=None,
            deposition="cic",
        )

    # can't use CIC deposition with accumulation or fractional
    with assert_raises(RuntimeError):
        yt.create_profile(
            ad,
            [("all", "particle_position_x"), ("all", "particle_position_y")],
            ("all", "particle_velocity_x"),
            logs={
                ("all", "particle_position_x"): False,
                ("all", "particle_position_y"): False,
                ("all", "particle_position_z"): False,
            },
            weight_field=None,
            deposition="cic",
            accumulation=True,
            fractional=True,
        )


def test_profile_zero_weight():
    def DMparticles(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 1
        return filter

    def DM_in_cell_mass(field, data):
        return data["deposit", "DM_density"] * data["index", "cell_volume"]

    add_particle_filter(
        "DM", function=DMparticles, filtered_type="io", requires=["particle_type"]
    )

    _fields = (
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_mass",
        "particle_velocity_x",
        "particle_velocity_y",
        "particle_velocity_z",
        "particle_type",
    )
    _units = ("cm", "cm", "cm", "g", "cm/s", "cm/s", "cm/s", "dimensionless")
    ds = fake_random_ds(
        32, particle_fields=_fields, particle_field_units=_units, particles=16
    )

    ds.add_particle_filter("DM")
    ds.add_field(
        ("gas", "DM_cell_mass"),
        units="g",
        function=DM_in_cell_mass,
        sampling_type="cell",
    )

    sp = ds.sphere(ds.domain_center, (10, "kpc"))

    profile = yt.create_profile(
        sp,
        [("gas", "density")],
        [("gas", "radial_velocity")],
        weight_field=("gas", "DM_cell_mass"),
    )

    assert not np.any(np.isnan(profile["gas", "radial_velocity"]))


def test_profile_sph_data():
    ds = fake_sph_orientation_ds()
    # test we create a profile without raising YTIllDefinedProfile
    yt.create_profile(
        ds.all_data(),
        [("gas", "density"), ("gas", "temperature")],
        [("gas", "kinetic_energy_density")],
        weight_field=None,
    )


def test_profile_override_limits():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, units=_units)

    sp = ds.sphere(ds.domain_center, (10, "kpc"))
    obins = np.linspace(-5, 5, 10)
    profile = yt.create_profile(
        sp,
        [("gas", "density")],
        [("gas", "temperature")],
        override_bins={("gas", "density"): (obins, "g/cm**3")},
    )
    assert_equal(ds.arr(obins, "g/cm**3"), profile.x_bins)

    profile = yt.create_profile(
        sp,
        [("gas", "density"), ("stream", "dinosaurs")],
        [("gas", "temperature")],
        override_bins={
            ("gas", "density"): (obins, "g/cm**3"),
            ("stream", "dinosaurs"): obins,
        },
    )
    assert_equal(ds.arr(obins, "g/cm**3"), profile.x_bins)
    assert_equal(ds.arr(obins, "dyne"), profile.y_bins)

    profile = yt.create_profile(
        sp,
        [("gas", "density"), (("stream", "dinosaurs")), ("stream", "tribbles")],
        [("gas", "temperature")],
        override_bins={
            ("gas", "density"): (obins, "g/cm**3"),
            (("stream", "dinosaurs")): obins,
            ("stream", "tribbles"): (obins, "erg"),
        },
    )
    assert_equal(ds.arr(obins, "g/cm**3"), profile.x_bins)
    assert_equal(ds.arr(obins, "dyne"), profile.y_bins)
    assert_equal(ds.arr(obins, "erg"), profile.z_bins)


class TestBadProfiles(unittest.TestCase):

    tmpdir = None
    curdir = None

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        # clean up
        shutil.rmtree(self.tmpdir)

    @requires_module("h5py")
    def test_unequal_data_shape_profile(self):
        density = np.random.random(128)
        temperature = np.random.random(128)
        mass = np.random.random((128, 128))

        my_data = {
            ("gas", "density"): density,
            ("gas", "temperature"): temperature,
            ("gas", "mass"): mass,
        }
        fake_ds_med = {"current_time": yt.YTQuantity(10, "Myr")}
        field_types = {field: "gas" for field in my_data.keys()}
        yt.save_as_dataset(fake_ds_med, "mydata.h5", my_data, field_types=field_types)

        ds = yt.load("mydata.h5")

        with assert_raises(YTProfileDataShape):
            yt.PhasePlot(
                ds.data,
                ("gas", "temperature"),
                ("gas", "density"),
                ("gas", "mass"),
            )

    @requires_module("h5py")
    def test_unequal_bin_field_profile(self):
        density = np.random.random(128)
        temperature = np.random.random(127)
        mass = np.random.random((128, 128))

        my_data = {
            ("gas", "density"): density,
            ("gas", "temperature"): temperature,
            ("gas", "mass"): mass,
        }
        fake_ds_med = {"current_time": yt.YTQuantity(10, "Myr")}
        field_types = {field: "gas" for field in my_data.keys()}
        yt.save_as_dataset(fake_ds_med, "mydata.h5", my_data, field_types=field_types)

        ds = yt.load("mydata.h5")

        with assert_raises(YTProfileDataShape):
            yt.PhasePlot(
                ds.data,
                ("gas", "temperature"),
                ("gas", "density"),
                ("gas", "mass"),
            )

    def test_set_linear_scaling_for_none_extrema(self):
        # See Issue #3431
        # Ensures that extrema are calculated in the same way on subsequent passes
        # through the PhasePlot machinery.
        ds = fake_sph_orientation_ds()
        p = yt.PhasePlot(
            ds,
            ("all", "particle_position_spherical_theta"),
            ("all", "particle_position_spherical_radius"),
            ("all", "particle_mass"),
            weight_field=None,
        )
        p.set_log(("all", "particle_position_spherical_theta"), False)
        p.save()


def test_index_field_units():
    # see #1849
    ds = fake_random_ds(16, length_unit=2)
    ad = ds.all_data()
    icv_units = ad["index", "cell_volume"].units
    assert str(icv_units) == "code_length**3"
    gcv_units = ad["gas", "cell_volume"].units
    assert str(gcv_units) == "cm**3"
    prof = ad.profile(
        [("gas", "density"), ("gas", "velocity_x")],
        [("gas", "cell_volume"), ("index", "cell_volume")],
        weight_field=None,
    )
    assert str(prof["index", "cell_volume"].units) == "code_length**3"
    assert str(prof["gas", "cell_volume"].units) == "cm**3"


@requires_module("astropy")
def test_export_astropy():
    from yt.units.yt_array import YTArray

    ds = fake_random_ds(64)
    ad = ds.all_data()
    prof = ad.profile(
        ("index", "radius"),
        [("gas", "density"), ("gas", "velocity_x")],
        weight_field=("index", "ones"),
        n_bins=32,
    )
    # export to AstroPy table
    at1 = prof.to_astropy_table()
    assert "radius" in at1.colnames
    assert "density" in at1.colnames
    assert "velocity_x" in at1.colnames
    assert_equal(prof.x.d, at1["radius"].value)
    assert_equal(prof[("gas", "density")].d, at1["density"].value)
    assert_equal(prof[("gas", "velocity_x")].d, at1["velocity_x"].value)
    assert prof.x.units == YTArray.from_astropy(at1["radius"]).units
    assert prof[("gas", "density")].units == YTArray.from_astropy(at1["density"]).units
    assert (
        prof[("gas", "velocity_x")].units
        == YTArray.from_astropy(at1["velocity_x"]).units
    )
    assert np.all(at1.mask["density"] == prof.used)
    at2 = prof.to_astropy_table(fields=("gas", "density"), only_used=True)
    assert "radius" in at2.colnames
    assert "velocity_x" not in at2.colnames
    assert_equal(prof.x.d[prof.used], at2["radius"].value)
    assert_equal(prof[("gas", "density")].d[prof.used], at2["density"].value)


@requires_module("pandas")
def test_export_pandas():
    ds = fake_random_ds(64)
    ad = ds.all_data()
    prof = ad.profile(
        "radius",
        [("gas", "density"), ("gas", "velocity_x")],
        weight_field=("index", "ones"),
        n_bins=32,
    )
    # export to pandas DataFrame
    df1 = prof.to_dataframe()
    assert "radius" in df1.columns
    assert "density" in df1.columns
    assert "velocity_x" in df1.columns
    assert_equal(prof.x.d, df1["radius"])
    assert_equal(prof[("gas", "density")].d, np.nan_to_num(df1["density"]))
    assert_equal(prof["velocity_x"].d, np.nan_to_num(df1["velocity_x"]))
    df2 = prof.to_dataframe(fields=("gas", "density"), only_used=True)
    assert "radius" in df2.columns
    assert "velocity_x" not in df2.columns
    assert_equal(prof.x.d[prof.used], df2["radius"])
    assert_equal(prof[("gas", "density")].d[prof.used], df2["density"])
