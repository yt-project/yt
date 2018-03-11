import yt
import numpy as np

from yt.data_objects.particle_filters import \
    add_particle_filter
from yt.data_objects.profiles import \
    Profile1D, \
    Profile2D, \
    Profile3D, \
    create_profile
from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_raises, \
    assert_rel_equal, \
    requires_file
from yt.utilities.exceptions import \
    YTIllDefinedProfile
from yt.visualization.profile_plotter import \
    ProfilePlot, \
    PhasePlot

_fields = ("density", "temperature", "dinosaurs", "tribbles")
_units = ("g/cm**3", "K", "dyne", "erg")


def test_profiles():
    ds = fake_random_ds(64, nprocs = 8, fields = _fields, units = _units)
    nv = ds.domain_dimensions.prod()
    dd = ds.all_data()
    (rmi, rma), (tmi, tma), (dmi, dma) = dd.quantities["Extrema"](
        ["density", "temperature", "dinosaurs"])
    rt, tt, dt = dd.quantities["TotalQuantity"](
        ["density", "temperature", "dinosaurs"])

    e1, e2 = 0.9, 1.1
    for nb in [8, 16, 32, 64]:
        for input_units in ['mks', 'cgs']:
            for ex in [rmi, rma, tmi, tma, dmi, dma]:
                getattr(ex, 'convert_to_%s' % input_units)()
            # We log all the fields or don't log 'em all.  No need to do them
            # individually.
            for lf in [True, False]:
                direct_profile = Profile1D(
                    dd, "density", nb, rmi*e1, rma*e2, lf, weight_field=None)
                direct_profile.add_fields(["ones", "temperature"])

                indirect_profile_s = create_profile(
                    dd, "density", ["ones", "temperature"], n_bins=nb,
                    extrema={'density': (rmi*e1, rma*e2)}, logs={'density': lf},
                    weight_field=None)

                indirect_profile_t = create_profile(
                    dd, ("gas", "density"),
                    [("index", "ones"), ("gas", "temperature")], n_bins=nb,
                    extrema={'density': (rmi*e1, rma*e2)}, logs={'density': lf},
                    weight_field=None)

                for p1d in [direct_profile, indirect_profile_s,
                            indirect_profile_t]:
                    assert_equal(p1d["index", "ones"].sum(), nv)
                    assert_rel_equal(tt, p1d["gas", "temperature"].sum(), 7)

                p2d = Profile2D(
                    dd,
                    "density",     nb, rmi*e1, rma*e2, lf,
                    "temperature", nb, tmi*e1, tma*e2, lf,
                    weight_field=None)
                p2d.add_fields(["ones", "temperature"])
                assert_equal(p2d["ones"].sum(), nv)
                assert_rel_equal(tt, p2d["temperature"].sum(), 7)

                p3d = Profile3D(
                    dd,
                    "density",     nb, rmi*e1, rma*e2, lf,
                    "temperature", nb, tmi*e1, tma*e2, lf,
                    "dinosaurs",   nb, dmi*e1, dma*e2, lf,
                    weight_field=None)
                p3d.add_fields(["ones", "temperature"])
                assert_equal(p3d["ones"].sum(), nv)
                assert_rel_equal(tt, p3d["temperature"].sum(), 7)

        p1d = Profile1D(dd, "x", nb, 0.0, 1.0, False,
                        weight_field = None)
        p1d.add_fields("ones")
        av = nv / nb
        assert_equal(p1d["ones"], np.ones(nb)*av)

        # We re-bin ones with a weight now
        p1d = Profile1D(dd, "x", nb, 0.0, 1.0, False,
                        weight_field = "temperature")
        p1d.add_fields(["ones"])
        assert_equal(p1d["ones"], np.ones(nb))

        # Verify we can access "ones" after adding a new field
        # See issue 988
        p1d.add_fields(["density"])
        assert_equal(p1d["ones"], np.ones(nb))

        p2d = Profile2D(dd, "x", nb, 0.0, 1.0, False,
                            "y", nb, 0.0, 1.0, False,
                            weight_field = None)
        p2d.add_fields("ones")
        av = nv / nb**2
        assert_equal(p2d["ones"], np.ones((nb, nb))*av)

        # We re-bin ones with a weight now
        p2d = Profile2D(dd, "x", nb, 0.0, 1.0, False,
                            "y", nb, 0.0, 1.0, False,
                            weight_field = "temperature")
        p2d.add_fields(["ones"])
        assert_equal(p2d["ones"], np.ones((nb, nb)))

        p3d = Profile3D(dd, "x", nb, 0.0, 1.0, False,
                            "y", nb, 0.0, 1.0, False,
                            "z", nb, 0.0, 1.0, False,
                            weight_field = None)
        p3d.add_fields("ones")
        av = nv / nb**3
        assert_equal(p3d["ones"], np.ones((nb, nb, nb))*av)

        # We re-bin ones with a weight now
        p3d = Profile3D(dd, "x", nb, 0.0, 1.0, False,
                            "y", nb, 0.0, 1.0, False,
                            "z", nb, 0.0, 1.0, False,
                            weight_field = "temperature")
        p3d.add_fields(["ones"])
        assert_equal(p3d["ones"], np.ones((nb,nb,nb)))

        p2d = create_profile(dd, ('gas', 'density'), ('gas', 'temperature'),
                             weight_field=('gas', 'cell_mass'),
                             extrema={'density': (None, rma*e2)})
        assert_equal(p2d.x_bins[0], rmi - np.spacing(rmi))
        assert_equal(p2d.x_bins[-1], rma*e2)

        p2d = create_profile(dd, ('gas', 'density'), ('gas', 'temperature'),
                             weight_field=('gas', 'cell_mass'),
                             extrema={'density': (rmi*e2, None)})
        assert_equal(p2d.x_bins[0], rmi*e2)
        assert_equal(p2d.x_bins[-1], rma + np.spacing(rma))


extrema_s = {'particle_position_x': (0, 1)}
logs_s = {'particle_position_x': False}

extrema_t = {('all', 'particle_position_x'): (0, 1)}
logs_t = {('all', 'particle_position_x'): False}

def test_particle_profiles():
    for nproc in [1, 2, 4, 8]:
        ds = fake_random_ds(32, nprocs=nproc, particles = 32**3)
        dd = ds.all_data()

        p1d = Profile1D(dd, "particle_position_x", 128,
                        0.0, 1.0, False, weight_field = None)
        p1d.add_fields(["particle_ones"])
        assert_equal(p1d["particle_ones"].sum(), 32**3)

        p1d = create_profile(dd, ["particle_position_x"], ["particle_ones"],
                             weight_field=None, n_bins=128, extrema=extrema_s,
                             logs=logs_s)
        assert_equal(p1d["particle_ones"].sum(), 32**3)

        p1d = create_profile(dd,
                             [("all", "particle_position_x")],
                             [("all", "particle_ones")],
                             weight_field=None, n_bins=128, extrema=extrema_t,
                             logs=logs_t)
        assert_equal(p1d["particle_ones"].sum(), 32**3)

        p2d = Profile2D(dd, "particle_position_x", 128, 0.0, 1.0, False,
                            "particle_position_y", 128, 0.0, 1.0, False,
                        weight_field = None)
        p2d.add_fields(["particle_ones"])
        assert_equal(p2d["particle_ones"].sum(), 32**3)

        p3d = Profile3D(dd, "particle_position_x", 128, 0.0, 1.0, False,
                            "particle_position_y", 128, 0.0, 1.0, False,
                            "particle_position_z", 128, 0.0, 1.0, False,
                        weight_field = None)
        p3d.add_fields(["particle_ones"])
        assert_equal(p3d["particle_ones"].sum(), 32**3)

def test_mixed_particle_mesh_profiles():
    ds = fake_random_ds(32, particles=10)
    ad = ds.all_data()
    assert_raises(
        YTIllDefinedProfile, ProfilePlot, ad, 'radius', 'particle_mass')
    assert_raises(
        YTIllDefinedProfile, ProfilePlot, ad, 'radius',
        ['particle_mass', 'particle_ones'])
    assert_raises(
        YTIllDefinedProfile, ProfilePlot, ad, 'radius',
        ['particle_mass', 'ones'])
    assert_raises(
        YTIllDefinedProfile, ProfilePlot, ad, 'particle_radius', 'particle_mass',
        'cell_mass')
    assert_raises(
        YTIllDefinedProfile, ProfilePlot, ad, 'radius', 'cell_mass',
        'particle_ones')

    assert_raises(
        YTIllDefinedProfile, PhasePlot, ad, 'radius', 'particle_mass',
        'velocity_x')
    assert_raises(
        YTIllDefinedProfile, PhasePlot, ad, 'particle_radius', 'particle_mass',
        'cell_mass')
    assert_raises(
        YTIllDefinedProfile, PhasePlot, ad, 'radius', 'cell_mass',
        'particle_ones')
    assert_raises(
        YTIllDefinedProfile, PhasePlot, ad, 'particle_radius', 'particle_mass',
        'particle_ones')

def test_particle_profile_negative_field():
    # see Issue #1340
    n_particles = int(1e4)

    ppx, ppy, ppz = np.random.normal(size=[3, n_particles])
    pvx, pvy, pvz = - np.ones((3, n_particles))

    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_velocity_x': pvx,
            'particle_velocity_y': pvy,
            'particle_velocity_z': pvz}

    bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])
    ds = yt.load_particles(data, bbox=bbox)
    ad = ds.all_data()

    profile = yt.create_profile(
        ad,
        ["particle_position_x", "particle_position_y"],
        "particle_velocity_x",
        logs = {'particle_position_x': True,
                'particle_position_y': True,
                'particle_position_z': True},
        weight_field=None)
    assert profile['particle_velocity_x'].min() < 0
    assert profile.x_bins.min() > 0
    assert profile.y_bins.min() > 0

    profile = yt.create_profile(
        ad,
        ["particle_position_x", "particle_position_y"],
        "particle_velocity_x",
        weight_field=None)
    assert profile['particle_velocity_x'].min() < 0
    assert profile.x_bins.min() < 0
    assert profile.y_bins.min() < 0

    # can't use CIC deposition with log-scaled bin fields
    with assert_raises(RuntimeError):
        yt.create_profile(
            ad,
            ["particle_position_x", "particle_position_y"],
            "particle_velocity_x",
            logs = {'particle_position_x': True,
                    'particle_position_y': False,
                    'particle_position_z': False},
            weight_field=None, deposition='cic')

    # can't use CIC deposition with accumulation or fractional
    with assert_raises(RuntimeError):
        yt.create_profile(
            ad,
            ["particle_position_x", "particle_position_y"],
            "particle_velocity_x",
            logs = {'particle_position_x': False,
                    'particle_position_y': False,
                    'particle_position_z': False},
            weight_field=None, deposition='cic',
            accumulation=True, fractional=True)

@requires_file("IsolatedGalaxy/galaxy0030/galaxy0030")
def test_profile_zero_weight():
    def DMparticles(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 1
        return filter

    def DM_in_cell_mass(field, data):
        return data['deposit', 'DM_density']*data['index', 'cell_volume']

    add_particle_filter("DM", function=DMparticles,
                        filtered_type='io', requires=["particle_type"])

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    ds.add_particle_filter('DM')

    ds.add_field(("gas", "DM_cell_mass"), units="g", function=DM_in_cell_mass,
                 sampling_type='cell')

    sp = ds.sphere(ds.domain_center, (10, 'kpc'))

    profile = yt.create_profile(sp,
                                [("gas", "density")],
                                [("gas", "temperature")],
                                weight_field=("gas", "DM_cell_mass"))

    assert not np.any(np.isnan(profile['gas', 'temperature']))
