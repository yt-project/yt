import numpy as np
import pytest
from numpy.testing import assert_equal

import yt
import yt.utilities.initial_conditions as ic
from yt.loaders import load_amr_grids, load_particles, load_uniform_grid
from yt.testing import fake_particle_ds, fake_sph_orientation_ds

# Field information


def test_stream_particles():
    num_particles = 100000
    domain_dims = (64, 64, 64)
    dens = np.random.random(domain_dims)
    x = np.random.uniform(size=num_particles)
    y = np.random.uniform(size=num_particles)
    z = np.random.uniform(size=num_particles)
    m = np.ones(num_particles)

    # Field operators and cell flagging methods

    fo = []
    fo.append(ic.TopHatSphere(0.1, [0.2, 0.3, 0.4], {"density": 2.0}))
    fo.append(ic.TopHatSphere(0.05, [0.7, 0.4, 0.75], {"density": 20.0}))

    # Add particles

    fields1 = {
        "density": dens,
        "particle_position_x": x,
        "particle_position_y": y,
        "particle_position_z": z,
        "particle_mass": m,
    }

    fields2 = fields1.copy()

    ug1 = load_uniform_grid(fields1, domain_dims, 1.0)
    ug2 = load_uniform_grid(fields2, domain_dims, 1.0, nprocs=8)

    # Check to make sure the number of particles is the same

    number_of_particles1 = np.sum([grid.NumberOfParticles for grid in ug1.index.grids])
    number_of_particles2 = np.sum([grid.NumberOfParticles for grid in ug2.index.grids])

    assert_equal(number_of_particles1, num_particles)
    assert_equal(number_of_particles1, number_of_particles2)

    for grid in ug2.index.grids:
        tot_parts = grid["io", "particle_position_x"].size
        tot_all_parts = grid["all", "particle_position_x"].size
        assert tot_parts == grid.NumberOfParticles
        assert tot_all_parts == grid.NumberOfParticles

    # Check to make sure the fields have been defined correctly

    for ptype in ("all", "io"):
        assert (
            ug1._get_field_info((ptype, "particle_position_x")).sampling_type
            == "particle"
        )
        assert (
            ug1._get_field_info((ptype, "particle_position_y")).sampling_type
            == "particle"
        )
        assert (
            ug1._get_field_info((ptype, "particle_position_z")).sampling_type
            == "particle"
        )
        assert ug1._get_field_info((ptype, "particle_mass")).sampling_type == "particle"
    assert not ug1._get_field_info(("gas", "density")).sampling_type == "particle"

    for ptype in ("all", "io"):
        assert (
            ug2._get_field_info((ptype, "particle_position_x")).sampling_type
            == "particle"
        )
        assert (
            ug2._get_field_info((ptype, "particle_position_y")).sampling_type
            == "particle"
        )
        assert (
            ug2._get_field_info((ptype, "particle_position_z")).sampling_type
            == "particle"
        )
        assert ug2._get_field_info((ptype, "particle_mass")).sampling_type == "particle"
    assert not ug2._get_field_info(("gas", "density")).sampling_type == "particle"

    # Now perform similar checks, but with multiple particle types

    num_dm_particles = 30000
    xd = np.random.uniform(size=num_dm_particles)
    yd = np.random.uniform(size=num_dm_particles)
    zd = np.random.uniform(size=num_dm_particles)
    md = np.ones(num_dm_particles)

    num_star_particles = 20000
    xs = np.random.uniform(size=num_star_particles)
    ys = np.random.uniform(size=num_star_particles)
    zs = np.random.uniform(size=num_star_particles)
    ms = 2.0 * np.ones(num_star_particles)

    dens = np.random.random(domain_dims)

    fields3 = {
        "density": dens,
        ("dm", "particle_position_x"): xd,
        ("dm", "particle_position_y"): yd,
        ("dm", "particle_position_z"): zd,
        ("dm", "particle_mass"): md,
        ("star", "particle_position_x"): xs,
        ("star", "particle_position_y"): ys,
        ("star", "particle_position_z"): zs,
        ("star", "particle_mass"): ms,
    }

    fields4 = fields3.copy()

    ug3 = load_uniform_grid(fields3, domain_dims, 1.0)
    ug4 = load_uniform_grid(fields4, domain_dims, 1.0, nprocs=8)

    # Check to make sure the number of particles is the same

    number_of_particles3 = np.sum([grid.NumberOfParticles for grid in ug3.index.grids])
    number_of_particles4 = np.sum([grid.NumberOfParticles for grid in ug4.index.grids])

    assert_equal(number_of_particles3, num_dm_particles + num_star_particles)
    assert_equal(number_of_particles3, number_of_particles4)

    for grid in ug4.index.grids:
        tot_parts = grid["dm", "particle_position_x"].size
        tot_parts += grid["star", "particle_position_x"].size
        tot_all_parts = grid["all", "particle_position_x"].size
        assert tot_parts == grid.NumberOfParticles
        assert tot_all_parts == grid.NumberOfParticles

    # Check to make sure the fields have been defined correctly

    for ptype in ("dm", "star"):
        assert (
            ug3._get_field_info((ptype, "particle_position_x")).sampling_type
            == "particle"
        )
        assert (
            ug3._get_field_info((ptype, "particle_position_y")).sampling_type
            == "particle"
        )
        assert (
            ug3._get_field_info((ptype, "particle_position_z")).sampling_type
            == "particle"
        )
        assert ug3._get_field_info((ptype, "particle_mass")).sampling_type == "particle"
        assert (
            ug4._get_field_info((ptype, "particle_position_x")).sampling_type
            == "particle"
        )
        assert (
            ug4._get_field_info((ptype, "particle_position_y")).sampling_type
            == "particle"
        )
        assert (
            ug4._get_field_info((ptype, "particle_position_z")).sampling_type
            == "particle"
        )
        assert ug4._get_field_info((ptype, "particle_mass")).sampling_type == "particle"


def test_load_particles_types():
    num_particles = 10000

    data1 = {
        "particle_position_x": np.random.random(size=num_particles),
        "particle_position_y": np.random.random(size=num_particles),
        "particle_position_z": np.random.random(size=num_particles),
        "particle_mass": np.ones(num_particles),
    }

    ds1 = load_particles(data1)
    ds1.index

    assert set(ds1.particle_types) == {"all", "io", "nbody"}

    dd = ds1.all_data()

    for ax in "xyz":
        assert dd["io", f"particle_position_{ax}"].size == num_particles
        assert dd["all", f"particle_position_{ax}"].size == num_particles
        assert dd["nbody", f"particle_position_{ax}"].size == num_particles

    num_dm_particles = 10000
    num_star_particles = 50000
    num_tot_particles = num_dm_particles + num_star_particles

    data2 = {
        ("dm", "particle_position_x"): np.random.random(size=num_dm_particles),
        ("dm", "particle_position_y"): np.random.random(size=num_dm_particles),
        ("dm", "particle_position_z"): np.random.random(size=num_dm_particles),
        ("dm", "particle_mass"): np.ones(num_dm_particles),
        ("star", "particle_position_x"): np.random.random(size=num_star_particles),
        ("star", "particle_position_y"): np.random.random(size=num_star_particles),
        ("star", "particle_position_z"): np.random.random(size=num_star_particles),
        ("star", "particle_mass"): 2.0 * np.ones(num_star_particles),
    }

    ds2 = load_particles(data2)
    ds2.index

    # We use set here because we don't care about the order and we just need
    # the elements to be correct
    assert set(ds2.particle_types) == {"all", "star", "dm", "nbody"}

    dd = ds2.all_data()

    for ax in "xyz":
        npart = 0
        for ptype in ds2.particle_types_raw:
            npart += dd[ptype, f"particle_position_{ax}"].size
        assert npart == num_tot_particles
        assert dd["all", f"particle_position_{ax}"].size == num_tot_particles


def test_load_particles_sph_types():
    num_particles = 10000

    data = {
        ("gas", "particle_position_x"): np.random.random(size=num_particles),
        ("gas", "particle_position_y"): np.random.random(size=num_particles),
        ("gas", "particle_position_z"): np.random.random(size=num_particles),
        ("gas", "particle_velocity_x"): np.random.random(size=num_particles),
        ("gas", "particle_velocity_y"): np.random.random(size=num_particles),
        ("gas", "particle_velocity_z"): np.random.random(size=num_particles),
        ("gas", "particle_mass"): np.ones(num_particles),
        ("gas", "density"): np.ones(num_particles),
        ("gas", "smoothing_length"): np.ones(num_particles),
        ("dm", "particle_position_x"): np.random.random(size=num_particles),
        ("dm", "particle_position_y"): np.random.random(size=num_particles),
        ("dm", "particle_position_z"): np.random.random(size=num_particles),
        ("dm", "particle_velocity_x"): np.random.random(size=num_particles),
        ("dm", "particle_velocity_y"): np.random.random(size=num_particles),
        ("dm", "particle_velocity_z"): np.random.random(size=num_particles),
        ("dm", "particle_mass"): np.ones(num_particles),
    }

    ds = load_particles(data)

    assert set(ds.particle_types) == {"gas", "dm"}
    assert ds._sph_ptypes == ("gas",)

    data.update(
        {
            ("cr_gas", "particle_position_x"): np.random.random(size=num_particles),
            ("cr_gas", "particle_position_y"): np.random.random(size=num_particles),
            ("cr_gas", "particle_position_z"): np.random.random(size=num_particles),
            ("cr_gas", "particle_velocity_x"): np.random.random(size=num_particles),
            ("cr_gas", "particle_velocity_y"): np.random.random(size=num_particles),
            ("cr_gas", "particle_velocity_z"): np.random.random(size=num_particles),
            ("cr_gas", "particle_mass"): np.ones(num_particles),
            ("cr_gas", "density"): np.ones(num_particles),
            ("cr_gas", "smoothing_length"): np.ones(num_particles),
        }
    )

    with pytest.raises(
        ValueError, match="Multiple SPH particle types are currently not supported!"
    ):
        load_particles(data)


def test_load_particles_with_data_source():
    ds1 = fake_particle_ds()

    # Load from dataset
    ad = ds1.all_data()
    fields = [("all", "particle_mass")]
    fields += [("all", f"particle_position_{ax}") for ax in "xyz"]
    data = {field: ad[field] for field in fields}
    ds2 = load_particles(data, data_source=ad)

    def in_cgs(quan):
        return quan.in_cgs().v

    # Test bbox is parsed correctly
    for attr in ["domain_left_edge", "domain_right_edge"]:
        assert np.allclose(in_cgs(getattr(ds1, attr)), in_cgs(getattr(ds2, attr)))

    # Test sim_time is parsed correctly
    assert in_cgs(ds1.current_time) == in_cgs(ds2.current_time)

    # Test code units are parsed correctly
    def get_cu(ds, dim):
        return ds.quan(1, "code_" + dim)

    for dim in ["length", "mass", "time", "velocity", "magnetic"]:
        assert in_cgs(get_cu(ds1, dim)) == in_cgs(get_cu(ds2, dim))


def test_add_sph_fields():
    ds = fake_particle_ds()
    ds.index
    assert set(ds.particle_types) == {"io", "all", "nbody"}

    ds.add_sph_fields()
    assert set(ds.particle_types) == {"io", "all"}
    assert ("io", "smoothing_length") in ds.field_list
    assert ("io", "density") in ds.field_list


def test_particles_outside_domain():
    np.random.seed(0x4D3D3D3)
    posx_arr = np.random.uniform(low=-1.6, high=1.5, size=1000)
    posy_arr = np.random.uniform(low=-1.5, high=1.5, size=1000)
    posz_arr = np.random.uniform(low=-1.5, high=1.5, size=1000)
    dens_arr = np.random.random((16, 16, 16))
    data = {
        "density": dens_arr,
        "particle_position_x": posx_arr,
        "particle_position_y": posy_arr,
        "particle_position_z": posz_arr,
    }
    bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
    ds = load_uniform_grid(data, (16, 16, 16), bbox=bbox, nprocs=4)
    wh = (posx_arr < bbox[0, 0]).nonzero()[0]
    assert wh.size == 1000 - ds.particle_type_counts["io"]
    ad = ds.all_data()
    assert ds.particle_type_counts["io"] == ad["all", "particle_position_x"].size


def test_stream_sph_projection():
    ds = fake_sph_orientation_ds()
    proj = ds.proj(("gas", "density"), 2)
    frb = proj.to_frb(ds.domain_width[0], (256, 256))
    image = frb["gas", "density"]
    assert image.max() > 0
    assert image.shape == (256, 256)


@pytest.mark.parametrize("loader", (load_uniform_grid, load_amr_grids))
def test_stream_non_cartesian_particles(loader):
    eps = 1e-6
    r, theta, phi = np.mgrid[
        0.0 : 1.0 - eps : 64j, 0.0 : np.pi - eps : 64j, 0.0 : 2.0 * np.pi - eps : 64j
    ]
    np.random.seed(0x4D3D3D3)
    ind = np.random.randint(0, 64 * 64 * 64, size=1000)

    particle_position_r = r.ravel()[ind]
    particle_position_theta = theta.ravel()[ind]
    particle_position_phi = phi.ravel()[ind]

    ds = load_uniform_grid(
        {
            "density": r,
            "temperature": phi,
            "entropy": phi,
            "particle_position_r": particle_position_r,
            "particle_position_theta": particle_position_theta,
            "particle_position_phi": particle_position_phi,
        },
        (64, 64, 64),
        bbox=np.array([[0.0, 1.0], [0.0, np.pi], [0.0, 2.0 * np.pi]]),
        geometry="spherical",
    )

    dd = ds.all_data()
    assert_equal(dd["all", "particle_position_r"].v, particle_position_r)
    assert_equal(dd["all", "particle_position_phi"].v, particle_position_phi)
    assert_equal(dd["all", "particle_position_theta"].v, particle_position_theta)


def test_stream_non_cartesian_particles_amr():
    eps = 1e-6
    r, theta, phi = np.mgrid[
        0.0 : 1.0 - eps : 64j, 0.0 : np.pi - eps : 64j, 0.0 : 2.0 * np.pi - eps : 64j
    ]
    np.random.seed(0x4D3D3D3)
    ind = np.random.randint(0, 64 * 64 * 64, size=1000)

    particle_position_r = r.ravel()[ind]
    particle_position_theta = theta.ravel()[ind]
    particle_position_phi = phi.ravel()[ind]

    ds = load_amr_grids(
        [
            {
                "density": r,
                "temperature": phi,
                "entropy": phi,
                "particle_position_r": particle_position_r,
                "particle_position_theta": particle_position_theta,
                "particle_position_phi": particle_position_phi,
                "dimensions": [64, 64, 64],
                "level": 0,
                "left_edge": [0.0, 0.0, 0.0],
                "right_edge": [1.0, np.pi, 2.0 * np.pi],
            }
        ],
        (64, 64, 64),
        bbox=np.array([[0.0, 1.0], [0.0, np.pi], [0.0, 2.0 * np.pi]]),
        geometry="spherical",
    )

    dd = ds.all_data()
    assert_equal(dd["all", "particle_position_r"].v, particle_position_r)
    assert_equal(dd["all", "particle_position_phi"].v, particle_position_phi)
    assert_equal(dd["all", "particle_position_theta"].v, particle_position_theta)


@pytest.fixture
def sph_dataset_with_integer_index():
    num_particles = 100

    data = {
        ("gas", "particle_position_x"): np.linspace(0.0, 1.0, num_particles),
        ("gas", "particle_position_y"): np.linspace(0.0, 1.0, num_particles),
        ("gas", "particle_position_z"): np.linspace(0.0, 1.0, num_particles),
        ("gas", "particle_mass"): np.ones(num_particles),
        ("gas", "density"): np.ones(num_particles),
        ("gas", "smoothing_length"): np.ones(num_particles) * 0.1,
        ("gas", "particle_index"): np.arange(0, num_particles),
    }

    ds = load_particles(data)
    return ds


def test_particle_dtypes_selection(sph_dataset_with_integer_index):
    # these operations will preserve data type
    ds = sph_dataset_with_integer_index
    ad = ds.all_data()
    assert ad["gas", "particle_index"].dtype == np.int64

    min_max = ad.quantities.extrema(("gas", "particle_index"))
    assert min_max.dtype == np.int64

    # check that subselections preserve type
    le = ds.domain_center - ds.domain_width / 10.0
    re = ds.domain_center + ds.domain_width / 10.0
    reg = ds.region(ds.domain_center, le, re)
    assert reg["gas", "particle_index"].dtype == np.int64

    vals = ds.slice(0, ds.domain_center[0])["gas", "particle_index"]
    assert vals.max() > 0
    assert vals.dtype == np.int64


def test_particle_dtypes_operations(sph_dataset_with_integer_index):
    # these operations will not preserve dtype (will be cast to float64).
    # note that the numerical outputs of these operations are not
    # physical (projecting the particle index does not make any physical
    # sense), but they do make sure the methods run in case any frontends
    # start setting physical fields with different data types.

    ds = sph_dataset_with_integer_index

    field = ("gas", "particle_index")
    frb = ds.proj(field, 2).to_frb(ds.domain_width[0], (64, 64))
    image = frb["gas", "particle_index"]
    assert image.max() > 0

    off_axis_prj = yt.off_axis_projection(
        ds,
        ds.domain_center,
        [0.5, 0.5, 0.5],
        ds.domain_width,
        (64, 64),
        ("gas", "particle_index"),
        weight=None,
    )
    assert off_axis_prj.max() > 0

    source = ds.all_data()
    custom_bins = np.linspace(-0.5, 99.5, 101)
    profile = source.profile(
        [("gas", "particle_index")],
        [
            ("gas", "mass"),
            ("gas", "particle_index"),
        ],
        weight_field=None,
        override_bins={("gas", "particle_index"): custom_bins},
        logs={("gas", "particle_index"): False},
    )
    assert profile.x.min() == 0.0
    assert profile.x.max() == 99.0
    assert np.all(profile["gas", "mass"] == 1.0)
    assert np.all(profile["gas", "particle_index"] == profile.x)

    # check a particle filter on the index
    p_index_min = 50

    def index_filter(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_index")] >= p_index_min
        return filter

    yt.add_particle_filter(
        "index_filter",
        function=index_filter,
        filtered_type="all",
        requires=["particle_index"],
    )

    ds.add_particle_filter("index_filter")
    assert ds.all_data()["index_filter", "particle_index"].min() == p_index_min
