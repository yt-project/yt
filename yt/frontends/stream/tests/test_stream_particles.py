import numpy as np
import pytest

import yt.utilities.initial_conditions as ic
from yt.loaders import load_amr_grids, load_particles, load_uniform_grid
from yt.testing import assert_equal, fake_particle_ds, fake_sph_orientation_ds

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
            ug1._get_field_info(ptype, "particle_position_x").sampling_type
            == "particle"
        )
        assert (
            ug1._get_field_info(ptype, "particle_position_y").sampling_type
            == "particle"
        )
        assert (
            ug1._get_field_info(ptype, "particle_position_z").sampling_type
            == "particle"
        )
        assert ug1._get_field_info(ptype, "particle_mass").sampling_type == "particle"
    assert not ug1._get_field_info("gas", "density").sampling_type == "particle"

    for ptype in ("all", "io"):
        assert (
            ug2._get_field_info(ptype, "particle_position_x").sampling_type
            == "particle"
        )
        assert (
            ug2._get_field_info(ptype, "particle_position_y").sampling_type
            == "particle"
        )
        assert (
            ug2._get_field_info(ptype, "particle_position_z").sampling_type
            == "particle"
        )
        assert ug2._get_field_info(ptype, "particle_mass").sampling_type == "particle"
    assert not ug2._get_field_info("gas", "density").sampling_type == "particle"

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
            ug3._get_field_info(ptype, "particle_position_x").sampling_type
            == "particle"
        )
        assert (
            ug3._get_field_info(ptype, "particle_position_y").sampling_type
            == "particle"
        )
        assert (
            ug3._get_field_info(ptype, "particle_position_z").sampling_type
            == "particle"
        )
        assert ug3._get_field_info(ptype, "particle_mass").sampling_type == "particle"
        assert (
            ug4._get_field_info(ptype, "particle_position_x").sampling_type
            == "particle"
        )
        assert (
            ug4._get_field_info(ptype, "particle_position_y").sampling_type
            == "particle"
        )
        assert (
            ug4._get_field_info(ptype, "particle_position_z").sampling_type
            == "particle"
        )
        assert ug4._get_field_info(ptype, "particle_mass").sampling_type == "particle"


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

    assert set(ds2.particle_types) == {"all", "star", "dm", "nbody"}

    dd = ds2.all_data()

    for ax in "xyz":
        npart = 0
        for ptype in ds2.particle_types_raw:
            npart += dd[ptype, f"particle_position_{ax}"].size
        assert npart == num_tot_particles
        assert dd["all", f"particle_position_{ax}"].size == num_tot_particles


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
    data = dict(
        density=dens_arr,
        particle_position_x=posx_arr,
        particle_position_y=posy_arr,
        particle_position_z=posz_arr,
    )
    bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
    ds = load_uniform_grid(data, (16, 16, 16), bbox=bbox, nprocs=4)
    wh = (posx_arr < bbox[0, 0]).nonzero()[0]
    assert wh.size == 1000 - ds.particle_type_counts["io"]
    ad = ds.all_data()
    assert ds.particle_type_counts["io"] == ad[("all", "particle_position_x")].size


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
