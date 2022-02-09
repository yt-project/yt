import numpy as np

import yt
from yt.testing import assert_equal, fake_random_ds


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_point_creation():
    ds = fake_random_ds(16)
    p1 = ds.point(ds.domain_center)
    p2 = ds.point([0.5, 0.5, 0.5])
    p3 = ds.point([0.5, 0.5, 0.5] * yt.units.cm)

    # ensure all three points are really at the same position
    for fname in "xyz":
        assert_equal(p1["index", fname], p2["index", fname])
        assert_equal(p1["index", fname], p3["index", fname])


def test_domain_point():
    nparticles = 3
    ds = fake_random_ds(16, particles=nparticles)
    p = ds.point(ds.domain_center)

    # ensure accessing one field works, store for comparison later
    point_den = p[("gas", "density")]
    point_vel = p[("gas", "velocity_x")]

    ad = ds.all_data()
    ppos = ad["all", "particle_position"]

    fpoint_den = ds.find_field_values_at_point(("gas", "density"), ds.domain_center)

    fpoint_den_vel = ds.find_field_values_at_point(
        [("gas", "density"), ("gas", "velocity_x")], ds.domain_center
    )

    assert_equal(point_den, fpoint_den)
    assert_equal(point_den, fpoint_den_vel[0])
    assert_equal(point_vel, fpoint_den_vel[1])

    ppos_den = ds.find_field_values_at_points(("gas", "density"), ppos)
    ppos_vel = ds.find_field_values_at_points(("gas", "velocity_x"), ppos)
    ppos_den_vel = ds.find_field_values_at_points(
        [("gas", "density"), ("gas", "velocity_x")], ppos
    )

    assert_equal(ppos_den.shape, (nparticles,))
    assert_equal(ppos_vel.shape, (nparticles,))
    assert_equal(len(ppos_den_vel), 2)
    assert_equal(ppos_den_vel[0], ppos_den)
    assert_equal(ppos_den_vel[1], ppos_vel)


def test_fast_find_field_values_at_points():
    ds = fake_random_ds(64, nprocs=8, particles=16**3)
    ad = ds.all_data()
    # right now this is slow for large numbers of particles, so randomly
    # sample 100 particles
    nparticles = 100
    ppos = ad["all", "particle_position"]
    ppos = ppos[np.random.randint(low=0, high=len(ppos), size=nparticles)]

    ppos_den = ds.find_field_values_at_points(("gas", "density"), ppos)
    ppos_vel = ds.find_field_values_at_points(("gas", "velocity_x"), ppos)
    ppos_den_vel = ds.find_field_values_at_points(
        [("gas", "density"), ("gas", "velocity_x")], ppos
    )

    assert_equal(ppos_den.shape, (nparticles,))
    assert_equal(ppos_vel.shape, (nparticles,))
    assert_equal(len(ppos_den_vel), 2)
    assert_equal(ppos_den_vel[0], ppos_den)
    assert_equal(ppos_den_vel[1], ppos_vel)
