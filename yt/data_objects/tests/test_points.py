import numpy as np
import yt

from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    requires_file

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_point_creation():
    ds = fake_random_ds(16)
    p1 = ds.point(ds.domain_center)
    p2 = ds.point([0.5, 0.5, 0.5])
    p3 = ds.point([0.5, 0.5, 0.5]*yt.units.cm)

    # ensure all three points are really at the same position
    for fname in 'xyz':
        assert_equal(p1[fname], p2[fname])
        assert_equal(p1[fname], p3[fname])

def test_domain_point():
    nparticles = 3
    ds = fake_random_ds(16, particles=nparticles)
    p = ds.point(ds.domain_center)

    # ensure accessing one field works, store for comparison later
    point_den = p['density']
    point_vel = p['velocity_x']

    ad = ds.all_data()
    ppos = ad['all', 'particle_position']

    fpoint_den = ds.find_field_values_at_point('density', ds.domain_center)

    fpoint_den_vel = ds.find_field_values_at_point(
        ['density', 'velocity_x'], ds.domain_center)

    assert_equal(point_den, fpoint_den)
    assert_equal(point_den, fpoint_den_vel[0])
    assert_equal(point_vel, fpoint_den_vel[1])

    ppos_den = ds.find_field_values_at_points('density', ppos)
    ppos_vel = ds.find_field_values_at_points('velocity_x', ppos)
    ppos_den_vel = ds.find_field_values_at_points(
        ['density', 'velocity_x'], ppos)

    assert_equal(ppos_den.shape, (nparticles,))
    assert_equal(ppos_vel.shape, (nparticles,))
    assert_equal(len(ppos_den_vel), 2)
    assert_equal(ppos_den_vel[0], ppos_den)
    assert_equal(ppos_den_vel[1], ppos_vel)

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"

@requires_file(g30)
def test_fast_find_field_values_at_points():
    ds = yt.load(g30)
    ad = ds.all_data()
    # right now this is slow for large numbers of particles, so randomly
    # sample 100 particles
    nparticles = 100
    ppos = ad['all', 'particle_position']
    ppos = ppos[np.random.random_integers(0, len(ppos), size=nparticles)]

    ppos_den = ds.find_field_values_at_points('density', ppos)
    ppos_vel = ds.find_field_values_at_points('velocity_x', ppos)
    ppos_den_vel = ds.find_field_values_at_points(
        ['density', 'velocity_x'], ppos)

    assert_equal(ppos_den.shape, (nparticles,))
    assert_equal(ppos_vel.shape, (nparticles,))
    assert_equal(len(ppos_den_vel), 2)
    assert_equal(ppos_den_vel[0], ppos_den)
    assert_equal(ppos_den_vel[1], ppos_vel)
