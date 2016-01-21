from yt.testing import fake_random_ds, assert_equal

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

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
