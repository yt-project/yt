import numpy as np

from yt.testing import \
    assert_equal, \
    fake_amr_ds, \
    fake_particle_ds, \
    fake_random_ds

# This will test the "dataset access" method.

def test_box_creation():
    ds = fake_random_ds(32, length_unit=2)
    left_edge = ds.arr([0.2, 0.2, 0.2], 'cm')
    right_edge = ds.arr([0.6, 0.6, 0.6], 'cm')
    center = (left_edge + right_edge)/2

    boxes = [
        ds.box(left_edge, right_edge),
        ds.box(0.5*np.array(left_edge), 0.5*np.array(right_edge)),
        ds.box((0.5*left_edge).tolist(), (0.5*right_edge).tolist())
    ]

    region = ds.region(center, left_edge, right_edge)

    for b in boxes:
        assert_equal(b.left_edge, region.left_edge)
        assert_equal(b.right_edge, region.right_edge)
        assert_equal(b.center, region.center)

def test_region_from_d():
    ds = fake_amr_ds(fields=["density"])
    # We'll do a couple here

    # First, no string units
    reg1 = ds.r[0.2:0.3,0.4:0.6,:]
    reg2 = ds.region([0.25, 0.5, 0.5], [0.2, 0.4, 0.0], [0.3, 0.6, 1.0])
    yield assert_equal, reg1["density"], reg2["density"]

    # Now, string units in some -- 1.0 == cm
    reg1 = ds.r[(0.1, 'cm'):(0.5, 'cm'), :, (0.25, 'cm'): (0.35, 'cm')]
    reg2 = ds.region([0.3, 0.5, 0.3], [0.1, 0.0, 0.25], [0.5, 1.0, 0.35])
    yield assert_equal, reg1["density"], reg2["density"]

    # Now, string units in some -- 1.0 == cm
    reg1 = ds.r[(0.1, 'cm'):(0.5, 'cm'), :, 0.25:0.35]
    reg2 = ds.region([0.3, 0.5, 0.3], [0.1, 0.0, 0.25], [0.5, 1.0, 0.35])
    yield assert_equal, reg1["density"], reg2["density"]

    # And, lots of : usage!
    reg1 = ds.r[:, :, :]
    reg2 = ds.all_data()
    yield assert_equal, reg1["density"], reg2["density"]

def test_accessing_all_data():
    # This will test first that we can access all_data, and next that we can
    # access it multiple times and get the *same object*.
    ds = fake_amr_ds(fields=["density"])
    dd = ds.all_data()
    yield assert_equal, ds.r["density"], dd["density"]
    # Now let's assert that it's the same object
    rho = ds.r["density"]
    rho *= 2.0
    yield assert_equal, dd["density"]*2.0, ds.r["density"]
    yield assert_equal, dd["gas", "density"]*2.0, ds.r["gas", "density"]

def test_particle_counts():
    ds = fake_random_ds(16, particles=100)
    assert ds.particle_type_counts == {'io': 100}

    pds = fake_particle_ds(npart=128)
    assert pds.particle_type_counts == {'io': 128}
