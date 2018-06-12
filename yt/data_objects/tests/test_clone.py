from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_array_equal

def test_clone_sphere():
    # Now we test that we can get different radial velocities based on field
    # parameters.

    # Get the first sphere
    ds = fake_random_ds(16, fields = ("density",
      "velocity_x", "velocity_y", "velocity_z"))
    sp0 = ds.sphere(ds.domain_center, 0.25)

    assert_equal(list(sp0.keys()), [])

    sp1 = sp0.clone()
    sp0["density"]
    assert_equal(list(sp0.keys()), (("gas","density"),))
    assert_equal(list(sp1.keys()), [])

    sp1["density"]

    assert_array_equal(sp0["density"], sp1["density"])

def test_clone_cut_region():
    ds = fake_random_ds(64, nprocs=4, fields=("density", "temperature"))
    dd = ds.all_data()
    reg1 = dd.cut_region(["obj['temperature'] > 0.5", "obj['density'] < 0.75"])
    reg2 = reg1.clone()
    assert_array_equal(reg1['density'], reg2['density'])
