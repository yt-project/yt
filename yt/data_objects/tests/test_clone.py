from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_array_equal

def test_clone_sphere()
    # Now we test that we can get different radial velocities based on field
    # parameters.

    # Get the first sphere
    ds = fake_random_ds(16, fields = ("density",
      "velocity_x", "velocity_y", "velocity_z"))
    sp0 = ds.sphere(ds.domain_center, 0.25)

    assert_equal(sp0.keys(), [])

    sp1 = sp0.clone()
    sp0["density"]
    assert_equal(sp0.keys(), (("gas","density"),))
    assert_equal(sp1.keys(), [])

    sp1["density"]

    assert_array_equal(sp0["density"], sp1["density"])
