from yt.testing import assert_array_equal, assert_equal, fake_random_ds


def test_clone_sphere():
    # Now we test that we can get different radial velocities based on field
    # parameters.
    fields = ("density", "velocity_x", "velocity_y", "velocity_z")
    units = ("g/cm**3", "cm/s", "cm/s", "cm/s")
    # Get the first sphere
    ds = fake_random_ds(16, fields=fields, units=units)
    sp0 = ds.sphere(ds.domain_center, 0.25)

    assert_equal(list(sp0.keys()), [])

    sp1 = sp0.clone()
    sp0[("gas", "density")]
    assert_equal(list(sp0.keys()), (("gas", "density"),))
    assert_equal(list(sp1.keys()), [])

    sp1[("gas", "density")]

    assert_array_equal(sp0[("gas", "density")], sp1[("gas", "density")])


def test_clone_cut_region():
    fields = ("density", "temperature")
    units = ("g/cm**3", "K")
    ds = fake_random_ds(64, nprocs=4, fields=fields, units=units)
    dd = ds.all_data()
    reg1 = dd.cut_region(
        ["obj[('gas', 'temperature')] > 0.5", "obj[('gas', 'density')] < 0.75"]
    )
    reg2 = reg1.clone()
    assert_array_equal(reg1[("gas", "density")], reg2[("gas", "density")])
