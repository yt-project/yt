from yt.testing import assert_equal, fake_amr_ds, fake_particle_ds, fake_random_ds


def test_center_squeeze():
    # checks that the center is reshaped correctly

    # create and test amr, random and particle data
    check_single_ds(fake_amr_ds(fields=("Density",), units=("g/cm**3",)))
    check_single_ds(fake_random_ds(16, fields=("Density",), units=("g/cm**3",)))
    check_single_ds(fake_particle_ds(npart=100))


def check_single_ds(ds):
    # checks that the center
    center = ds.domain_center  # reference center value
    for test_shape in [(1, 3), (1, 1, 3)]:
        new_center = center.reshape(test_shape)
        assert_equal(ds.sphere(new_center, 0.25).center, center)
        assert_equal(ds.slice(0, 0.25, center=new_center).center, center)
        assert_equal(
            ds.region(new_center, [-0.25, -0.25, -0.25], [0.25, 0.25, 0.25]).center,
            center,
        )
