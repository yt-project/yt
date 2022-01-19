import numpy as np

from yt.loaders import load_uniform_grid
from yt.testing import assert_equal, fake_random_ds


def test_exclude_above():
    the_ds = fake_random_ds(ndims=3)
    all_data = the_ds.all_data()
    new_ds = all_data.exclude_above(("gas", "density"), 1)
    assert_equal(new_ds[("gas", "density")], all_data[("gas", "density")])
    new_ds = all_data.exclude_above(("gas", "density"), 1e6, "g/m**3")
    assert_equal(new_ds[("gas", "density")], all_data[("gas", "density")])
    new_ds = all_data.exclude_above(("gas", "density"), 0)
    assert_equal(new_ds[("gas", "density")], [])


def test_exclude_below():
    the_ds = fake_random_ds(ndims=3)
    all_data = the_ds.all_data()
    new_ds = all_data.exclude_below(("gas", "density"), 1)
    assert_equal(new_ds[("gas", "density")], [])
    new_ds = all_data.exclude_below(("gas", "density"), 1e6, "g/m**3")
    assert_equal(new_ds[("gas", "density")], [])
    new_ds = all_data.exclude_below(("gas", "density"), 0)
    assert_equal(new_ds[("gas", "density")], all_data[("gas", "density")])


def test_exclude_nan():
    test_array = np.nan * np.ones((10, 10, 10))
    test_array[1, 1, :] = 1
    data = dict(density=test_array)
    ds = load_uniform_grid(data, test_array.shape, length_unit="cm", nprocs=1)
    ad = ds.all_data()
    no_nan_ds = ad.exclude_nan(("gas", "density"))
    assert_equal(no_nan_ds[("gas", "density")], np.array(np.ones(10)))


def test_equal():
    test_array = np.ones((10, 10, 10))
    test_array[1, 1, :] = 2.0
    test_array[2, 1, :] = 3.0
    data = dict(density=test_array)
    ds = load_uniform_grid(data, test_array.shape, length_unit="cm", nprocs=1)
    ad = ds.all_data()
    no_ones = ad.exclude_equal(("gas", "density"), 1.0)
    assert np.all(no_ones[("gas", "density")] != 1.0)
    only_ones = ad.include_equal(("gas", "density"), 1.0)
    assert np.all(only_ones[("gas", "density")] == 1.0)


def test_inside_outside():
    test_array = np.ones((10, 10, 10))
    test_array[1, 1, :] = 2.0
    test_array[2, 1, :] = 3.0
    data = dict(density=test_array)
    ds = load_uniform_grid(data, test_array.shape, length_unit="cm", nprocs=1)
    ad = ds.all_data()

    only_ones_and_twos = ad.include_inside(("gas", "density"), 0.9, 2.1)
    assert np.all(only_ones_and_twos[("gas", "density")] != 3.0)
    assert len(only_ones_and_twos[("gas", "density")]) == 990

    only_ones_and_twos = ad.exclude_outside(("gas", "density"), 0.9, 2.1)
    assert len(only_ones_and_twos[("gas", "density")]) == 990
    assert np.all(only_ones_and_twos[("gas", "density")] != 3.0)

    only_threes = ad.include_outside(("gas", "density"), 0.9, 2.1)
    assert np.all(only_threes[("gas", "density")] == 3)
    assert len(only_threes[("gas", "density")]) == 10

    only_threes = ad.include_outside(("gas", "density"), 0.9, 2.1)
    assert np.all(only_threes[("gas", "density")] == 3)
    assert len(only_threes[("gas", "density")]) == 10

    # Repeat, but convert units to g/m**3
    only_ones_and_twos = ad.include_inside(("gas", "density"), 0.9e6, 2.1e6, "g/m**3")
    assert np.all(only_ones_and_twos[("gas", "density")] != 3.0)
    assert len(only_ones_and_twos[("gas", "density")]) == 990

    only_ones_and_twos = ad.exclude_outside(("gas", "density"), 0.9e6, 2.1e6, "g/m**3")
    assert len(only_ones_and_twos[("gas", "density")]) == 990
    assert np.all(only_ones_and_twos[("gas", "density")] != 3.0)

    only_threes = ad.include_outside(("gas", "density"), 0.9e6, 2.1e6, "g/m**3")
    assert np.all(only_threes[("gas", "density")] == 3)
    assert len(only_threes[("gas", "density")]) == 10

    only_threes = ad.include_outside(("gas", "density"), 0.9e6, 2.1e6, "g/m**3")
    assert np.all(only_threes[("gas", "density")] == 3)
    assert len(only_threes[("gas", "density")]) == 10
