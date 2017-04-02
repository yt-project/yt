from yt.testing import \
    assert_array_equal, \
    fake_random_ds
from yt.units import cm

def test_box_creation():

    # test that creating a region with left and right edge
    # with units works
    ds = fake_random_ds(32, length_unit=2)
    reg = ds.box([0, 0, 0]*cm, [2, 2, 2]*cm)
    dens_units = reg['density']

    reg = ds.box([0, 0, 0], [1, 1, 1])
    dens_no_units = reg['density']

    assert_array_equal(dens_units, dens_no_units)
