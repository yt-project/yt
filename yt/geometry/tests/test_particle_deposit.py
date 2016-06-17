from yt.utilities.exceptions import \
    YTBoundsDefinitionError

from yt.testing import \
    fake_random_ds
from numpy.testing import \
    assert_raises

def test_cic_deposit():
    ds = fake_random_ds(64, nprocs = 8, particles=64**3)
    my_reg = ds.arbitrary_grid(ds.domain_left_edge, ds.domain_right_edge,
            dims=[1, 800, 800])
    f = ("deposit", "all_cic")
    assert_raises(YTBoundsDefinitionError, my_reg.__getitem__, f)
