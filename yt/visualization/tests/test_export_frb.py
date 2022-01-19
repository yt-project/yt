import numpy as np

from yt.testing import assert_allclose_units, assert_equal, fake_random_ds


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_export_frb():
    test_ds = fake_random_ds(128)
    slc = test_ds.slice(0, 0.5)
    frb = slc.to_frb((0.5, "unitary"), 64)
    frb_ds = frb.export_dataset(fields=[("gas", "density")], nprocs=8)
    dd_frb = frb_ds.all_data()

    assert_equal(frb_ds.domain_left_edge.v, np.array([0.25, 0.25, 0.0]))
    assert_equal(frb_ds.domain_right_edge.v, np.array([0.75, 0.75, 1.0]))
    assert_equal(frb_ds.domain_width.v, np.array([0.5, 0.5, 1.0]))
    assert_equal(frb_ds.domain_dimensions, np.array([64, 64, 1], dtype="int64"))
    assert_allclose_units(
        frb[("gas", "density")].sum(),
        dd_frb.quantities.total_quantity(("gas", "density")),
    )
    assert_equal(frb_ds.index.num_grids, 8)
