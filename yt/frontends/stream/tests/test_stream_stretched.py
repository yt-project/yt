import numpy as np

from yt import load_uniform_grid
from yt.testing import assert_almost_equal, assert_equal


def test_variable_dx():
    np.random.seed(0x4D3D3D3)

    data = {"density": np.random.random((128, 128, 128))}

    cell_widths = []
    for _ in range(3):
        cw = np.random.random(128)
        cw /= cw.sum()
        cell_widths.append(cw)

    ds = load_uniform_grid(
        data,
        [128, 128, 128],
        bbox=np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]),
        cell_widths=cell_widths,
    )

    # We now check that we get all of our original cell widths back out, and
    # only those cell widths

    assert_equal(np.unique(ds.index.grids[0]["index", "dx"]).size, 128)
    assert_equal(ds.index.grids[0]["index", "dx"][:, 0, 0], cell_widths[0])

    assert_equal(np.unique(ds.index.grids[0]["index", "dx"]).size, 128)
    assert_equal(ds.index.grids[0]["index", "dy"][0, :, 0], cell_widths[1])

    assert_equal(np.unique(ds.index.grids[0]["index", "dx"]).size, 128)
    assert_equal(ds.index.grids[0]["index", "dz"][0, 0, :], cell_widths[2])

    assert_equal(np.unique(ds.index.grids[0]["index", "x"]).size, 128)
    center_x = np.add.accumulate(cell_widths[0]) - 0.5 * cell_widths[0]
    assert_equal(center_x, ds.index.grids[0]["index", "x"][:, 0, 0])

    assert_equal(np.unique(ds.index.grids[0]["index", "y"]).size, 128)
    center_y = np.add.accumulate(cell_widths[1]) - 0.5 * cell_widths[1]
    assert_equal(center_y, ds.index.grids[0]["index", "y"][0, :, 0])

    assert_equal(np.unique(ds.index.grids[0]["index", "z"]).size, 128)
    center_z = np.add.accumulate(cell_widths[2]) - 0.5 * cell_widths[2]
    assert_equal(center_z, ds.index.grids[0]["index", "z"][0, 0, :])

    assert_almost_equal(ds.r[:].sum(("index", "cell_volume")), ds.domain_width.prod())

    for ax in "xyz":
        dd = ds.all_data()
        p = dd.integrate("ones", axis=ax)
        assert_almost_equal(p["index", "ones"].max().d, 1.0)
        assert_almost_equal(p["index", "ones"].min().d, 1.0)
