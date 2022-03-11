import numpy as np

from yt import load_uniform_grid
from yt.testing import assert_equal


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
