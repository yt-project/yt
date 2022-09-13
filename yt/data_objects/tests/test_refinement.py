import numpy as np

import yt
from yt.testing import assert_array_equal, assert_equal


def setup_fake_refby():
    refine_by = np.array([5, 1, 1])
    top_grid_dim = [100, 10, 2]
    n1 = 100
    n2 = 10
    n3 = 2

    grid_data = [
        dict(
            left_edge=[0.0, 0.0, 0.0],
            right_edge=[1.0, np.pi, np.pi * 2.0],
            level=0,
            dimensions=np.array([n1, n2, n3]),
        ),
        dict(
            left_edge=[0.0, 0.0, 0.0],
            right_edge=[0.5, np.pi, np.pi * 2.0],
            level=1,
            dimensions=refine_by * [n1 / 2.0, n2, n3],
        ),
    ]

    for g in grid_data:
        g["density"] = (np.random.random(g["dimensions"].astype("i8")), "g/cm**3")
    bbox = np.array([[0.0, 1.0], [0.0, np.pi], [0.0, np.pi * 2]])

    ds = yt.load_amr_grids(
        grid_data,
        top_grid_dim,
        bbox=bbox,
        geometry="spherical",
        refine_by=refine_by,
        length_unit="kpc",
    )
    return ds


def test_refine_by():
    ds = setup_fake_refby()
    dd = ds.all_data()
    # This checks that we always refine_by 1 in dimensions 2 and 3
    dims = ds.domain_dimensions * ds.refine_by**ds.max_level
    for i in range(1, 3):
        # Check the refine_by == 1
        ncoords = np.unique(dd.icoords[:, i]).size
        assert_equal(ncoords, dims[i])
    for g in ds.index.grids:
        dims = ds.domain_dimensions * ds.refine_by**g.Level
        # Now we can check converting back to the reference space
        v = ((g.icoords + 1) / dims.astype("f8")).max(axis=0)
        v *= ds.domain_width
        assert_array_equal(v, g.RightEdge.d)
