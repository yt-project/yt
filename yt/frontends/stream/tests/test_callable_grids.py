import numpy as np

from yt import load_amr_grids, load_hdf5_file
from yt.testing import (
    _amr_grid_index,
    assert_almost_equal,
    assert_equal,
    requires_file,
    requires_module,
)

turb_vels = "UnigridData/turb_vels.h5"

_existing_fields = (
    "Bx",
    "By",
    "Bz",
    "Density",
    "MagneticEnergy",
    "Temperature",
    "turb_x-velocity",
    "turb_y-velocity",
    "turb_z-velocity",
    "x-velocity",
    "y-velocity",
    "z-velocity",
)


@requires_file(turb_vels)
@requires_module("h5py")
def test_load_hdf5_file():
    ds1 = load_hdf5_file(turb_vels)
    assert_equal(ds1.domain_dimensions, [256, 256, 256])
    for field_name in _existing_fields:
        assert ("stream", field_name) in ds1.field_list
    assert_equal(ds1.r[:]["ones"].size, 256 * 256 * 256)
    assert_equal(ds1.r[:]["Density"].size, 256 * 256 * 256)
    # Now we test that we get the same results regardless of our decomp
    ds2 = load_hdf5_file(turb_vels, nchunks=19)
    assert_equal(ds2.domain_dimensions, [256, 256, 256])
    assert_equal(ds2.r[:]["ones"].size, 256 * 256 * 256)
    assert_equal(ds2.r[:]["Density"].size, 256 * 256 * 256)
    assert_almost_equal(ds2.r[:]["Density"].min(), ds1.r[:]["Density"].min())
    assert_almost_equal(ds2.r[:]["Density"].max(), ds1.r[:]["Density"].max())
    assert_almost_equal(ds2.r[:]["Density"].std(), ds1.r[:]["Density"].std())


_x_coefficients = (100, 50, 30, 10, 20)
_y_coefficients = (20, 90, 80, 30, 30)
_z_coefficients = (50, 10, 90, 40, 40)


def _grid_data_function(grid, field_name):
    # We want N points from the cell-center to the cell-center on the other side
    x, y, z = (
        np.linspace(
            grid.LeftEdge[i] + grid.dds[i] / 2,
            grid.RightEdge[i] - grid.dds[i] / 2,
            grid.ActiveDimensions[i],
        )
        for i in (0, 1, 2)
    )
    r = np.sqrt(
        ((x.d - 0.5) ** 2)[:, None, None]
        + ((y.d - 0.5) ** 2)[None, :, None]
        + ((z.d - 0.5) ** 2)[None, None, :]
    )
    atten = np.exp(-20 * (1.1 * r**2))
    xv = sum(
        c * np.sin(2 ** (1 + i) * (x.d * np.pi * 2))
        for i, c in enumerate(_x_coefficients)
    )
    yv = sum(
        c * np.sin(2 ** (1 + i) * (y.d * np.pi * 2))
        for i, c in enumerate(_y_coefficients)
    )
    zv = sum(
        c * np.sin(2 ** (1 + i) * (z.d * np.pi * 2))
        for i, c in enumerate(_z_coefficients)
    )
    return atten * (xv[:, None, None] * yv[None, :, None] * zv[None, None, :])


def test_load_callable():
    grid_data = []
    for level, le, re, dims in _amr_grid_index:
        grid_data.append(
            {
                "level": level,
                "left_edge": le,
                "right_edge": re,
                "dimensions": dims,
                "density": _grid_data_function,
            }
        )
    ds = load_amr_grids(
        grid_data, [32, 32, 32], bbox=np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    )
    assert_equal(ds.r[:].sum("cell_volume"), ds.domain_width.prod())
    assert_almost_equal(ds.r[:].max("density").d, 2660218.62833899)
    assert_almost_equal(ds.r[:].min("density").d, -2660218.62833899)
