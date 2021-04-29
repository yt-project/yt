import numpy as np

from yt import SlicePlot
from yt.frontends.stream.data_structures import hexahedral_connectivity
from yt.loaders import load_hexahedral_mesh
from yt.testing import assert_almost_equal, assert_equal

# Field information


def test_stream_hexahedral():
    np.random.seed(0x4D3D3D3)
    Nx, Ny, Nz = 32, 18, 24
    # Note what we're doing here -- we are creating a randomly spaced mesh, but
    # because of how the accumulate operation works, we also reset the leftmost
    # cell boundary to 0.0.
    cell_x = np.random.random(Nx + 1)
    cell_x /= cell_x.sum()
    cell_x = np.add.accumulate(cell_x)
    cell_x[0] = 0.0

    cell_y = np.random.random(Ny + 1)
    cell_y /= cell_y.sum()
    cell_y = np.add.accumulate(cell_y)
    cell_y[0] = 0.0

    cell_z = np.random.random(Nz + 1)
    cell_z /= cell_z.sum()
    cell_z = np.add.accumulate(cell_z)
    cell_z[0] = 0.0

    coords, conn = hexahedral_connectivity(cell_x, cell_y, cell_z)
    data = {"random_field": np.random.random((Nx, Ny, Nz))}
    bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    ds = load_hexahedral_mesh(data, conn, coords, bbox=bbox)
    dd = ds.all_data()
    # raise RuntimeError
    assert_almost_equal(float(dd[("gas", "cell_volume")].sum(dtype="float64")), 1.0)
    assert_equal(dd[("index", "ones")].size, Nx * Ny * Nz)
    # Now we try it with a standard mesh
    cell_x = np.linspace(0.0, 1.0, Nx + 1)
    cell_y = np.linspace(0.0, 1.0, Ny + 1)
    cell_z = np.linspace(0.0, 1.0, Nz + 1)
    coords, conn = hexahedral_connectivity(cell_x, cell_y, cell_z)
    data = {"random_field": np.random.random((Nx, Ny, Nz))}
    bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    ds = load_hexahedral_mesh(data, conn, coords, bbox=bbox)
    dd = ds.all_data()
    assert_almost_equal(float(dd[("gas", "cell_volume")].sum(dtype="float64")), 1.0)
    assert_equal(dd[("index", "ones")].size, Nx * Ny * Nz)
    assert_almost_equal(dd[("index", "dx")].to_ndarray(), 1.0 / Nx)
    assert_almost_equal(dd[("index", "dy")].to_ndarray(), 1.0 / Ny)
    assert_almost_equal(dd[("index", "dz")].to_ndarray(), 1.0 / Nz)

    s = SlicePlot(ds, "x", "random_field")
    s._setup_plots()
    s.frb[("stream", "random_field")]
