import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

import yt.utilities.linear_interpolators as lin
from yt.testing import fake_random_ds
from yt.utilities.lib.interpolators import ghost_zone_interpolate


def test_linear_interpolator_1d():
    random_data = np.random.random(64)
    fv = {"x": np.mgrid[0.0:1.0:64j]}
    # evenly spaced bins
    ufi = lin.UnilinearFieldInterpolator(random_data, (0.0, 1.0), "x", truncate=True)
    assert_array_equal(ufi(fv), random_data)

    # randomly spaced bins
    size = 64
    shift = (1.0 / size) * np.random.random(size) - (0.5 / size)
    fv["x"] += shift
    ufi = lin.UnilinearFieldInterpolator(
        random_data, np.linspace(0.0, 1.0, size) + shift, "x", truncate=True
    )
    assert_array_almost_equal(ufi(fv), random_data, 15)


def test_linear_interpolator_2d():
    random_data = np.random.random((64, 64))
    # evenly spaced bins
    fv = dict(zip("xy", np.mgrid[0.0:1.0:64j, 0.0:1.0:64j], strict=True))
    fv = dict(zip("xyz", np.mgrid[0.0:1.0:64j, 0.0:1.0:64j], strict=False))
    bfi = lin.BilinearFieldInterpolator(
        random_data, (0.0, 1.0, 0.0, 1.0), "xy", truncate=True
    )
    assert_array_equal(bfi(fv), random_data)

    # randomly spaced bins
    size = 64
    bins = np.linspace(0.0, 1.0, size)
    shifts = {ax: (1.0 / size) * np.random.random(size) - (0.5 / size) for ax in "xy"}
    fv["x"] += shifts["x"][:, np.newaxis]
    fv["y"] += shifts["y"]
    bfi = lin.BilinearFieldInterpolator(
        random_data, (bins + shifts["x"], bins + shifts["y"]), "xy", truncate=True
    )
    assert_array_almost_equal(bfi(fv), random_data, 15)


def test_linear_interpolator_3d():
    random_data = np.random.random((64, 64, 64))
    # evenly spaced bins
    fv = dict(zip("xyz", np.mgrid[0.0:1.0:64j, 0.0:1.0:64j, 0.0:1.0:64j], strict=True))
    tfi = lin.TrilinearFieldInterpolator(
        random_data, (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), "xyz", truncate=True
    )
    assert_array_almost_equal(tfi(fv), random_data)

    # randomly spaced bins
    size = 64
    bins = np.linspace(0.0, 1.0, size)
    shifts = {ax: (1.0 / size) * np.random.random(size) - (0.5 / size) for ax in "xyz"}
    fv["x"] += shifts["x"][:, np.newaxis, np.newaxis]
    fv["y"] += shifts["y"][:, np.newaxis]
    fv["z"] += shifts["z"]
    tfi = lin.TrilinearFieldInterpolator(
        random_data,
        (bins + shifts["x"], bins + shifts["y"], bins + shifts["z"]),
        "xyz",
        True,
    )
    assert_array_almost_equal(tfi(fv), random_data, 15)


def test_linear_interpolator_4d():
    size = 32
    random_data = np.random.random((size,) * 4)
    # evenly spaced bins
    step = complex(0, size)
    fv = dict(
        zip(
            "xyzw",
            np.mgrid[0.0:1.0:step, 0.0:1.0:step, 0.0:1.0:step, 0.0:1.0:step],
            strict=True,
        )
    )
    tfi = lin.QuadrilinearFieldInterpolator(
        random_data, (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0), "xyzw", truncate=True
    )
    assert_array_almost_equal(tfi(fv), random_data)

    # randomly spaced bins
    bins = np.linspace(0.0, 1.0, size)
    shifts = {ax: (1.0 / size) * np.random.random(size) - (0.5 / size) for ax in "xyzw"}
    fv["x"] += shifts["x"][:, np.newaxis, np.newaxis, np.newaxis]
    fv["y"] += shifts["y"][:, np.newaxis, np.newaxis]
    fv["z"] += shifts["z"][:, np.newaxis]
    fv["w"] += shifts["w"]
    tfi = lin.QuadrilinearFieldInterpolator(
        random_data,
        (
            bins + shifts["x"],
            bins + shifts["y"],
            bins + shifts["z"],
            bins + shifts["w"],
        ),
        "xyzw",
        truncate=True,
    )
    assert_array_almost_equal(tfi(fv), random_data, 15)


def test_ghost_zone_extrapolation():
    ds = fake_random_ds(16)

    g = ds.index.grids[0]
    vec = g.get_vertex_centered_data(
        [("index", "x"), ("index", "y"), ("index", "z")], no_ghost=True
    )
    for i, ax in enumerate("xyz"):
        xc = g["index", ax]

        tf = lin.TrilinearFieldInterpolator(
            xc,
            (
                g.LeftEdge[0] + g.dds[0] / 2.0,
                g.RightEdge[0] - g.dds[0] / 2.0,
                g.LeftEdge[1] + g.dds[1] / 2.0,
                g.RightEdge[1] - g.dds[1] / 2.0,
                g.LeftEdge[2] + g.dds[2] / 2.0,
                g.RightEdge[2] - g.dds[2] / 2.0,
            ),
            ["x", "y", "z"],
            truncate=True,
        )

        lx, ly, lz = np.mgrid[
            g.LeftEdge[0] : g.RightEdge[0] : (g.ActiveDimensions[0] + 1) * 1j,
            g.LeftEdge[1] : g.RightEdge[1] : (g.ActiveDimensions[1] + 1) * 1j,
            g.LeftEdge[2] : g.RightEdge[2] : (g.ActiveDimensions[2] + 1) * 1j,
        ]
        xi = tf({"x": lx, "y": ly, "z": lz})

        xz = np.zeros(g.ActiveDimensions + 1)
        ghost_zone_interpolate(
            1,
            xc,
            np.array([0.5, 0.5, 0.5], dtype="f8"),
            xz,
            np.array([0.0, 0.0, 0.0], dtype="f8"),
        )

        ii = (lx, ly, lz)[i]
        assert_array_equal(ii, vec["index", ax])
        assert_array_equal(ii, xi)
        assert_array_equal(ii, xz)


def test_get_vertex_centered_data():
    ds = fake_random_ds(16)
    g = ds.index.grids[0]
    g.get_vertex_centered_data([("gas", "density")], no_ghost=True)


_lin_interpolators_by_dim = {
    1: lin.UnilinearFieldInterpolator,
    2: lin.BilinearFieldInterpolator,
    3: lin.TrilinearFieldInterpolator,
    4: lin.QuadrilinearFieldInterpolator,
}


@pytest.mark.parametrize("ndim", list(_lin_interpolators_by_dim.keys()))
def test_table_override(ndim):
    sz = 8

    random_data = np.random.random((sz,) * ndim)

    field_names = "xyzw"[:ndim]
    slc = slice(0.0, 1.0, complex(0, sz))
    fv = dict(zip(field_names, np.mgrid[(slc,) * ndim], strict=False))
    boundaries = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)[: ndim * 2]

    interp_class = _lin_interpolators_by_dim[ndim]
    # get interpolator with default behavior (a copy of the table is stored)
    interpolator = interp_class(random_data, boundaries, field_names, truncate=True)
    # evaluate at table grid nodes (should return the input table exactly)
    assert_array_almost_equal(interpolator(fv), random_data)
    # check that we can change the table used after interpolator initialization
    table_2 = random_data * 2
    assert_array_almost_equal(interpolator(fv, table=table_2), table_2)

    # check that we can use an interpolator without storing the initial table
    interpolator = interp_class(
        random_data, boundaries, field_names, truncate=True, store_table=False
    )
    assert_array_almost_equal(interpolator(fv, table=table_2), table_2)

    with pytest.raises(
        ValueError, match="You must either store the table used when initializing"
    ):
        interpolator(fv)


@pytest.mark.parametrize("ndim", list(_lin_interpolators_by_dim.keys()))
def test_bin_validation(ndim):
    interp_class = _lin_interpolators_by_dim[ndim]

    sz = 8
    random_data = np.random.random((sz,) * ndim)
    field_names = "xyzw"[:ndim]

    bad_bounds = np.linspace(0.0, 1.0, sz - 1)
    good_bounds = np.linspace(0.0, 1.0, sz)
    if ndim == 1:
        bounds = bad_bounds
    else:
        bounds = [
            good_bounds,
        ] * ndim
        bounds[0] = bad_bounds
        bounds = tuple(bounds)

    with pytest.raises(ValueError, match=f"{field_names[0]} bins array not"):
        interp_class(random_data, bounds, field_names)
