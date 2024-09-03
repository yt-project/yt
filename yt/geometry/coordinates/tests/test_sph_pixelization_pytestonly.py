import numpy as np
import pytest
import unyt

import yt
from yt.data_objects.selection_objects.region import YTRegion
from yt.testing import (
    assert_rel_equal,
    cubicspline_python,
    distancematrix,
    fake_random_sph_ds,
    fake_sph_flexible_grid_ds,
    integrate_kernel,
)


@pytest.mark.parametrize("weighted", [True, False])
@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("depth", [None, (1.0, "cm")])
@pytest.mark.parametrize("shiftcenter", [False, True])
@pytest.mark.parametrize("axis", [0, 1, 2])
def test_sph_proj_general_alongaxes(
    axis: int,
    shiftcenter: bool,
    depth: float | None,
    periodic: bool,
    weighted: bool,
) -> None:
    """
    The previous projection tests were for a specific issue.
    Here, we test more functionality of the projections.
    We just send lines of sight through pixel centers for convenience.
    Particles at [0.5, 1.5, 2.5] (in each coordinate)
    smoothing lengths 0.25
    all particles have mass 1., density 1.5,
    except the single center particle, with mass 2., density 3.

    Parameters:
    -----------
    axis: {0, 1, 2}
        projection axis (aligned with sim. axis)
    shiftcenter: bool
        shift the coordinates to center the projection on.
        (The grid is offset to this same center)
    depth: float or None
        depth of the projection slice
    periodic: bool
        assume periodic boundary conditions, or not
    weighted: bool
        make a weighted projection (density-weighted density), or not

    Returns:
    --------
    None
    """
    if shiftcenter:
        center = unyt.unyt_array(np.array((0.625, 0.625, 0.625)), "cm")
    else:
        center = unyt.unyt_array(np.array((1.5, 1.5, 1.5)), "cm")
    bbox = unyt.unyt_array(np.array([[0.0, 3.0], [0.0, 3.0], [0.0, 3.0]]), "cm")
    hsml_factor = 0.5
    unitrho = 1.5

    # test correct centering, particle selection
    def makemasses(i, j, k):
        if i == j == k == 1:
            return 2.0
        else:
            return 1.0

    # m / rho, factor 1. / hsml**2 is included in the kernel integral
    # (density is adjusted, so same for center particle)
    prefactor = 1.0 / unitrho  # / (0.5 * 0.5)**2
    dl_cen = integrate_kernel(cubicspline_python, 0.0, 0.25)

    # result shouldn't depend explicitly on the center if we re-center
    # the data, unless we get cut-offs in the non-periodic case
    ds = fake_sph_flexible_grid_ds(
        hsml_factor=hsml_factor,
        nperside=3,
        periodic=periodic,
        offsets=np.full(3, 0.5),
        massgenerator=makemasses,
        unitrho=unitrho,
        bbox=bbox.v,
        recenter=center.v,
    )
    if depth is None:
        source = ds.all_data()
    else:
        depth = unyt.unyt_quantity(*depth)
        le = np.array(ds.domain_left_edge)
        re = np.array(ds.domain_right_edge)
        le[axis] = center[axis] - 0.5 * depth
        re[axis] = center[axis] + 0.5 * depth
        cen = 0.5 * (le + re)
        reg = YTRegion(center=cen, left_edge=le, right_edge=re, ds=ds)
        source = reg

    # we don't actually want a plot, it's just a straightforward,
    # common way to get an frb / image array
    if weighted:
        toweight_field = ("gas", "density")
    else:
        toweight_field = None
    prj = yt.ProjectionPlot(
        ds,
        axis,
        ("gas", "density"),
        width=(2.5, "cm"),
        weight_field=toweight_field,
        buff_size=(5, 5),
        center=center,
        data_source=source,
    )
    img = prj.frb.data[("gas", "density")]
    if weighted:
        expected_out = np.zeros(
            (
                5,
                5,
            ),
            dtype=img.v.dtype,
        )
        expected_out[::2, ::2] = unitrho
        if depth is None:
            ## during shift, particle coords do wrap around edges
            # if (not periodic) and shiftcenter:
            #    # weight 1. for unitrho, 2. for 2. * untrho
            #    expected_out[2, 2] *= 5. / 3.
            # else:
            # weight (2 * 1.) for unitrho, (1 * 2.) for 2. * unitrho
            expected_out[2, 2] *= 1.5
        else:
            # only 2 * unitrho element included
            expected_out[2, 2] *= 2.0
    else:
        expected_out = np.zeros(
            (
                5,
                5,
            ),
            dtype=img.v.dtype,
        )
        expected_out[::2, ::2] = dl_cen * prefactor * unitrho
        if depth is None:
            # 3 particles per l.o.s., including the denser one
            expected_out *= 3.0
            expected_out[2, 2] *= 4.0 / 3.0
        else:
            # 1 particle per l.o.s., including the denser one
            expected_out[2, 2] *= 2.0
    # grid is shifted to the left -> 'missing' stuff at the left
    if (not periodic) and shiftcenter:
        expected_out[:1, :] = 0.0
        expected_out[:, :1] = 0.0
    # print(axis, shiftcenter, depth, periodic, weighted)
    # print(expected_out)
    # print(img.v)
    assert_rel_equal(expected_out, img.v, 5)


@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("shiftcenter", [False, True])
@pytest.mark.parametrize("zoff", [0.0, 0.1, 0.5, 1.0])
@pytest.mark.parametrize("axis", [0, 1, 2])
def test_sph_slice_general_alongaxes(
    axis: int,
    shiftcenter: bool,
    periodic: bool,
    zoff: float,
) -> None:
    """
    Particles at [0.5, 1.5, 2.5] (in each coordinate)
    smoothing lengths 0.25
    all particles have mass 1., density 1.5,
    except the single center particle, with mass 2., density 3.

    Parameters:
    -----------
    axis: {0, 1, 2}
        projection axis (aligned with sim. axis)
    northvector: tuple
        y-axis direction in the final plot (direction vector)
    shiftcenter: bool
        shift the coordinates to center the projection on.
        (The grid is offset to this same center)
    zoff: float
        offset of the slice plane from the SPH particle center plane
    periodic: bool
        assume periodic boundary conditions, or not

    Returns:
    --------
    None
    """
    if shiftcenter:
        center = unyt.unyt_array(np.array((0.625, 0.625, 0.625)), "cm")
    else:
        center = unyt.unyt_array(np.array((1.5, 1.5, 1.5)), "cm")
    bbox = unyt.unyt_array(np.array([[0.0, 3.0], [0.0, 3.0], [0.0, 3.0]]), "cm")
    hsml_factor = 0.5
    unitrho = 1.5

    # test correct centering, particle selection
    def makemasses(i, j, k):
        if i == j == k == 1:
            return 2.0
        elif i == j == k == 2:
            return 3.0
        else:
            return 1.0

    # result shouldn't depend explicitly on the center if we re-center
    # the data, unless we get cut-offs in the non-periodic case
    ds = fake_sph_flexible_grid_ds(
        hsml_factor=hsml_factor,
        nperside=3,
        periodic=periodic,
        offsets=np.full(3, 0.5),
        massgenerator=makemasses,
        unitrho=unitrho,
        bbox=bbox.v,
        recenter=center.v,
    )
    ad = ds.all_data()
    # print(ad[('gas', 'position')])
    outgridsize = 10
    width = 2.5
    _center = center.to("cm").v.copy()
    _center[axis] += zoff

    # we don't actually want a plot, it's just a straightforward,
    # common way to get an frb / image array
    slc = yt.SlicePlot(
        ds,
        axis,
        ("gas", "density"),
        width=(width, "cm"),
        buff_size=(outgridsize,) * 2,
        center=(_center, "cm"),
    )
    img = slc.frb.data[("gas", "density")]

    # center is same in non-projection coords
    if axis == 0:
        ci = 1
    else:
        ci = 0
    gridcens = (
        _center[ci]
        - 0.5 * width
        + 0.5 * width / outgridsize
        + np.arange(outgridsize) * width / outgridsize
    )
    xgrid = np.repeat(gridcens, outgridsize)
    ygrid = np.tile(gridcens, outgridsize)
    zgrid = np.full(outgridsize**2, _center[axis])
    gridcoords = np.empty((outgridsize**2, 3), dtype=xgrid.dtype)
    if axis == 2:
        gridcoords[:, 0] = xgrid
        gridcoords[:, 1] = ygrid
        gridcoords[:, 2] = zgrid
    elif axis == 0:
        gridcoords[:, 0] = zgrid
        gridcoords[:, 1] = xgrid
        gridcoords[:, 2] = ygrid
    elif axis == 1:
        gridcoords[:, 0] = ygrid
        gridcoords[:, 1] = zgrid
        gridcoords[:, 2] = xgrid
    ad = ds.all_data()
    sphcoords = np.array(
        [
            (ad[("gas", "x")]).to("cm"),
            (ad[("gas", "y")]).to("cm"),
            (ad[("gas", "z")]).to("cm"),
        ]
    ).T
    # print("sphcoords:")
    # print(sphcoords)
    # print("gridcoords:")
    # print(gridcoords)
    dists = distancematrix(
        gridcoords,
        sphcoords,
        periodic=(periodic,) * 3,
        periods=np.array([3.0, 3.0, 3.0]),
    )
    # print("dists <= 1:")
    # print(dists <= 1)
    sml = (ad[("gas", "smoothing_length")]).to("cm")
    normkern = cubicspline_python(dists / sml.v[np.newaxis, :])
    sphcontr = normkern / sml[np.newaxis, :] ** 3 * ad[("gas", "mass")]
    contsum = np.sum(sphcontr, axis=1)
    sphweights = (
        normkern
        / sml[np.newaxis, :] ** 3
        * ad[("gas", "mass")]
        / ad[("gas", "density")]
    )
    weights = np.sum(sphweights, axis=1)
    nzeromask = np.logical_not(weights == 0)
    expected = np.zeros(weights.shape, weights.dtype)
    expected[nzeromask] = contsum[nzeromask] / weights[nzeromask]
    expected = expected.reshape((outgridsize, outgridsize))
    # expected[np.isnan(expected)] = 0.0  # convention in the slices

    # print("expected:\n", expected)
    # print("recovered:\n", img.v)
    assert_rel_equal(expected, img.v, 5)


@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("shiftcenter", [False, True])
@pytest.mark.parametrize("northvector", [None, (1.0e-4, 1.0, 0.0)])
@pytest.mark.parametrize("zoff", [0.0, 0.1, 0.5, 1.0])
def test_sph_slice_general_offaxis(
    northvector: tuple[float, float, float] | None,
    shiftcenter: bool,
    zoff: float,
    periodic: bool,
) -> None:
    """
    Same as the on-axis slices, but we rotate the basis vectors
    to test whether roations are handled ok. the rotation is chosen to
    be small so that in/exclusion of particles within bboxes, etc.
    works out the same way.
    Particles at [0.5, 1.5, 2.5] (in each coordinate)
    smoothing lengths 0.25
    all particles have mass 1., density 1.5,
    except the single center particle, with mass 2., density 3.

    Parameters:
    -----------
    northvector: tuple
        y-axis direction in the final plot (direction vector)
    shiftcenter: bool
        shift the coordinates to center the projection on.
        (The grid is offset to this same center)
    zoff: float
        offset of the slice plane from the SPH particle center plane
    periodic: bool
        assume periodic boundary conditions, or not

    Returns:
    --------
    None
    """
    if shiftcenter:
        center = np.array((0.625, 0.625, 0.625))  # cm
    else:
        center = np.array((1.5, 1.5, 1.5))  # cm
    bbox = unyt.unyt_array(np.array([[0.0, 3.0], [0.0, 3.0], [0.0, 3.0]]), "cm")
    hsml_factor = 0.5
    unitrho = 1.5

    # test correct centering, particle selection
    def makemasses(i, j, k):
        if i == j == k == 1:
            return 2.0
        else:
            return 1.0

    # try to make sure dl differences from periodic wrapping are small
    epsilon = 1e-4
    projaxis = np.array([epsilon, 0.00, np.sqrt(1.0 - epsilon**2)])
    e1dir = projaxis / np.sqrt(np.sum(projaxis**2))
    if northvector is None:
        e2dir = np.array([0.0, 1.0, 0.0])
    else:
        e2dir = np.asarray(northvector)
    e2dir = e2dir - np.sum(e1dir * e2dir) * e2dir  # orthonormalize
    e2dir /= np.sqrt(np.sum(e2dir**2))
    e3dir = np.cross(e2dir, e1dir)

    outgridsize = 10
    width = 2.5
    _center = center.copy()
    _center += zoff * e1dir

    ds = fake_sph_flexible_grid_ds(
        hsml_factor=hsml_factor,
        nperside=3,
        periodic=periodic,
        offsets=np.full(3, 0.5),
        massgenerator=makemasses,
        unitrho=unitrho,
        bbox=bbox.v,
        recenter=center,
        e1hat=e1dir,
        e2hat=e2dir,
        e3hat=e3dir,
    )

    # source = ds.all_data()
    # couple to dataset -> right unit registry
    center = ds.arr(center, "cm")
    # print("position:\n", source["gas", "position"])
    slc = yt.SlicePlot(
        ds,
        e1dir,
        ("gas", "density"),
        width=(width, "cm"),
        buff_size=(outgridsize,) * 2,
        center=(_center, "cm"),
        north_vector=e2dir,
    )
    img = slc.frb.data[("gas", "density")]

    # center is same in x/y (e3dir/e2dir)
    gridcenx = (
        np.dot(_center, e3dir)
        - 0.5 * width
        + 0.5 * width / outgridsize
        + np.arange(outgridsize) * width / outgridsize
    )
    gridceny = (
        np.dot(_center, e2dir)
        - 0.5 * width
        + 0.5 * width / outgridsize
        + np.arange(outgridsize) * width / outgridsize
    )
    xgrid = np.repeat(gridcenx, outgridsize)
    ygrid = np.tile(gridceny, outgridsize)
    zgrid = np.full(outgridsize**2, np.dot(_center, e1dir))
    gridcoords = (
        xgrid[:, np.newaxis] * e3dir[np.newaxis, :]
        + ygrid[:, np.newaxis] * e2dir[np.newaxis, :]
        + zgrid[:, np.newaxis] * e1dir[np.newaxis, :]
    )
    # print("gridcoords:")
    # print(gridcoords)
    ad = ds.all_data()
    sphcoords = np.array(
        [
            (ad[("gas", "x")]).to("cm"),
            (ad[("gas", "y")]).to("cm"),
            (ad[("gas", "z")]).to("cm"),
        ]
    ).T
    dists = distancematrix(
        gridcoords,
        sphcoords,
        periodic=(periodic,) * 3,
        periods=np.array([3.0, 3.0, 3.0]),
    )
    sml = (ad[("gas", "smoothing_length")]).to("cm")
    normkern = cubicspline_python(dists / sml.v[np.newaxis, :])
    sphcontr = normkern / sml[np.newaxis, :] ** 3 * ad[("gas", "mass")]
    contsum = np.sum(sphcontr, axis=1)
    sphweights = (
        normkern
        / sml[np.newaxis, :] ** 3
        * ad[("gas", "mass")]
        / ad[("gas", "density")]
    )
    weights = np.sum(sphweights, axis=1)
    nzeromask = np.logical_not(weights == 0)
    expected = np.zeros(weights.shape, weights.dtype)
    expected[nzeromask] = contsum[nzeromask] / weights[nzeromask]
    expected = expected.reshape((outgridsize, outgridsize))
    expected = expected.T  # transposed for image plotting
    # expected[np.isnan(expected)] = 0.0  # convention in the slices

    # print(axis, shiftcenter, depth, periodic, weighted)
    # print("expected:\n", expected)
    # print("recovered:\n", img.v)
    assert_rel_equal(expected, img.v, 4)


# only axis-aligned; testing YTArbitraryGrid, YTCoveringGrid
@pytest.mark.parametrize("periodic", [True, False, (True, True, False)])
@pytest.mark.parametrize("wholebox", [True, False])
def test_sph_grid(
    periodic: bool | tuple[bool, bool, bool],
    wholebox: bool,
) -> None:
    bbox = np.array([[-1.0, 3.0], [1.0, 5.2], [-1.0, 3.0]])
    ds = fake_random_sph_ds(50, bbox, periodic=periodic)

    if not hasattr(periodic, "__len__"):
        periodic = (periodic,) * 3

    if wholebox:
        left = bbox[:, 0].copy()
        level = 2
        ncells = np.array([2**level] * 3)
        # print("left: ", left)
        # print("ncells: ", ncells)
        resgrid = ds.covering_grid(level, tuple(left), ncells)
        right = bbox[:, 1].copy()
        xedges = np.linspace(left[0], right[0], ncells[0] + 1)
        yedges = np.linspace(left[1], right[1], ncells[1] + 1)
        zedges = np.linspace(left[2], right[2], ncells[2] + 1)
    else:
        left = np.array([-1.0, 1.8, -1.0])
        right = np.array([2.5, 5.2, 2.5])
        ncells = np.array([3, 4, 4])
        resgrid = ds.arbitrary_grid(left, right, dims=ncells)
        xedges = np.linspace(left[0], right[0], ncells[0] + 1)
        yedges = np.linspace(left[1], right[1], ncells[1] + 1)
        zedges = np.linspace(left[2], right[2], ncells[2] + 1)
    res = resgrid["gas", "density"]
    xcens = 0.5 * (xedges[:-1] + xedges[1:])
    ycens = 0.5 * (yedges[:-1] + yedges[1:])
    zcens = 0.5 * (zedges[:-1] + zedges[1:])

    ad = ds.all_data()
    sphcoords = np.array(
        [
            (ad[("gas", "x")]).to("cm"),
            (ad[("gas", "y")]).to("cm"),
            (ad[("gas", "z")]).to("cm"),
        ]
    ).T
    gridx, gridy, gridz = np.meshgrid(xcens, ycens, zcens, indexing="ij")
    outshape = gridx.shape
    gridx = gridx.flatten()
    gridy = gridy.flatten()
    gridz = gridz.flatten()
    gridcoords = np.array([gridx, gridy, gridz]).T
    periods = bbox[:, 1] - bbox[:, 0]
    dists = distancematrix(gridcoords, sphcoords, periodic=periodic, periods=periods)
    sml = (ad[("gas", "smoothing_length")]).to("cm")
    normkern = cubicspline_python(dists / sml.v[np.newaxis, :])
    sphcontr = normkern / sml[np.newaxis, :] ** 3 * ad[("gas", "mass")]
    contsum = np.sum(sphcontr, axis=1)
    sphweights = (
        normkern
        / sml[np.newaxis, :] ** 3
        * ad[("gas", "mass")]
        / ad[("gas", "density")]
    )
    weights = np.sum(sphweights, axis=1)
    nzeromask = np.logical_not(weights == 0)
    expected = np.zeros(weights.shape, weights.dtype)
    expected[nzeromask] = contsum[nzeromask] / weights[nzeromask]
    expected = expected.reshape(outshape)
    # expected[np.isnan(expected)] = 0.0  # convention in the slices

    # print(axis, shiftcenter, depth, periodic, weighted)
    # print("expected:\n", expected)
    # print("recovered:\n", res.v)
    assert_rel_equal(expected, res.v, 4)
