import hashlib
import itertools as it
import os
import pickle
import shutil
import sys
import tempfile
import unittest
from functools import wraps
from importlib.util import find_spec
from shutil import which
from unittest import SkipTest

import matplotlib
import numpy as np
import pytest
from more_itertools import always_iterable
from numpy.random import RandomState
from unyt.exceptions import UnitOperationError

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.funcs import is_sequence
from yt.loaders import load
from yt.units.yt_array import YTArray, YTQuantity

# we import this in a weird way from numpy.testing to avoid triggering
# flake8 errors from the unused imports. These test functions are imported
# elsewhere in yt from here so we want them to be imported here.
from numpy.testing import assert_array_equal, assert_almost_equal  # NOQA isort:skip
from numpy.testing import assert_equal, assert_array_less  # NOQA isort:skip
from numpy.testing import assert_string_equal  # NOQA isort:skip
from numpy.testing import assert_array_almost_equal_nulp  # isort:skip
from numpy.testing import assert_allclose, assert_raises  # NOQA isort:skip
from numpy.testing import assert_approx_equal  # NOQA isort:skip
from numpy.testing import assert_array_almost_equal  # NOQA isort:skip

ANSWER_TEST_TAG = "answer_test"
# Expose assert_true and assert_less_equal from unittest.TestCase
# this is adopted from nose. Doing this here allows us to avoid importing
# nose at the top level.
class _Dummy(unittest.TestCase):
    def nop(self):
        pass


_t = _Dummy("nop")

assert_true = _t.assertTrue
assert_less_equal = _t.assertLessEqual


def assert_rel_equal(a1, a2, decimals, err_msg="", verbose=True):
    # We have nan checks in here because occasionally we have fields that get
    # weighted without non-zero weights.  I'm looking at you, particle fields!
    if isinstance(a1, np.ndarray):
        assert a1.size == a2.size
        # Mask out NaNs
        assert (np.isnan(a1) == np.isnan(a2)).all()
        a1[np.isnan(a1)] = 1.0
        a2[np.isnan(a2)] = 1.0
        # Mask out 0
        ind1 = np.array(np.abs(a1) < np.finfo(a1.dtype).eps)
        ind2 = np.array(np.abs(a2) < np.finfo(a2.dtype).eps)
        assert (ind1 == ind2).all()
        a1[ind1] = 1.0
        a2[ind2] = 1.0
    elif np.any(np.isnan(a1)) and np.any(np.isnan(a2)):
        return True
    if not isinstance(a1, np.ndarray) and a1 == a2 == 0.0:
        # NANS!
        a1 = a2 = 1.0
    return assert_almost_equal(
        np.array(a1) / np.array(a2), 1.0, decimals, err_msg=err_msg, verbose=verbose
    )


def amrspace(extent, levels=7, cells=8):
    """Creates two numpy arrays representing the left and right bounds of
    an AMR grid as well as an array for the AMR level of each cell.

    Parameters
    ----------
    extent : array-like
        This a sequence of length 2*ndims that is the bounds of each dimension.
        For example, the 2D unit square would be given by [0.0, 1.0, 0.0, 1.0].
        A 3D cylindrical grid may look like [0.0, 2.0, -1.0, 1.0, 0.0, 2*np.pi].
    levels : int or sequence of ints, optional
        This is the number of AMR refinement levels.  If given as a sequence (of
        length ndims), then each dimension will be refined down to this level.
        All values in this array must be the same or zero.  A zero valued dimension
        indicates that this dim should not be refined.  Taking the 3D cylindrical
        example above if we don't want refine theta but want r and z at 5 we would
        set levels=(5, 5, 0).
    cells : int, optional
        This is the number of cells per refinement level.

    Returns
    -------
    left : float ndarray, shape=(npoints, ndims)
        The left AMR grid points.
    right : float ndarray, shape=(npoints, ndims)
        The right AMR grid points.
    level : int ndarray, shape=(npoints,)
        The AMR level for each point.

    Examples
    --------
    >>> l, r, lvl = amrspace([0.0, 2.0, 1.0, 2.0, 0.0, 3.14], levels=(3, 3, 0), cells=2)
    >>> print(l)
    [[ 0.     1.     0.   ]
     [ 0.25   1.     0.   ]
     [ 0.     1.125  0.   ]
     [ 0.25   1.125  0.   ]
     [ 0.5    1.     0.   ]
     [ 0.     1.25   0.   ]
     [ 0.5    1.25   0.   ]
     [ 1.     1.     0.   ]
     [ 0.     1.5    0.   ]
     [ 1.     1.5    0.   ]]

    """
    extent = np.asarray(extent, dtype="f8")
    dextent = extent[1::2] - extent[::2]
    ndims = len(dextent)

    if isinstance(levels, int):
        minlvl = maxlvl = levels
        levels = np.array([levels] * ndims, dtype="int32")
    else:
        levels = np.asarray(levels, dtype="int32")
        minlvl = levels.min()
        maxlvl = levels.max()
        if minlvl != maxlvl and (minlvl != 0 or {minlvl, maxlvl} != set(levels)):
            raise ValueError("all levels must have the same value or zero.")
    dims_zero = levels == 0
    dims_nonzero = ~dims_zero
    ndims_nonzero = dims_nonzero.sum()

    npoints = (cells**ndims_nonzero - 1) * maxlvl + 1
    left = np.empty((npoints, ndims), dtype="float64")
    right = np.empty((npoints, ndims), dtype="float64")
    level = np.empty(npoints, dtype="int32")

    # fill zero dims
    left[:, dims_zero] = extent[::2][dims_zero]
    right[:, dims_zero] = extent[1::2][dims_zero]

    # fill non-zero dims
    dcell = 1.0 / cells
    left_slice = tuple(
        slice(extent[2 * n], extent[2 * n + 1], extent[2 * n + 1])
        if dims_zero[n]
        else slice(0.0, 1.0, dcell)
        for n in range(ndims)
    )
    right_slice = tuple(
        slice(extent[2 * n + 1], extent[2 * n], -extent[2 * n + 1])
        if dims_zero[n]
        else slice(dcell, 1.0 + dcell, dcell)
        for n in range(ndims)
    )
    left_norm_grid = np.reshape(np.mgrid[left_slice].T.flat[ndims:], (-1, ndims))
    lng_zero = left_norm_grid[:, dims_zero]
    lng_nonzero = left_norm_grid[:, dims_nonzero]

    right_norm_grid = np.reshape(np.mgrid[right_slice].T.flat[ndims:], (-1, ndims))
    rng_zero = right_norm_grid[:, dims_zero]
    rng_nonzero = right_norm_grid[:, dims_nonzero]

    level[0] = maxlvl
    left[0, :] = extent[::2]
    right[0, dims_zero] = extent[1::2][dims_zero]
    right[0, dims_nonzero] = (dcell**maxlvl) * dextent[dims_nonzero] + extent[::2][
        dims_nonzero
    ]
    for i, lvl in enumerate(range(maxlvl, 0, -1)):
        start = (cells**ndims_nonzero - 1) * i + 1
        stop = (cells**ndims_nonzero - 1) * (i + 1) + 1
        dsize = dcell ** (lvl - 1) * dextent[dims_nonzero]
        level[start:stop] = lvl
        left[start:stop, dims_zero] = lng_zero
        left[start:stop, dims_nonzero] = lng_nonzero * dsize + extent[::2][dims_nonzero]
        right[start:stop, dims_zero] = rng_zero
        right[start:stop, dims_nonzero] = (
            rng_nonzero * dsize + extent[::2][dims_nonzero]
        )

    return left, right, level


def _check_field_unit_args_helper(args: dict, default_args: dict):
    values = list(args.values())
    keys = list(args.keys())
    if all(v is None for v in values):
        for key in keys:
            args[key] = default_args[key]
    elif None in values:
        raise ValueError(
            "Error in creating a fake dataset:"
            f" either all or none of the following arguments need to specified: {keys}."
        )
    elif any(len(v) != len(values[0]) for v in values):
        raise ValueError(
            "Error in creating a fake dataset:"
            f" all the following arguments must have the same length: {keys}."
        )
    return list(args.values())


_fake_random_ds_default_fields = ("density", "velocity_x", "velocity_y", "velocity_z")
_fake_random_ds_default_units = ("g/cm**3", "cm/s", "cm/s", "cm/s")
_fake_random_ds_default_negative = (False, False, False, False)


def fake_random_ds(
    ndims,
    peak_value=1.0,
    fields=None,
    units=None,
    particle_fields=None,
    particle_field_units=None,
    negative=False,
    nprocs=1,
    particles=0,
    length_unit=1.0,
    unit_system="cgs",
    bbox=None,
    default_species_fields=None,
):
    from yt.loaders import load_uniform_grid

    prng = RandomState(0x4D3D3D3)
    if not is_sequence(ndims):
        ndims = [ndims, ndims, ndims]
    else:
        assert len(ndims) == 3
    if not is_sequence(negative):
        if fields:
            negative = [negative for f in fields]
        else:
            negative = None

    fields, units, negative = _check_field_unit_args_helper(
        {
            "fields": fields,
            "units": units,
            "negative": negative,
        },
        {
            "fields": _fake_random_ds_default_fields,
            "units": _fake_random_ds_default_units,
            "negative": _fake_random_ds_default_negative,
        },
    )

    offsets = []
    for n in negative:
        if n:
            offsets.append(0.5)
        else:
            offsets.append(0.0)
    data = {}
    for field, offset, u in zip(fields, offsets, units):
        v = (prng.random_sample(ndims) - offset) * peak_value
        if field[0] == "all":
            v = v.ravel()
        data[field] = (v, u)
    if particles:
        if particle_fields is not None:
            for field, unit in zip(particle_fields, particle_field_units):
                if field in ("particle_position", "particle_velocity"):
                    data["io", field] = (prng.random_sample((int(particles), 3)), unit)
                else:
                    data["io", field] = (prng.random_sample(size=int(particles)), unit)
        else:
            for f in (f"particle_position_{ax}" for ax in "xyz"):
                data["io", f] = (prng.random_sample(size=particles), "code_length")
            for f in (f"particle_velocity_{ax}" for ax in "xyz"):
                data["io", f] = (prng.random_sample(size=particles) - 0.5, "cm/s")
            data["io", "particle_mass"] = (prng.random_sample(particles), "g")
    ug = load_uniform_grid(
        data,
        ndims,
        length_unit=length_unit,
        nprocs=nprocs,
        unit_system=unit_system,
        bbox=bbox,
        default_species_fields=default_species_fields,
    )
    return ug


_geom_transforms = {
    # These are the bounds we want.  Cartesian we just assume goes 0 .. 1.
    "cartesian": ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
    "spherical": ((0.0, 0.0, 0.0), (1.0, np.pi, 2 * np.pi)),
    "cylindrical": ((0.0, 0.0, 0.0), (1.0, 1.0, 2.0 * np.pi)),  # rzt
    "polar": ((0.0, 0.0, 0.0), (1.0, 2.0 * np.pi, 1.0)),  # rtz
    "geographic": ((-90.0, -180.0, 0.0), (90.0, 180.0, 1000.0)),  # latlonalt
    "internal_geographic": ((-90.0, -180.0, 0.0), (90.0, 180.0, 1000.0)),  # latlondep
}


_fake_amr_ds_default_fields = ("Density",)
_fake_amr_ds_default_units = ("g/cm**3",)


def fake_amr_ds(
    fields=None, units=None, geometry="cartesian", particles=0, length_unit=None
):
    from yt.loaders import load_amr_grids

    fields, units = _check_field_unit_args_helper(
        {
            "fields": fields,
            "units": units,
        },
        {
            "fields": _fake_amr_ds_default_fields,
            "units": _fake_amr_ds_default_units,
        },
    )

    prng = RandomState(0x4D3D3D3)
    LE, RE = _geom_transforms[geometry]
    LE = np.array(LE)
    RE = np.array(RE)
    data = []
    for gspec in _amr_grid_index:
        level, left_edge, right_edge, dims = gspec
        left_edge = left_edge * (RE - LE) + LE
        right_edge = right_edge * (RE - LE) + LE
        gdata = dict(
            level=level, left_edge=left_edge, right_edge=right_edge, dimensions=dims
        )
        for f, u in zip(fields, units):
            gdata[f] = (prng.random_sample(dims), u)
        if particles:
            for i, f in enumerate(f"particle_position_{ax}" for ax in "xyz"):
                pdata = prng.random_sample(particles)
                pdata /= right_edge[i] - left_edge[i]
                pdata += left_edge[i]
                gdata["io", f] = (pdata, "code_length")
            for f in (f"particle_velocity_{ax}" for ax in "xyz"):
                gdata["io", f] = (prng.random_sample(particles) - 0.5, "cm/s")
            gdata["io", "particle_mass"] = (prng.random_sample(particles), "g")
        data.append(gdata)
    bbox = np.array([LE, RE]).T
    return load_amr_grids(
        data, [32, 32, 32], geometry=geometry, bbox=bbox, length_unit=length_unit
    )


_fake_particle_ds_default_fields = (
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_mass",
    "particle_velocity_x",
    "particle_velocity_y",
    "particle_velocity_z",
)
_fake_particle_ds_default_units = ("cm", "cm", "cm", "g", "cm/s", "cm/s", "cm/s")
_fake_particle_ds_default_negative = (False, False, False, False, True, True, True)


def fake_particle_ds(
    fields=None,
    units=None,
    negative=None,
    npart=16**3,
    length_unit=1.0,
    data=None,
):
    from yt.loaders import load_particles

    prng = RandomState(0x4D3D3D3)
    if negative is not None and not is_sequence(negative):
        negative = [negative for f in fields]

    fields, units, negative = _check_field_unit_args_helper(
        {
            "fields": fields,
            "units": units,
            "negative": negative,
        },
        {
            "fields": _fake_particle_ds_default_fields,
            "units": _fake_particle_ds_default_units,
            "negative": _fake_particle_ds_default_negative,
        },
    )

    offsets = []
    for n in negative:
        if n:
            offsets.append(0.5)
        else:
            offsets.append(0.0)
    data = data if data else {}
    for field, offset, u in zip(fields, offsets, units):
        if field in data:
            v = data[field]
            continue
        if "position" in field:
            v = prng.normal(loc=0.5, scale=0.25, size=npart)
            np.clip(v, 0.0, 1.0, v)
        v = prng.random_sample(npart) - offset
        data[field] = (v, u)
    bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    ds = load_particles(data, 1.0, bbox=bbox)
    return ds


def fake_tetrahedral_ds():
    from yt.frontends.stream.sample_data.tetrahedral_mesh import (
        _connectivity,
        _coordinates,
    )
    from yt.loaders import load_unstructured_mesh

    prng = RandomState(0x4D3D3D3)

    # the distance from the origin
    node_data = {}
    dist = np.sum(_coordinates**2, 1)
    node_data[("connect1", "test")] = dist[_connectivity]

    # each element gets a random number
    elem_data = {}
    elem_data[("connect1", "elem")] = prng.rand(_connectivity.shape[0])

    ds = load_unstructured_mesh(
        _connectivity, _coordinates, node_data=node_data, elem_data=elem_data
    )
    return ds


def fake_hexahedral_ds(fields=None):
    from yt.frontends.stream.sample_data.hexahedral_mesh import (
        _connectivity,
        _coordinates,
    )
    from yt.loaders import load_unstructured_mesh

    prng = RandomState(0x4D3D3D3)
    # the distance from the origin
    node_data = {}
    dist = np.sum(_coordinates**2, 1)
    node_data[("connect1", "test")] = dist[_connectivity - 1]

    for field in always_iterable(fields):
        node_data[("connect1", field)] = dist[_connectivity - 1]

    # each element gets a random number
    elem_data = {}
    elem_data[("connect1", "elem")] = prng.rand(_connectivity.shape[0])

    ds = load_unstructured_mesh(
        _connectivity - 1, _coordinates, node_data=node_data, elem_data=elem_data
    )
    return ds


def small_fake_hexahedral_ds():
    from yt.loaders import load_unstructured_mesh

    _coordinates = np.array(
        [
            [-1.0, -1.0, -1.0],
            [0.0, -1.0, -1.0],
            [-0.0, 0.0, -1.0],
            [-1.0, -0.0, -1.0],
            [-1.0, -1.0, 0.0],
            [-0.0, -1.0, 0.0],
            [-0.0, 0.0, -0.0],
            [-1.0, 0.0, -0.0],
        ]
    )
    _connectivity = np.array([[1, 2, 3, 4, 5, 6, 7, 8]])

    # the distance from the origin
    node_data = {}
    dist = np.sum(_coordinates**2, 1)
    node_data[("connect1", "test")] = dist[_connectivity - 1]

    ds = load_unstructured_mesh(_connectivity - 1, _coordinates, node_data=node_data)
    return ds


def fake_stretched_ds(N=16):
    from yt.loaders import load_uniform_grid

    np.random.seed(0x4D3D3D3)
    data = {"density": np.random.random((N, N, N))}

    cell_widths = []
    for _ in range(3):
        cw = np.random.random(N)
        cw /= cw.sum()
        cell_widths.append(cw)
    return load_uniform_grid(
        data,
        [N, N, N],
        bbox=np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]),
        cell_widths=cell_widths,
    )


def fake_vr_orientation_test_ds(N=96, scale=1):
    """
    create a toy dataset that puts a sphere at (0,0,0), a single cube
    on +x, two cubes on +y, and three cubes on +z in a domain from
    [-1*scale,1*scale]**3.  The lower planes
    (x = -1*scale, y = -1*scale, z = -1*scale) are also given non-zero
    values.

    This dataset allows you to easily explore orientations and
    handiness in VR and other renderings

    Parameters
    ----------

    N : integer
       The number of cells along each direction

    scale : float
       A spatial scale, the domain boundaries will be multiplied by scale to
       test datasets that have spatial different scales (e.g. data in CGS units)

    """
    from yt.loaders import load_uniform_grid

    xmin = ymin = zmin = -1.0 * scale
    xmax = ymax = zmax = 1.0 * scale

    dcoord = (xmax - xmin) / N

    arr = np.zeros((N, N, N), dtype=np.float64)
    arr[:, :, :] = 1.0e-4

    bbox = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    # coordinates -- in the notation data[i, j, k]
    x = (np.arange(N) + 0.5) * dcoord + xmin
    y = (np.arange(N) + 0.5) * dcoord + ymin
    z = (np.arange(N) + 0.5) * dcoord + zmin

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    # sphere at the origin
    c = np.array([0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)])
    r = np.sqrt((x3d - c[0]) ** 2 + (y3d - c[1]) ** 2 + (z3d - c[2]) ** 2)
    arr[r < 0.05] = 1.0

    arr[abs(x3d - xmin) < 2 * dcoord] = 0.3
    arr[abs(y3d - ymin) < 2 * dcoord] = 0.3
    arr[abs(z3d - zmin) < 2 * dcoord] = 0.3

    # single cube on +x
    xc = 0.75 * scale
    dx = 0.05 * scale
    idx = np.logical_and(
        np.logical_and(x3d > xc - dx, x3d < xc + dx),
        np.logical_and(
            np.logical_and(y3d > -dx, y3d < dx), np.logical_and(z3d > -dx, z3d < dx)
        ),
    )
    arr[idx] = 1.0

    # two cubes on +y
    dy = 0.05 * scale
    for yc in [0.65 * scale, 0.85 * scale]:
        idx = np.logical_and(
            np.logical_and(y3d > yc - dy, y3d < yc + dy),
            np.logical_and(
                np.logical_and(x3d > -dy, x3d < dy), np.logical_and(z3d > -dy, z3d < dy)
            ),
        )
        arr[idx] = 0.8

    # three cubes on +z
    dz = 0.05 * scale
    for zc in [0.5 * scale, 0.7 * scale, 0.9 * scale]:
        idx = np.logical_and(
            np.logical_and(z3d > zc - dz, z3d < zc + dz),
            np.logical_and(
                np.logical_and(x3d > -dz, x3d < dz), np.logical_and(y3d > -dz, y3d < dz)
            ),
        )
        arr[idx] = 0.6

    data = dict(density=(arr, "g/cm**3"))
    ds = load_uniform_grid(data, arr.shape, bbox=bbox)
    return ds


def fake_sph_orientation_ds():
    """Returns an in-memory SPH dataset useful for testing

    This dataset should have one particle at the origin, one more particle
    along the x axis, two along y, and three along z. All particles will
    have non-overlapping smoothing regions with a radius of 0.25, masses of 1,
    and densities of 1, and zero velocity.
    """
    from yt import load_particles

    npart = 7

    # one particle at the origin, one particle along x-axis, two along y,
    # three along z
    data = {
        "particle_position_x": (np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), "cm"),
        "particle_position_y": (np.array([0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0]), "cm"),
        "particle_position_z": (np.array([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0]), "cm"),
        "particle_mass": (np.ones(npart), "g"),
        "particle_velocity_x": (np.zeros(npart), "cm/s"),
        "particle_velocity_y": (np.zeros(npart), "cm/s"),
        "particle_velocity_z": (np.zeros(npart), "cm/s"),
        "smoothing_length": (0.25 * np.ones(npart), "cm"),
        "density": (np.ones(npart), "g/cm**3"),
        "temperature": (np.ones(npart), "K"),
    }

    bbox = np.array([[-4, 4], [-4, 4], [-4, 4]])

    return load_particles(data=data, length_unit=1.0, bbox=bbox)


def fake_sph_grid_ds(hsml_factor=1.0):
    """Returns an in-memory SPH dataset useful for testing

    This dataset should have 27 particles with the particles arranged uniformly
    on a 3D grid. The bottom left corner is (0.5,0.5,0.5) and the top right
    corner is (2.5,2.5,2.5). All particles will have non-overlapping smoothing
    regions with a radius of 0.05, masses of 1, and densities of 1, and zero
    velocity.
    """
    from yt import load_particles

    npart = 27

    x = np.empty(npart)
    y = np.empty(npart)
    z = np.empty(npart)

    tot = 0
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                x[tot] = i + 0.5
                y[tot] = j + 0.5
                z[tot] = k + 0.5
                tot += 1

    data = {
        "particle_position_x": (x, "cm"),
        "particle_position_y": (y, "cm"),
        "particle_position_z": (z, "cm"),
        "particle_mass": (np.ones(npart), "g"),
        "particle_velocity_x": (np.zeros(npart), "cm/s"),
        "particle_velocity_y": (np.zeros(npart), "cm/s"),
        "particle_velocity_z": (np.zeros(npart), "cm/s"),
        "smoothing_length": (0.05 * np.ones(npart) * hsml_factor, "cm"),
        "density": (np.ones(npart), "g/cm**3"),
        "temperature": (np.ones(npart), "K"),
    }

    bbox = np.array([[0, 3], [0, 3], [0, 3]])

    return load_particles(data=data, length_unit=1.0, bbox=bbox)


def construct_octree_mask(prng=RandomState(0x1D3D3D3), refined=None):  # noqa B008
    # Implementation taken from url:
    # http://docs.hyperion-rt.org/en/stable/advanced/indepth_oct.html

    if refined in (None, True):
        refined = [True]
    if not refined:
        refined = [False]
        return refined

    # Loop over subcells
    for _ in range(8):
        # Insert criterion for whether cell should be sub-divided. Here we
        # just use a random number to demonstrate.
        divide = prng.random_sample() < 0.12

        # Append boolean to overall list
        refined.append(divide)

        # If the cell is sub-divided, recursively divide it further
        if divide:
            construct_octree_mask(prng, refined)
    return refined


def fake_octree_ds(
    prng=RandomState(0x4D3D3D3),  # noqa B008
    refined=None,
    quantities=None,
    bbox=None,
    sim_time=0.0,
    length_unit=None,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    num_zones=2,
    partial_coverage=1,
    unit_system="cgs",
):
    from yt.loaders import load_octree

    octree_mask = np.asarray(
        construct_octree_mask(prng=prng, refined=refined), dtype=np.uint8
    )
    particles = np.sum(np.invert(octree_mask))

    if quantities is None:
        quantities = {}
        quantities[("gas", "density")] = prng.random_sample((particles, 1))
        quantities[("gas", "velocity_x")] = prng.random_sample((particles, 1))
        quantities[("gas", "velocity_y")] = prng.random_sample((particles, 1))
        quantities[("gas", "velocity_z")] = prng.random_sample((particles, 1))

    ds = load_octree(
        octree_mask=octree_mask,
        data=quantities,
        bbox=bbox,
        sim_time=sim_time,
        length_unit=length_unit,
        mass_unit=mass_unit,
        time_unit=time_unit,
        velocity_unit=velocity_unit,
        magnetic_unit=magnetic_unit,
        periodicity=periodicity,
        partial_coverage=partial_coverage,
        num_zones=num_zones,
        unit_system=unit_system,
    )
    return ds


def add_noise_fields(ds):
    """Add 4 classes of noise fields to a dataset"""
    prng = RandomState(0x4D3D3D3)

    def _binary_noise(field, data):
        """random binary data"""
        return prng.randint(low=0, high=2, size=data.size).astype("float64")

    def _positive_noise(field, data):
        """random strictly positive data"""
        return prng.random_sample(data.size) + 1e-16

    def _negative_noise(field, data):
        """random negative data"""
        return -prng.random_sample(data.size)

    def _even_noise(field, data):
        """random data with mixed signs"""
        return 2 * prng.random_sample(data.size) - 1

    ds.add_field(("gas", "noise0"), _binary_noise, sampling_type="cell")
    ds.add_field(("gas", "noise1"), _positive_noise, sampling_type="cell")
    ds.add_field(("gas", "noise2"), _negative_noise, sampling_type="cell")
    ds.add_field(("gas", "noise3"), _even_noise, sampling_type="cell")


def expand_keywords(keywords, full=False):
    """
    expand_keywords is a means for testing all possible keyword
    arguments in the nosetests.  Simply pass it a dictionary of all the
    keyword arguments and all of the values for these arguments that you
    want to test.

    It will return a list of kwargs dicts containing combinations of
    the various kwarg values you passed it.  These can then be passed
    to the appropriate function in nosetests.

    If full=True, then every possible combination of keywords is produced,
    otherwise, every keyword option is included at least once in the output
    list.  Be careful, by using full=True, you may be in for an exponentially
    larger number of tests!

    Parameters
    ----------

    keywords : dict
        a dictionary where the keys are the keywords for the function,
        and the values of each key are the possible values that this key
        can take in the function

    full : bool
        if set to True, every possible combination of given keywords is
        returned

    Returns
    -------

    array of dicts
        An array of dictionaries to be individually passed to the appropriate
        function matching these kwargs.

    Examples
    --------

    >>> keywords = {}
    >>> keywords["dpi"] = (50, 100, 200)
    >>> keywords["cmap"] = ("cmyt.arbre", "cmyt.kelp")
    >>> list_of_kwargs = expand_keywords(keywords)
    >>> print(list_of_kwargs)

    array([{'cmap': 'cmyt.arbre', 'dpi': 50},
           {'cmap': 'cmyt.kelp', 'dpi': 100},
           {'cmap': 'cmyt.arbre', 'dpi': 200}], dtype=object)

    >>> list_of_kwargs = expand_keywords(keywords, full=True)
    >>> print(list_of_kwargs)

    array([{'cmap': 'cmyt.arbre', 'dpi': 50},
           {'cmap': 'cmyt.arbre', 'dpi': 100},
           {'cmap': 'cmyt.arbre', 'dpi': 200},
           {'cmap': 'cmyt.kelp', 'dpi': 50},
           {'cmap': 'cmyt.kelp', 'dpi': 100},
           {'cmap': 'cmyt.kelp', 'dpi': 200}], dtype=object)

    >>> for kwargs in list_of_kwargs:
    ...     write_projection(*args, **kwargs)
    """

    # if we want every possible combination of keywords, use iter magic
    if full:
        keys = sorted(keywords)
        list_of_kwarg_dicts = np.array(
            [
                dict(zip(keys, prod))
                for prod in it.product(*(keywords[key] for key in keys))
            ]
        )

    # if we just want to probe each keyword, but not necessarily every
    # combination
    else:
        # Determine the maximum number of values any of the keywords has
        num_lists = 0
        for val in keywords.values():
            if isinstance(val, str):
                num_lists = max(1.0, num_lists)
            else:
                num_lists = max(len(val), num_lists)

        # Construct array of kwargs dicts, each element of the list is a different
        # **kwargs dict.  each kwargs dict gives a different combination of
        # the possible values of the kwargs

        # initialize array
        list_of_kwarg_dicts = np.array([dict() for x in range(num_lists)])

        # fill in array
        for i in np.arange(num_lists):
            list_of_kwarg_dicts[i] = {}
            for key in keywords.keys():
                # if it's a string, use it (there's only one)
                if isinstance(keywords[key], str):
                    list_of_kwarg_dicts[i][key] = keywords[key]
                # if there are more options, use the i'th val
                elif i < len(keywords[key]):
                    list_of_kwarg_dicts[i][key] = keywords[key][i]
                # if there are not more options, use the 0'th val
                else:
                    list_of_kwarg_dicts[i][key] = keywords[key][0]

    return list_of_kwarg_dicts


def skip(reason: str):
    # a drop-in replacement for pytest.mark.skip decorator with nose-compatibility
    def dec(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            raise SkipTest(reason)

        return wrapper

    return dec


def skipif(condition: bool, reason: str):
    # a drop-in replacement for pytest.mark.skipif decorator with nose-compatibility
    def dec(func):
        if condition:
            return skip(reason)(func)
        else:
            return func

    return dec


def requires_module(module):
    """
    Decorator that takes a module name as an argument and tries to import it.
    If the module imports without issue, the function is returned, but if not,
    a null function is returned. This is so tests that depend on certain modules
    being imported will not fail if the module is not installed on the testing
    platform.
    """
    return skipif(find_spec(module) is None, reason=f"Missing required module {module}")


def requires_module_pytest(*module_names):
    """
    This is a replacement for yt.testing.requires_module that's
    compatible with pytest, and accepts an arbitrary number of requirements to
    avoid stacking decorators

    Important: this is meant to decorate test functions only, it won't work as a
    decorator to fixture functions.
    It's meant to be imported as
    >>> from yt.testing import requires_module_pytest as requires_module

    So that it can be later renamed to `requires_module`.
    """

    def deco(func):
        missing = [name for name in module_names if find_spec(name) is None]

        # note that order between these two decorators matters
        @pytest.mark.skipif(
            missing,
            reason=f"missing requirement(s): {', '.join(missing)}",
        )
        @wraps(func)
        def inner_func(*args, **kwargs):
            return func(*args, **kwargs)

        return inner_func

    return deco


def requires_file(req_file):
    condition = (
        not os.path.exists(req_file)
        and not os.path.exists(os.path.join(ytcfg.get("yt", "test_data_dir"), req_file))
        and not ytcfg.get("yt", "internals", "strict_requires")
    )
    return skipif(condition, reason=f"Missing required file {req_file}")


def disable_dataset_cache(func):
    @wraps(func)
    def newfunc(*args, **kwargs):
        restore_cfg_state = False
        if not ytcfg.get("yt", "skip_dataset_cache"):
            ytcfg["yt", "skip_dataset_cache"] = True
            restore_cfg_state = True
        rv = func(*args, **kwargs)
        if restore_cfg_state:
            ytcfg["yt", "skip_dataset_cache"] = False
        return rv

    return newfunc


@disable_dataset_cache
def units_override_check(fn):
    units_list = ["length", "time", "mass", "velocity", "magnetic", "temperature"]
    ds1 = load(fn)
    units_override = {}
    attrs1 = []
    attrs2 = []
    for u in units_list:
        unit_attr = getattr(ds1, f"{u}_unit", None)
        if unit_attr is not None:
            attrs1.append(unit_attr)
            units_override[f"{u}_unit"] = (unit_attr.v, unit_attr.units)
    del ds1
    ds2 = load(fn, units_override=units_override)
    assert len(ds2.units_override) > 0
    for u in units_list:
        unit_attr = getattr(ds2, f"{u}_unit", None)
        if unit_attr is not None:
            attrs2.append(unit_attr)
    assert_equal(attrs1, attrs2)


# This is an export of the 40 grids in IsolatedGalaxy that are of level 4 or
# lower.  It's just designed to give a sample AMR index to deal with.
_amr_grid_index = [
    [0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [32, 32, 32]],
    [1, [0.25, 0.21875, 0.25], [0.5, 0.5, 0.5], [16, 18, 16]],
    [1, [0.5, 0.21875, 0.25], [0.75, 0.5, 0.5], [16, 18, 16]],
    [1, [0.21875, 0.5, 0.25], [0.5, 0.75, 0.5], [18, 16, 16]],
    [1, [0.5, 0.5, 0.25], [0.75, 0.75, 0.5], [16, 16, 16]],
    [1, [0.25, 0.25, 0.5], [0.5, 0.5, 0.75], [16, 16, 16]],
    [1, [0.5, 0.25, 0.5], [0.75, 0.5, 0.75], [16, 16, 16]],
    [1, [0.25, 0.5, 0.5], [0.5, 0.75, 0.75], [16, 16, 16]],
    [1, [0.5, 0.5, 0.5], [0.75, 0.75, 0.75], [16, 16, 16]],
    [2, [0.5, 0.5, 0.5], [0.71875, 0.71875, 0.71875], [28, 28, 28]],
    [3, [0.5, 0.5, 0.5], [0.6640625, 0.65625, 0.6796875], [42, 40, 46]],
    [4, [0.5, 0.5, 0.5], [0.59765625, 0.6015625, 0.6015625], [50, 52, 52]],
    [2, [0.28125, 0.5, 0.5], [0.5, 0.734375, 0.71875], [28, 30, 28]],
    [3, [0.3359375, 0.5, 0.5], [0.5, 0.671875, 0.6640625], [42, 44, 42]],
    [4, [0.40625, 0.5, 0.5], [0.5, 0.59765625, 0.59765625], [48, 50, 50]],
    [2, [0.5, 0.28125, 0.5], [0.71875, 0.5, 0.71875], [28, 28, 28]],
    [3, [0.5, 0.3359375, 0.5], [0.671875, 0.5, 0.6640625], [44, 42, 42]],
    [4, [0.5, 0.40625, 0.5], [0.6015625, 0.5, 0.59765625], [52, 48, 50]],
    [2, [0.28125, 0.28125, 0.5], [0.5, 0.5, 0.71875], [28, 28, 28]],
    [3, [0.3359375, 0.3359375, 0.5], [0.5, 0.5, 0.671875], [42, 42, 44]],
    [
        4,
        [0.46484375, 0.37890625, 0.50390625],
        [0.4765625, 0.390625, 0.515625],
        [6, 6, 6],
    ],
    [4, [0.40625, 0.40625, 0.5], [0.5, 0.5, 0.59765625], [48, 48, 50]],
    [2, [0.5, 0.5, 0.28125], [0.71875, 0.71875, 0.5], [28, 28, 28]],
    [3, [0.5, 0.5, 0.3359375], [0.6796875, 0.6953125, 0.5], [46, 50, 42]],
    [4, [0.5, 0.5, 0.40234375], [0.59375, 0.6015625, 0.5], [48, 52, 50]],
    [2, [0.265625, 0.5, 0.28125], [0.5, 0.71875, 0.5], [30, 28, 28]],
    [3, [0.3359375, 0.5, 0.328125], [0.5, 0.65625, 0.5], [42, 40, 44]],
    [4, [0.40234375, 0.5, 0.40625], [0.5, 0.60546875, 0.5], [50, 54, 48]],
    [2, [0.5, 0.265625, 0.28125], [0.71875, 0.5, 0.5], [28, 30, 28]],
    [3, [0.5, 0.3203125, 0.328125], [0.6640625, 0.5, 0.5], [42, 46, 44]],
    [4, [0.5, 0.3984375, 0.40625], [0.546875, 0.5, 0.5], [24, 52, 48]],
    [4, [0.546875, 0.41796875, 0.4453125], [0.5625, 0.4375, 0.5], [8, 10, 28]],
    [4, [0.546875, 0.453125, 0.41796875], [0.5546875, 0.48046875, 0.4375], [4, 14, 10]],
    [4, [0.546875, 0.4375, 0.4375], [0.609375, 0.5, 0.5], [32, 32, 32]],
    [4, [0.546875, 0.4921875, 0.41796875], [0.56640625, 0.5, 0.4375], [10, 4, 10]],
    [
        4,
        [0.546875, 0.48046875, 0.41796875],
        [0.5703125, 0.4921875, 0.4375],
        [12, 6, 10],
    ],
    [4, [0.55859375, 0.46875, 0.43359375], [0.5703125, 0.48046875, 0.4375], [6, 6, 2]],
    [2, [0.265625, 0.28125, 0.28125], [0.5, 0.5, 0.5], [30, 28, 28]],
    [3, [0.328125, 0.3359375, 0.328125], [0.5, 0.5, 0.5], [44, 42, 44]],
    [4, [0.4140625, 0.40625, 0.40625], [0.5, 0.5, 0.5], [44, 48, 48]],
]


def check_results(func):
    r"""This is a decorator for a function to verify that the (numpy ndarray)
    result of a function is what it should be.

    This function is designed to be used for very light answer testing.
    Essentially, it wraps around a larger function that returns a numpy array,
    and that has results that should not change.  It is not necessarily used
    inside the testing scripts themselves, but inside testing scripts written
    by developers during the testing of pull requests and new functionality.
    If a hash is specified, it "wins" and the others are ignored.  Otherwise,
    tolerance is 1e-8 (just above single precision.)

    The correct results will be stored if the command line contains
    --answer-reference , and otherwise it will compare against the results on
    disk.  The filename will be func_results_ref_FUNCNAME.cpkl where FUNCNAME
    is the name of the function being tested.

    If you would like more control over the name of the pickle file the results
    are stored in, you can pass the result_basename keyword argument to the
    function you are testing.  The check_results decorator will use the value
    of the keyword to construct the filename of the results data file.  If
    result_basename is not specified, the name of the testing function is used.

    This will raise an exception if the results are not correct.

    Examples
    --------

    >>> @check_results
    ... def my_func(ds):
    ...     return ds.domain_width

    >>> my_func(ds)

    >>> @check_results
    ... def field_checker(dd, field_name):
    ...     return dd[field_name]

    >>> field_checker(ds.all_data(), "density", result_basename="density")

    """

    def compute_results(func):
        @wraps(func)
        def _func(*args, **kwargs):
            name = kwargs.pop("result_basename", func.__name__)
            rv = func(*args, **kwargs)
            if hasattr(rv, "convert_to_base"):
                rv.convert_to_base()
                _rv = rv.ndarray_view()
            else:
                _rv = rv
            mi = _rv.min()
            ma = _rv.max()
            st = _rv.std(dtype="float64")
            su = _rv.sum(dtype="float64")
            si = _rv.size
            ha = hashlib.md5(_rv.tobytes()).hexdigest()
            fn = f"func_results_ref_{name}.cpkl"
            with open(fn, "wb") as f:
                pickle.dump((mi, ma, st, su, si, ha), f)
            return rv

        return _func

    from yt.mods import unparsed_args

    if "--answer-reference" in unparsed_args:
        return compute_results(func)

    def compare_results(func):
        @wraps(func)
        def _func(*args, **kwargs):
            name = kwargs.pop("result_basename", func.__name__)
            rv = func(*args, **kwargs)
            if hasattr(rv, "convert_to_base"):
                rv.convert_to_base()
                _rv = rv.ndarray_view()
            else:
                _rv = rv
            vals = (
                _rv.min(),
                _rv.max(),
                _rv.std(dtype="float64"),
                _rv.sum(dtype="float64"),
                _rv.size,
                hashlib.md5(_rv.tobytes()).hexdigest(),
            )
            fn = f"func_results_ref_{name}.cpkl"
            if not os.path.exists(fn):
                print("Answers need to be created with --answer-reference .")
                return False
            with open(fn, "rb") as f:
                ref = pickle.load(f)
            print(f"Sizes: {vals[4] == ref[4]} ({vals[4]}, {ref[4]})")
            assert_allclose(vals[0], ref[0], 1e-8, err_msg="min")
            assert_allclose(vals[1], ref[1], 1e-8, err_msg="max")
            assert_allclose(vals[2], ref[2], 1e-8, err_msg="std")
            assert_allclose(vals[3], ref[3], 1e-8, err_msg="sum")
            assert_equal(vals[4], ref[4])
            print("Hashes equal: %s" % (vals[-1] == ref[-1]))
            return rv

        return _func

    return compare_results(func)


def periodicity_cases(ds):
    # This is a generator that yields things near the corners.  It's good for
    # getting different places to check periodicity.
    yield (ds.domain_left_edge + ds.domain_right_edge) / 2.0
    dx = ds.domain_width / ds.domain_dimensions
    # We start one dx in, and only go to one in as well.
    for i in (1, ds.domain_dimensions[0] - 2):
        for j in (1, ds.domain_dimensions[1] - 2):
            for k in (1, ds.domain_dimensions[2] - 2):
                center = dx * np.array([i, j, k]) + ds.domain_left_edge
                yield center


def run_nose(
    verbose=False,
    run_answer_tests=False,
    answer_big_data=False,
    call_pdb=False,
    module=None,
):
    issue_deprecation_warning(
        "yt.run_nose (aka yt.testing.run_nose) is deprecated. "
        "Please do not rely on this function as it will be removed "
        "in the process of migrating yt tests from nose to pytest.",
        since="4.1.0",
    )

    from yt.utilities.logger import ytLogger as mylog
    from yt.utilities.on_demand_imports import _nose

    orig_level = mylog.getEffectiveLevel()
    mylog.setLevel(50)
    nose_argv = sys.argv
    nose_argv += ["--exclude=answer_testing", "--detailed-errors", "--exe"]
    if call_pdb:
        nose_argv += ["--pdb", "--pdb-failures"]
    if verbose:
        nose_argv.append("-v")
    if run_answer_tests:
        nose_argv.append("--with-answer-testing")
    if answer_big_data:
        nose_argv.append("--answer-big-data")
    if module:
        nose_argv.append(module)
    initial_dir = os.getcwd()
    yt_file = os.path.abspath(__file__)
    yt_dir = os.path.dirname(yt_file)
    if os.path.samefile(os.path.dirname(yt_dir), initial_dir):
        # Provide a nice error message to work around nose bug
        # see https://github.com/nose-devs/nose/issues/701
        raise RuntimeError(
            """
    The yt.run_nose function does not work correctly when invoked in
    the same directory as the installed yt package. Try starting
    a python session in a different directory before invoking yt.run_nose
    again. Alternatively, you can also run the "nosetests" executable in
    the current directory like so:

        $ nosetests
            """
        )
    os.chdir(yt_dir)
    try:
        _nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
        mylog.setLevel(orig_level)


def assert_allclose_units(actual, desired, rtol=1e-7, atol=0, **kwargs):
    """Raise an error if two objects are not equal up to desired tolerance

    This is a wrapper for :func:`numpy.testing.assert_allclose` that also
    verifies unit consistency

    Parameters
    ----------
    actual : array-like
        Array obtained (possibly with attached units)
    desired : array-like
        Array to compare with (possibly with attached units)
    rtol : float, optional
        Relative tolerance, defaults to 1e-7
    atol : float or quantity, optional
        Absolute tolerance. If units are attached, they must be consistent
        with the units of ``actual`` and ``desired``. If no units are attached,
        assumes the same units as ``desired``. Defaults to zero.

    Notes
    -----
    Also accepts additional keyword arguments accepted by
    :func:`numpy.testing.assert_allclose`, see the documentation of that
    function for details.

    """
    # Create a copy to ensure this function does not alter input arrays
    act = YTArray(actual)
    des = YTArray(desired)

    try:
        des = des.in_units(act.units)
    except UnitOperationError as e:
        raise AssertionError(
            "Units of actual (%s) and desired (%s) do not have "
            "equivalent dimensions" % (act.units, des.units)
        ) from e

    rt = YTArray(rtol)
    if not rt.units.is_dimensionless:
        raise AssertionError(f"Units of rtol ({rt.units}) are not dimensionless")

    if not isinstance(atol, YTArray):
        at = YTQuantity(atol, des.units)

    try:
        at = at.in_units(act.units)
    except UnitOperationError as e:
        raise AssertionError(
            "Units of atol (%s) and actual (%s) do not have "
            "equivalent dimensions" % (at.units, act.units)
        ) from e

    # units have been validated, so we strip units before calling numpy
    # to avoid spurious errors
    act = act.value
    des = des.value
    rt = rt.value
    at = at.value

    return assert_allclose(act, des, rt, at, **kwargs)


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    with open(fname, "rb") as fimg:
        data = fimg.read()
    image_type = ""

    # see http://www.w3.org/TR/PNG/#5PNG-file-signature
    if data.startswith(b"\211PNG\r\n\032\n"):
        image_type = ".png"
    # see http://www.mathguide.de/info/tools/media-types/image/jpeg
    elif data.startswith(b"\377\330"):
        image_type = ".jpeg"
    elif data.startswith(b"%!PS-Adobe"):
        data_str = data.decode("utf-8", "ignore")
        if "EPSF" in data_str[: data_str.index("\n")]:
            image_type = ".eps"
        else:
            image_type = ".ps"
    elif data.startswith(b"%PDF"):
        image_type = ".pdf"

    extension = os.path.splitext(fname)[1]

    assert (
        image_type == extension
    ), "Expected an image of type '{}' but '{}' is an image of type '{}'".format(
        extension,
        fname,
        image_type,
    )


def requires_backend(backend):
    """Decorator to check for a specified matplotlib backend.

    This decorator returns the decorated function if the specified `backend`
    is same as of `matplotlib.get_backend()`, otherwise returns null function.
    It could be used to execute function only when a particular `backend` of
    matplotlib is being used.

    Parameters
    ----------
    backend : String
        The value which is compared with the current matplotlib backend in use.

    """
    return skipif(
        backend.lower() != matplotlib.get_backend().lower(),
        reason=f"'{backend}' backend not in use",
    )


def requires_external_executable(*names):
    missing = [name for name in names if which(name) is None]
    return skipif(
        len(missing) > 0, reason=f"missing external executable(s): {', '.join(missing)}"
    )


class TempDirTest(unittest.TestCase):
    """
    A test class that runs in a temporary directory and
    removes it afterward.
    """

    def setUp(self):
        self.curdir = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)


class ParticleSelectionComparison:
    """
    This is a test helper class that takes a particle dataset, caches the
    particles it has on disk (manually reading them using lower-level IO
    routines) and then received a data object that it compares against manually
    running the data object's selection routines.  All supplied data objects
    must be created from the input dataset.
    """

    def __init__(self, ds):
        self.ds = ds
        # Construct an index so that we get all the data_files
        ds.index
        particles = {}
        # hsml is the smoothing length we use for radial selection
        hsml = {}
        for data_file in ds.index.data_files:
            for ptype, pos_arr in ds.index.io._yield_coordinates(data_file):
                particles.setdefault(ptype, []).append(pos_arr)
                if ptype in getattr(ds, "_sph_ptypes", ()):
                    hsml.setdefault(ptype, []).append(
                        ds.index.io._get_smoothing_length(
                            data_file, pos_arr.dtype, pos_arr.shape
                        )
                    )
        for ptype in particles:
            particles[ptype] = np.concatenate(particles[ptype])
            if ptype in hsml:
                hsml[ptype] = np.concatenate(hsml[ptype])
        self.particles = particles
        self.hsml = hsml

    def compare_dobj_selection(self, dobj):
        for ptype in sorted(self.particles):
            x, y, z = self.particles[ptype].T
            # Set our radii to zero for now, I guess?
            radii = self.hsml.get(ptype, 0.0)
            sel_index = dobj.selector.select_points(x, y, z, radii)
            if sel_index is None:
                sel_pos = np.empty((0, 3))
            else:
                sel_pos = self.particles[ptype][sel_index, :]

            obj_results = []
            for chunk in dobj.chunks([], "io"):
                obj_results.append(chunk[ptype, "particle_position"])
            if any(_.size > 0 for _ in obj_results):
                obj_results = np.concatenate(obj_results, axis=0)
            else:
                obj_results = np.empty((0, 3))
            # Sometimes we get unitary scaling or other floating point noise. 5
            # NULP should be OK.  This is mostly for stuff like Rockstar, where
            # the f32->f64 casting happens at different places depending on
            # which code path we use.
            assert_array_almost_equal_nulp(sel_pos, obj_results, 5)

    def run_defaults(self):
        """
        This runs lots of samples that touch different types of wraparounds.

        Specifically, it does:

            * sphere in center with radius 0.1 unitary
            * sphere in center with radius 0.2 unitary
            * sphere in each of the eight corners of the domain with radius 0.1 unitary
            * sphere in center with radius 0.5 unitary
            * box that covers 0.1 .. 0.9
            * box from 0.8 .. 0.85
            * box from 0.3..0.6, 0.2..0.8, 0.0..0.1
        """
        sp1 = self.ds.sphere("c", (0.1, "unitary"))
        self.compare_dobj_selection(sp1)

        sp2 = self.ds.sphere("c", (0.2, "unitary"))
        self.compare_dobj_selection(sp2)

        centers = [
            [0.04, 0.5, 0.5],
            [0.5, 0.04, 0.5],
            [0.5, 0.5, 0.04],
            [0.04, 0.04, 0.04],
            [0.96, 0.5, 0.5],
            [0.5, 0.96, 0.5],
            [0.5, 0.5, 0.96],
            [0.96, 0.96, 0.96],
        ]
        r = self.ds.quan(0.1, "unitary")
        for center in centers:
            c = self.ds.arr(center, "unitary") + self.ds.domain_left_edge.in_units(
                "unitary"
            )
            if not all(self.ds.periodicity):
                # filter out the periodic bits for non-periodic datasets
                if any(c - r < self.ds.domain_left_edge) or any(
                    c + r > self.ds.domain_right_edge
                ):
                    continue
            sp = self.ds.sphere(c, (0.1, "unitary"))
            self.compare_dobj_selection(sp)

        sp = self.ds.sphere("c", (0.5, "unitary"))
        self.compare_dobj_selection(sp)

        dd = self.ds.all_data()
        self.compare_dobj_selection(dd)

        # This is in raw numbers, so we can offset for the left edge
        LE = self.ds.domain_left_edge.in_units("unitary").d

        reg1 = self.ds.r[
            (0.1 + LE[0], "unitary") : (0.9 + LE[0], "unitary"),
            (0.1 + LE[1], "unitary") : (0.9 + LE[1], "unitary"),
            (0.1 + LE[2], "unitary") : (0.9 + LE[2], "unitary"),
        ]
        self.compare_dobj_selection(reg1)

        reg2 = self.ds.r[
            (0.8 + LE[0], "unitary") : (0.85 + LE[0], "unitary"),
            (0.8 + LE[1], "unitary") : (0.85 + LE[1], "unitary"),
            (0.8 + LE[2], "unitary") : (0.85 + LE[2], "unitary"),
        ]
        self.compare_dobj_selection(reg2)

        reg3 = self.ds.r[
            (0.3 + LE[0], "unitary") : (0.6 + LE[0], "unitary"),
            (0.2 + LE[1], "unitary") : (0.8 + LE[1], "unitary"),
            (0.0 + LE[2], "unitary") : (0.1 + LE[2], "unitary"),
        ]
        self.compare_dobj_selection(reg3)
