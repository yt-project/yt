"""
Utilities to aid testing.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import hashlib
import cPickle
import itertools as it
import numpy as np
import importlib
import os
from yt.funcs import *
from yt.config import ytcfg
from numpy.testing import assert_array_equal, assert_almost_equal, \
    assert_approx_equal, assert_array_almost_equal, assert_equal, \
    assert_array_less, assert_string_equal, assert_array_almost_equal_nulp,\
    assert_allclose, assert_raises
from yt.units.yt_array import uconcatenate
import yt.fields.api as field_api
from yt.convenience import load


def assert_rel_equal(a1, a2, decimals, err_msg='', verbose=True):
    # We have nan checks in here because occasionally we have fields that get
    # weighted without non-zero weights.  I'm looking at you, particle fields!
    if isinstance(a1, np.ndarray):
        assert(a1.size == a2.size)
        # Mask out NaNs
        assert((np.isnan(a1) == np.isnan(a2)).all())
        a1[np.isnan(a1)] = 1.0
        a2[np.isnan(a2)] = 1.0
        # Mask out 0
        ind1 = np.array(np.abs(a1) < np.finfo(a1.dtype).eps)
        ind2 = np.array(np.abs(a2) < np.finfo(a2.dtype).eps)
        assert((ind1 == ind2).all())
        a1[ind1] = 1.0
        a2[ind2] = 1.0
    elif np.any(np.isnan(a1)) and np.any(np.isnan(a2)):
        return True
    if a1 == a2 == 0.0:
        # NANS!
        a1 = a2 = 1.0
    return assert_almost_equal(np.array(a1)/np.array(a2), 1.0, decimals, err_msg=err_msg,
                               verbose=verbose)

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
    >>> l, r, lvl = amrspace([0.0, 2.0, 1.0, 2.0, 0.0, 3.14], levels=(3,3,0), cells=2)
    >>> print l
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
    extent = np.asarray(extent, dtype='f8')
    dextent = extent[1::2] - extent[::2]
    ndims = len(dextent)

    if isinstance(levels, int):
        minlvl = maxlvl = levels
        levels = np.array([levels]*ndims, dtype='int32')
    else:
        levels = np.asarray(levels, dtype='int32')
        minlvl = levels.min()
        maxlvl = levels.max()
        if minlvl != maxlvl and (minlvl != 0 or set([minlvl, maxlvl]) != set(levels)):
            raise ValueError("all levels must have the same value or zero.")
    dims_zero = (levels == 0)
    dims_nonzero = ~dims_zero
    ndims_nonzero = dims_nonzero.sum()

    npoints = (cells**ndims_nonzero - 1)*maxlvl + 1
    left = np.empty((npoints, ndims), dtype='float64')
    right = np.empty((npoints, ndims), dtype='float64')
    level = np.empty(npoints, dtype='int32')

    # fill zero dims
    left[:,dims_zero] = extent[::2][dims_zero]
    right[:,dims_zero] = extent[1::2][dims_zero]

    # fill non-zero dims
    dcell = 1.0 / cells
    left_slice =  tuple([slice(extent[2*n], extent[2*n+1], extent[2*n+1]) if \
        dims_zero[n] else slice(0.0,1.0,dcell) for n in range(ndims)])
    right_slice = tuple([slice(extent[2*n+1], extent[2*n], -extent[2*n+1]) if \
        dims_zero[n] else slice(dcell,1.0+dcell,dcell) for n in range(ndims)])
    left_norm_grid = np.reshape(np.mgrid[left_slice].T.flat[ndims:], (-1, ndims))
    lng_zero = left_norm_grid[:,dims_zero]
    lng_nonzero = left_norm_grid[:,dims_nonzero]

    right_norm_grid = np.reshape(np.mgrid[right_slice].T.flat[ndims:], (-1, ndims))
    rng_zero = right_norm_grid[:,dims_zero]
    rng_nonzero = right_norm_grid[:,dims_nonzero]

    level[0] = maxlvl
    left[0,:] = extent[::2]
    right[0,dims_zero] = extent[1::2][dims_zero]
    right[0,dims_nonzero] = (dcell**maxlvl)*dextent[dims_nonzero] + extent[::2][dims_nonzero]
    for i, lvl in enumerate(range(maxlvl, 0, -1)):
        start = (cells**ndims_nonzero - 1)*i + 1
        stop = (cells**ndims_nonzero - 1)*(i+1) + 1
        dsize = dcell**(lvl-1) * dextent[dims_nonzero]
        level[start:stop] = lvl
        left[start:stop,dims_zero] = lng_zero
        left[start:stop,dims_nonzero] = lng_nonzero*dsize + extent[::2][dims_nonzero]
        right[start:stop,dims_zero] = rng_zero
        right[start:stop,dims_nonzero] = rng_nonzero*dsize + extent[::2][dims_nonzero]

    return left, right, level

def fake_random_ds(
        ndims, peak_value = 1.0,
        fields = ("density", "velocity_x", "velocity_y", "velocity_z"),
        units = ('g/cm**3', 'cm/s', 'cm/s', 'cm/s'),
        negative = False, nprocs = 1, particles = 0, length_unit=1.0):
    from yt.data_objects.api import data_object_registry
    from yt.frontends.stream.api import load_uniform_grid
    if not iterable(ndims):
        ndims = [ndims, ndims, ndims]
    else:
        assert(len(ndims) == 3)
    if not iterable(negative):
        negative = [negative for f in fields]
    assert(len(fields) == len(negative))
    offsets = []
    for n in negative:
        if n:
            offsets.append(0.5)
        else:
            offsets.append(0.0)
    data = {}
    for field, offset, u in zip(fields, offsets, units):
        v = (np.random.random(ndims) - offset) * peak_value
        if field[0] == "all":
            data['number_of_particles'] = v.size
            v = v.ravel()
        data[field] = (v, u)
    if particles:
        for f in ('particle_position_%s' % ax for ax in 'xyz'):
            data[f] = (np.random.uniform(size = particles), 'code_length')
        for f in ('particle_velocity_%s' % ax for ax in 'xyz'):
            data[f] = (np.random.random(size = particles) - 0.5, 'cm/s')
        data['particle_mass'] = (np.random.random(particles), 'g')
        data['number_of_particles'] = particles
    ug = load_uniform_grid(data, ndims, length_unit=length_unit, nprocs=nprocs)
    return ug

_geom_transforms = {
    # These are the bounds we want.  Cartesian we just assume goes 0 .. 1.
    'cartesian'  : ( (0.0, 0.0, 0.0), (1.0, 1.0, 1.0) ),
    'spherical'  : ( (0.0, 0.0, 0.0), (1.0, np.pi, 2*np.pi) ),
    'cylindrical': ( (0.0, 0.0, 0.0), (1.0, 1.0, 2.0*np.pi) ), # rzt
    'polar'      : ( (0.0, 0.0, 0.0), (1.0, 2.0*np.pi, 1.0) ), # rtz
    'geographic' : ( (-90.0, -180.0, 0.0), (90.0, 180.0, 1000.0) ), # latlonalt
}

def fake_amr_ds(fields = ("Density",), geometry = "cartesian"):
    from yt.frontends.stream.api import load_amr_grids
    LE, RE = _geom_transforms[geometry]
    LE = np.array(LE)
    RE = np.array(RE)
    data = []
    for gspec in _amr_grid_index:
        level, left_edge, right_edge, dims = gspec
        left_edge = left_edge * (RE - LE) + LE
        right_edge = right_edge * (RE - LE) + LE
        gdata = dict(level = level,
                     left_edge = left_edge,
                     right_edge = right_edge,
                     dimensions = dims)
        for f in fields:
            gdata[f] = np.random.random(dims)
        data.append(gdata)
    bbox = np.array([LE, RE]).T
    return load_amr_grids(data, [32, 32, 32], 1.0, geometry=geometry, bbox=bbox)

def fake_particle_ds(
        fields = ("particle_position_x",
                  "particle_position_y",
                  "particle_position_z",
                  "particle_mass",
                  "particle_velocity_x",
                  "particle_velocity_y",
                  "particle_velocity_z"),
        units = ('cm', 'cm', 'cm', 'g', 'cm/s', 'cm/s', 'cm/s'),
        negative = (False, False, False, False, True, True, True),
        npart = 16**3, length_unit=1.0):
    from yt.frontends.stream.api import load_particles
    if not iterable(negative):
        negative = [negative for f in fields]
    assert(len(fields) == len(negative))
    offsets = []
    for n in negative:
        if n:
            offsets.append(0.5)
        else:
            offsets.append(0.0)
    data = {}
    for field, offset, u in zip(fields, offsets, units):
        if "position" in field:
            v = np.random.normal(npart, 0.5, 0.25)
            np.clip(v, 0.0, 1.0, v)
        v = (np.random.random(npart) - offset)
        data[field] = (v, u)
    bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    ds = load_particles(data, 1.0, bbox=bbox)
    return ds

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
    >>> keywords['dpi'] = (50, 100, 200)
    >>> keywords['cmap'] = ('algae', 'jet')
    >>> list_of_kwargs = expand_keywords(keywords)
    >>> print list_of_kwargs

    array([{'cmap': 'algae', 'dpi': 50},
           {'cmap': 'jet', 'dpi': 100},
           {'cmap': 'algae', 'dpi': 200}], dtype=object)

    >>> list_of_kwargs = expand_keywords(keywords, full=True)
    >>> print list_of_kwargs

    array([{'cmap': 'algae', 'dpi': 50},
           {'cmap': 'algae', 'dpi': 100},
           {'cmap': 'algae', 'dpi': 200},
           {'cmap': 'jet', 'dpi': 50},
           {'cmap': 'jet', 'dpi': 100},
           {'cmap': 'jet', 'dpi': 200}], dtype=object)

    >>> for kwargs in list_of_kwargs:
    ...     write_projection(*args, **kwargs)
    """

    # if we want every possible combination of keywords, use iter magic
    if full:
        keys = sorted(keywords)
        list_of_kwarg_dicts = np.array([dict(zip(keys, prod)) for prod in \
                              it.product(*(keywords[key] for key in keys))])

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

def requires_module(module):
    """
    Decorator that takes a module name as an argument and tries to import it.
    If the module imports without issue, the function is returned, but if not,
    a null function is returned. This is so tests that depend on certain modules
    being imported will not fail if the module is not installed on the testing
    platform.
    """
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    try:
        importlib.import_module(module)
    except ImportError:
        return ffalse
    else:
        return ftrue

def requires_file(req_file):
    path = ytcfg.get("yt", "test_data_dir")
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if os.path.exists(req_file):
        return ftrue
    else:
        if os.path.exists(os.path.join(path,req_file)):
            return ftrue
        else:
            return ffalse

def units_override_check(fn):
    ytcfg["yt","skip_dataset_cache"] = "True"
    units_list = ["length","time","mass","velocity",
                  "magnetic","temperature"]
    ds1 = load(fn)
    units_override = {}
    attrs1 = []
    attrs2 = []
    for u in units_list:
        unit_attr = getattr(ds1, "%s_unit" % u, None)
        if unit_attr is not None:
            attrs1.append(unit_attr)
            units_override["%s_unit" % u] = (unit_attr.v, str(unit_attr.units))
    del ds1
    ds2 = load(fn, units_override=units_override)
    ytcfg["yt","skip_dataset_cache"] = "False"
    assert(len(ds2.units_override) > 0)
    for u in units_list:
        unit_attr = getattr(ds2, "%s_unit" % u, None)
        if unit_attr is not None:
            attrs2.append(unit_attr)
    yield assert_equal, attrs1, attrs2

# This is an export of the 40 grids in IsolatedGalaxy that are of level 4 or
# lower.  It's just designed to give a sample AMR index to deal with.
_amr_grid_index = [
 [ 0,
  [0.0,0.0,0.0],
  [1.0,1.0,1.0],
  [32,32,32],
 ],
 [ 1,
  [0.25,0.21875,0.25],
  [0.5,0.5,0.5],
  [16,18,16],
 ],
 [ 1,
  [0.5,0.21875,0.25],
  [0.75,0.5,0.5],
  [16,18,16],
 ],
 [ 1,
  [0.21875,0.5,0.25],
  [0.5,0.75,0.5],
  [18,16,16],
 ],
 [ 1,
  [0.5,0.5,0.25],
  [0.75,0.75,0.5],
  [16,16,16],
 ],
 [ 1,
  [0.25,0.25,0.5],
  [0.5,0.5,0.75],
  [16,16,16],
 ],
 [ 1,
  [0.5,0.25,0.5],
  [0.75,0.5,0.75],
  [16,16,16],
 ],
 [ 1,
  [0.25,0.5,0.5],
  [0.5,0.75,0.75],
  [16,16,16],
 ],
 [ 1,
  [0.5,0.5,0.5],
  [0.75,0.75,0.75],
  [16,16,16],
 ],
 [ 2,
  [0.5,0.5,0.5],
  [0.71875,0.71875,0.71875],
  [28,28,28],
 ],
 [ 3,
  [0.5,0.5,0.5],
  [0.6640625,0.65625,0.6796875],
  [42,40,46],
 ],
 [ 4,
  [0.5,0.5,0.5],
  [0.59765625,0.6015625,0.6015625],
  [50,52,52],
 ],
 [ 2,
  [0.28125,0.5,0.5],
  [0.5,0.734375,0.71875],
  [28,30,28],
 ],
 [ 3,
  [0.3359375,0.5,0.5],
  [0.5,0.671875,0.6640625],
  [42,44,42],
 ],
 [ 4,
  [0.40625,0.5,0.5],
  [0.5,0.59765625,0.59765625],
  [48,50,50],
 ],
 [ 2,
  [0.5,0.28125,0.5],
  [0.71875,0.5,0.71875],
  [28,28,28],
 ],
 [ 3,
  [0.5,0.3359375,0.5],
  [0.671875,0.5,0.6640625],
  [44,42,42],
 ],
 [ 4,
  [0.5,0.40625,0.5],
  [0.6015625,0.5,0.59765625],
  [52,48,50],
 ],
 [ 2,
  [0.28125,0.28125,0.5],
  [0.5,0.5,0.71875],
  [28,28,28],
 ],
 [ 3,
  [0.3359375,0.3359375,0.5],
  [0.5,0.5,0.671875],
  [42,42,44],
 ],
 [ 4,
  [0.46484375,0.37890625,0.50390625],
  [0.4765625,0.390625,0.515625],
  [6,6,6],
 ],
 [ 4,
  [0.40625,0.40625,0.5],
  [0.5,0.5,0.59765625],
  [48,48,50],
 ],
 [ 2,
  [0.5,0.5,0.28125],
  [0.71875,0.71875,0.5],
  [28,28,28],
 ],
 [ 3,
  [0.5,0.5,0.3359375],
  [0.6796875,0.6953125,0.5],
  [46,50,42],
 ],
 [ 4,
  [0.5,0.5,0.40234375],
  [0.59375,0.6015625,0.5],
  [48,52,50],
 ],
 [ 2,
  [0.265625,0.5,0.28125],
  [0.5,0.71875,0.5],
  [30,28,28],
 ],
 [ 3,
  [0.3359375,0.5,0.328125],
  [0.5,0.65625,0.5],
  [42,40,44],
 ],
 [ 4,
  [0.40234375,0.5,0.40625],
  [0.5,0.60546875,0.5],
  [50,54,48],
 ],
 [ 2,
  [0.5,0.265625,0.28125],
  [0.71875,0.5,0.5],
  [28,30,28],
 ],
 [ 3,
  [0.5,0.3203125,0.328125],
  [0.6640625,0.5,0.5],
  [42,46,44],
 ],
 [ 4,
  [0.5,0.3984375,0.40625],
  [0.546875,0.5,0.5],
  [24,52,48],
 ],
 [ 4,
  [0.546875,0.41796875,0.4453125],
  [0.5625,0.4375,0.5],
  [8,10,28],
 ],
 [ 4,
  [0.546875,0.453125,0.41796875],
  [0.5546875,0.48046875,0.4375],
  [4,14,10],
 ],
 [ 4,
  [0.546875,0.4375,0.4375],
  [0.609375,0.5,0.5],
  [32,32,32],
 ],
 [ 4,
  [0.546875,0.4921875,0.41796875],
  [0.56640625,0.5,0.4375],
  [10,4,10],
 ],
 [ 4,
  [0.546875,0.48046875,0.41796875],
  [0.5703125,0.4921875,0.4375],
  [12,6,10],
 ],
 [ 4,
  [0.55859375,0.46875,0.43359375],
  [0.5703125,0.48046875,0.4375],
  [6,6,2],
 ],
 [ 2,
  [0.265625,0.28125,0.28125],
  [0.5,0.5,0.5],
  [30,28,28],
 ],
 [ 3,
  [0.328125,0.3359375,0.328125],
  [0.5,0.5,0.5],
  [44,42,44],
 ],
 [ 4,
  [0.4140625,0.40625,0.40625],
  [0.5,0.5,0.5],
  [44,48,48],
 ],
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

    @check_results
    def my_func(ds):
        return ds.domain_width

    my_func(ds)

    @check_results
    def field_checker(dd, field_name):
        return dd[field_name]

    field_cheker(ds.all_data(), 'density', result_basename='density')

    """
    def compute_results(func):
        def _func(*args, **kwargs):
            name = kwargs.pop("result_basename", func.func_name)
            rv = func(*args, **kwargs)
            if hasattr(rv, "convert_to_cgs"):
                rv.convert_to_cgs()
                _rv = rv.ndarray_view()
            else:
                _rv = rv
            mi = _rv.min()
            ma = _rv.max()
            st = _rv.std(dtype="float64")
            su = _rv.sum(dtype="float64")
            si = _rv.size
            ha = hashlib.md5(_rv.tostring()).hexdigest()
            fn = "func_results_ref_%s.cpkl" % (name)
            with open(fn, "wb") as f:
                cPickle.dump( (mi, ma, st, su, si, ha), f)
            return rv
        return _func
    from yt.mods import unparsed_args
    if "--answer-reference" in unparsed_args:
        return compute_results(func)

    def compare_results(func):
        def _func(*args, **kwargs):
            name = kwargs.pop("result_basename", func.func_name)
            rv = func(*args, **kwargs)
            if hasattr(rv, "convert_to_cgs"):
                rv.convert_to_cgs()
                _rv = rv.ndarray_view()
            else:
                _rv = rv
            vals = (_rv.min(),
                    _rv.max(),
                    _rv.std(dtype="float64"),
                    _rv.sum(dtype="float64"),
                    _rv.size,
                    hashlib.md5(_rv.tostring()).hexdigest() )
            fn = "func_results_ref_%s.cpkl" % (name)
            if not os.path.exists(fn):
                print "Answers need to be created with --answer-reference ."
                return False
            with open(fn, "rb") as f:
                ref = cPickle.load(f)
            print "Sizes: %s (%s, %s)" % (vals[4] == ref[4], vals[4], ref[4])
            assert_allclose(vals[0], ref[0], 1e-8, err_msg="min")
            assert_allclose(vals[1], ref[1], 1e-8, err_msg="max")
            assert_allclose(vals[2], ref[2], 1e-8, err_msg="std")
            assert_allclose(vals[3], ref[3], 1e-8, err_msg="sum")
            assert_equal(vals[4], ref[4])
            print "Hashes equal: %s" % (vals[-1] == ref[-1])
            return rv
        return _func
    return compare_results(func)

def periodicity_cases(ds):
    # This is a generator that yields things near the corners.  It's good for
    # getting different places to check periodicity.
    yield (ds.domain_left_edge + ds.domain_right_edge)/2.0
    dx = ds.domain_width / ds.domain_dimensions
    # We start one dx in, and only go to one in as well.
    for i in (1, ds.domain_dimensions[0] - 2):
        for j in (1, ds.domain_dimensions[1] - 2):
            for k in (1, ds.domain_dimensions[2] - 2):
                center = dx * np.array([i,j,k]) + ds.domain_left_edge
                yield center

def run_nose(verbose=False, run_answer_tests=False, answer_big_data=False):
    import nose, os, sys, yt
    from yt.funcs import mylog
    orig_level = mylog.getEffectiveLevel()
    mylog.setLevel(50)
    nose_argv = sys.argv
    nose_argv += ['--exclude=answer_testing','--detailed-errors', '--exe']
    if verbose:
        nose_argv.append('-v')
    if run_answer_tests:
        nose_argv.append('--with-answer-testing')
    if answer_big_data:
        nose_argv.append('--answer-big-data')
    initial_dir = os.getcwd()
    yt_file = os.path.abspath(yt.__file__)
    yt_dir = os.path.dirname(yt_file)
    os.chdir(yt_dir)
    try:
        nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
        mylog.setLevel(orig_level)
