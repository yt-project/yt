# This module is not part of the public namespace `yt.units`
# It is home to wrapper functions that are directly copied from unyt 2.9.2
# We vendor them as a transition step towards unyt 3.0 (in devlopment),
# where these wrapper functions are deprecated and are should be replaced with vanilla numpy API
# FUTURE:
# - require unyt>=3.0
# - deprecate these functions in yt too

from unyt import unyt_array, unyt_quantity
import numpy as np


def _validate_numpy_wrapper_units(v, arrs):
    if not any(isinstance(a, unyt_array) for a in arrs):
        return v
    if not all(isinstance(a, unyt_array) for a in arrs):
        raise RuntimeError("Not all of your arrays are unyt_arrays.")
    a1 = arrs[0]
    if not all(a.units == a1.units for a in arrs[1:]):
        raise RuntimeError("Your arrays must have identical units.")
    v.units = a1.units
    return v


def uconcatenate(arrs, axis=0):
    """Concatenate a sequence of arrays.

    This wrapper around numpy.concatenate preserves units. All input arrays
    must have the same units.  See the documentation of numpy.concatenate for
    full details.

    Examples
    --------
    >>> from unyt import cm
    >>> A = [1, 2, 3]*cm
    >>> B = [2, 3, 4]*cm
    >>> uconcatenate((A, B))
    unyt_array([1, 2, 3, 2, 3, 4], 'cm')

    """
    v = np.concatenate(arrs, axis=axis)
    v = _validate_numpy_wrapper_units(v, arrs)
    return v


def ucross(arr1, arr2, registry=None, axisa=-1, axisb=-1, axisc=-1, axis=None):
    """Applies the cross product to two YT arrays.

    This wrapper around numpy.cross preserves units.
    See the documentation of numpy.cross for full
    details.
    """
    v = np.cross(arr1, arr2, axisa=axisa, axisb=axisb, axisc=axisc, axis=axis)
    units = arr1.units * arr2.units
    arr = unyt_array(v, units, registry=registry)
    return arr


def uintersect1d(arr1, arr2, assume_unique=False):
    """Find the sorted unique elements of the two input arrays.

    A wrapper around numpy.intersect1d that preserves units.  All input arrays
    must have the same units.  See the documentation of numpy.intersect1d for
    full details.

    Examples
    --------
    >>> from unyt import cm
    >>> A = [1, 2, 3]*cm
    >>> B = [2, 3, 4]*cm
    >>> uintersect1d(A, B)
    unyt_array([2, 3], 'cm')

    """
    v = np.intersect1d(arr1, arr2, assume_unique=assume_unique)
    v = _validate_numpy_wrapper_units(v, [arr1, arr2])
    return v


def uunion1d(arr1, arr2):
    """Find the union of two arrays.

    A wrapper around numpy.intersect1d that preserves units.  All input arrays
    must have the same units.  See the documentation of numpy.intersect1d for
    full details.

    Examples
    --------
    >>> from unyt import cm
    >>> A = [1, 2, 3]*cm
    >>> B = [2, 3, 4]*cm
    >>> uunion1d(A, B)
    unyt_array([1, 2, 3, 4], 'cm')

    """
    v = np.union1d(arr1, arr2)
    v = _validate_numpy_wrapper_units(v, [arr1, arr2])
    return v


def unorm(data, ord=None, axis=None, keepdims=False):
    """Matrix or vector norm that preserves units

    This is a wrapper around np.linalg.norm that preserves units. See
    the documentation for that function for descriptions of the keyword
    arguments.

    Examples
    --------
    >>> from unyt import km
    >>> data = [1, 2, 3]*km
    >>> print(unorm(data))
    3.7416573867739413 km
    """
    norm = np.linalg.norm(data, ord=ord, axis=axis, keepdims=keepdims)
    if norm.shape == ():
        return unyt_quantity(norm, data.units)
    return unyt_array(norm, data.units)


def udot(op1, op2):
    """Matrix or vector dot product that preserves units

    This is a wrapper around np.dot that preserves units.

    Examples
    --------
    >>> from unyt import km, s
    >>> a = np.eye(2)*km
    >>> b = (np.ones((2, 2)) * 2)*s
    >>> print(udot(a, b))
    [[2. 2.]
     [2. 2.]] km*s
    """
    dot = np.dot(op1.d, op2.d)
    units = op1.units * op2.units
    if dot.shape == ():
        return unyt_quantity(dot, units)
    return unyt_array(dot, units)


def uvstack(arrs):
    """Stack arrays in sequence vertically (row wise) while preserving units

    This is a wrapper around np.vstack that preserves units.

    Examples
    --------
    >>> from unyt import km
    >>> a = [1, 2, 3]*km
    >>> b = [2, 3, 4]*km
    >>> print(uvstack([a, b]))
    [[1 2 3]
     [2 3 4]] km
    """
    v = np.vstack(arrs)
    v = _validate_numpy_wrapper_units(v, arrs)
    return v


def uhstack(arrs):
    """Stack arrays in sequence horizontally while preserving units

    This is a wrapper around np.hstack that preserves units.

    Examples
    --------
    >>> from unyt import km
    >>> a = [1, 2, 3]*km
    >>> b = [2, 3, 4]*km
    >>> print(uhstack([a, b]))
    [1 2 3 2 3 4] km
    >>> a = [[1],[2],[3]]*km
    >>> b = [[2],[3],[4]]*km
    >>> print(uhstack([a, b]))
    [[1 2]
     [2 3]
     [3 4]] km
    """
    v = np.hstack(arrs)
    v = _validate_numpy_wrapper_units(v, arrs)
    return v


def ustack(arrs, axis=0):
    """Join a sequence of arrays along a new axis while preserving units

    The axis parameter specifies the index of the new axis in the
    dimensions of the result. For example, if ``axis=0`` it will be the
    first dimension and if ``axis=-1`` it will be the last dimension.

    This is a wrapper around np.stack that preserves units. See the
    documentation for np.stack for full details.

    Examples
    --------
    >>> from unyt import km
    >>> a = [1, 2, 3]*km
    >>> b = [2, 3, 4]*km
    >>> print(ustack([a, b]))
    [[1 2 3]
     [2 3 4]] km
    """
    v = np.stack(arrs, axis=axis)
    v = _validate_numpy_wrapper_units(v, arrs)
    return v
