"""
Commonly used mathematical functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import math
from yt.units.yt_array import \
    YTArray

prec_accum = {
    np.int:                 np.int64,
    np.int8:                np.int64,
    np.int16:               np.int64,
    np.int32:               np.int64,
    np.int64:               np.int64,
    np.uint8:               np.uint64,
    np.uint16:              np.uint64,
    np.uint32:              np.uint64,
    np.uint64:              np.uint64,
    np.float:               np.float64,
    np.float16:             np.float64,
    np.float32:             np.float64,
    np.float64:             np.float64,
    np.complex:             np.complex128,
    np.complex64:           np.complex128,
    np.complex128:          np.complex128,
    np.dtype('int'):        np.int64,
    np.dtype('int8'):       np.int64,
    np.dtype('int16'):      np.int64,
    np.dtype('int32'):      np.int64,
    np.dtype('int64'):      np.int64,
    np.dtype('uint8'):      np.uint64,
    np.dtype('uint16'):     np.uint64,
    np.dtype('uint32'):     np.uint64,
    np.dtype('uint64'):     np.uint64,
    np.dtype('float'):      np.float64,
    np.dtype('float16'):    np.float64,
    np.dtype('float32'):    np.float64,
    np.dtype('float64'):    np.float64,
    np.dtype('complex'):    np.complex128,
    np.dtype('complex64'):  np.complex128,
    np.dtype('complex128'): np.complex128,
}

def periodic_position(pos, ds):
    r"""Assuming periodicity, find the periodic position within the domain.

    Parameters
    ----------
    pos : array
        An array of floats.

    ds : Dataset
        A simulation static output.

    Examples
    --------
    >>> a = np.array([1.1, 0.5, 0.5])
    >>> data = {'Density':np.ones([32,32,32])}
    >>> ds = load_uniform_grid(data, [32,32,32], 1.0)
    >>> ppos = periodic_position(a, ds)
    >>> ppos
    array([ 0.1,  0.5,  0.5])
    """

    off = (pos - ds.domain_left_edge) % ds.domain_width
    return ds.domain_left_edge + off

def periodic_dist(a, b, period, periodicity=(True, True, True)):
    r"""Find the Euclidean periodic distance between two sets of points.

    Parameters
    ----------
    a : array or list
        Either an ndim long list of coordinates corresponding to a single point
        or an (ndim, npoints) list of coordinates for many points in space.

    b : array of list
        Either an ndim long list of coordinates corresponding to a single point
        or an (ndim, npoints) list of coordinates for many points in space.

    period : float or array or list
        If the volume is symmetrically periodic, this can be a single float,
        otherwise an array or list of floats giving the periodic size of the
        volume for each dimension.

    periodicity : An ndim-element tuple of booleans
        If an entry is true, the domain is assumed to be periodic along
        that direction.

    Examples
    --------
    >>> a = [0.1, 0.1, 0.1]
    >>> b = [0.9, 0,9, 0.9]
    >>> period = 1.
    >>> dist = periodic_dist(a, b, 1.)
    >>> dist
    0.346410161514
    """
    a = np.array(a)
    b = np.array(b)
    period = np.array(period)

    if period.size == 1:
        period = np.array([period, period, period])

    if a.shape != b.shape:
        raise RuntimeError("Arrays must be the same shape.")

    if period.shape != a.shape and len(a.shape) > 1:
        n_tup = tuple([1 for i in range(a.ndim-1)])
        period = np.tile(np.reshape(period, (a.shape[0],)+n_tup), (1,)+a.shape[1:])
    elif len(a.shape) == 1:
        a = np.reshape(a, (a.shape[0],)+(1,1))
        b = np.reshape(b, (a.shape[0],)+(1,1))
        period = np.reshape(period, (a.shape[0],)+(1,1))

    c = np.empty((2,) + a.shape, dtype="float64")
    c[0,:] = np.abs(a - b)

    p_directions = [i for i,p in enumerate(periodicity) if p is True]
    np_directions = [i for i,p in enumerate(periodicity) if p is False]
    for d in p_directions:
        c[1,d,:] = period[d,:] - np.abs(a - b)[d,:]
    for d in np_directions:
        c[1,d,:] = c[0,d,:]

    d = np.amin(c, axis=0)**2
    r2 = d.sum(axis=0)
    if r2.size == 1:
        return np.sqrt(r2[0,0])
    return np.sqrt(r2)

def euclidean_dist(a, b):
    r"""Find the Euclidean distance between two points.

    Parameters
    ----------
    a : array or list
        Either an ndim long list of coordinates corresponding to a single point
        or an (ndim, npoints) list of coordinates for many points in space.

    b : array or list
        Either an ndim long list of coordinates corresponding to a single point
        or an (ndim, npoints) list of coordinates for many points in space.

    Examples
    --------
    >>> a = [0.1, 0.1, 0.1]
    >>> b = [0.9, 0,9, 0.9]
    >>> period = 1.
    >>> dist = euclidean_dist(a, b)
    >>> dist
    1.38564064606

    """
    a = np.array(a)
    b = np.array(b)
    if a.shape != b.shape: RuntimeError("Arrays must be the same shape.")
    c = a.copy()
    np.subtract(c, b, c)
    np.power(c, 2, c)
    c = c.sum(axis = 0)
    if isinstance(c, np.ndarray):
        np.sqrt(c, c)
    else:
        # This happens if a and b only have one entry.
        c = math.sqrt(c)
    return c

def rotate_vector_3D(a, dim, angle):
    r"""Rotates the elements of an array around an axis by some angle.

    Given an array of 3D vectors a, this rotates them around a coordinate axis
    by a clockwise angle. An alternative way to think about it is the
    coordinate axes are rotated counterclockwise, which changes the directions
    of the vectors accordingly.

    Parameters
    ----------
    a : array
        An array of 3D vectors with dimension Nx3.

    dim : integer
        A integer giving the axis around which the vectors will be rotated.
        (x, y, z) = (0, 1, 2).

    angle : float
        The angle in radians through which the vectors will be rotated
        clockwise.

    Examples
    --------
    >>> a = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], [3, 4, 5]])
    >>> b = rotate_vector_3D(a, 2, np.pi/2)
    >>> print b
    [[  1.00000000e+00  -1.00000000e+00   0.00000000e+00]
    [  6.12323400e-17  -1.00000000e+00   1.00000000e+00]
    [  1.00000000e+00   6.12323400e-17   1.00000000e+00]
    [  1.00000000e+00  -1.00000000e+00   1.00000000e+00]
    [  4.00000000e+00  -3.00000000e+00   5.00000000e+00]]

    """
    mod = False
    if len(a.shape) == 1:
        mod = True
        a = np.array([a])
    if a.shape[1] !=3:
        raise SyntaxError("The second dimension of the array a must be == 3!")
    if dim == 0:
        R = np.array([[1, 0,0],
            [0, np.cos(angle), np.sin(angle)],
            [0, -np.sin(angle), np.cos(angle)]])
    elif dim == 1:
        R = np.array([[np.cos(angle), 0, -np.sin(angle)],
            [0, 1, 0],
            [np.sin(angle), 0, np.cos(angle)]])
    elif dim == 2:
        R = np.array([[np.cos(angle), np.sin(angle), 0],
            [-np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]])
    else:
        raise SyntaxError("dim must be 0, 1, or 2!")
    if mod:
        return np.dot(R, a.T).T[0]
    else:
        return np.dot(R, a.T).T

def modify_reference_frame(CoM, L, P=None, V=None):
    r"""Rotates and translates data into a new reference frame to make
    calculations easier.

    This is primarily useful for calculations of halo data.
    The data is translated into the center of mass frame.
    Next, it is rotated such that the angular momentum vector for the data
    is aligned with the z-axis. Put another way, if one calculates the angular
    momentum vector on the data that comes out of this function, it will
    always be along the positive z-axis.
    If the center of mass is re-calculated, it will be at the origin.

    Parameters
    ----------
    CoM : array
        The center of mass in 3D.

    L : array
        The angular momentum vector.

    Optional
    --------

    P : array
        The positions of the data to be modified (i.e. particle or grid cell
        postions). The array should be Nx3.

    V : array
        The velocities of the data to be modified (i.e. particle or grid cell
        velocities). The array should be Nx3.

    Returns
    -------
    L : array
        The angular momentum vector equal to [0, 0, 1] modulo machine error.

    P : array
        The modified positional data. Only returned if P is not None

    V : array
        The modified velocity data. Only returned if V is not None

    Examples
    --------
    >>> CoM = np.array([0.5, 0.5, 0.5])
    >>> L = np.array([1, 0, 0])
    >>> P = np.array([[1, 0.5, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0.5], [0, 0, 0]])
    >>> V = p.copy()
    >>> LL, PP, VV = modify_reference_frame(CoM, L, P, V)
    >>> LL
    array([  6.12323400e-17,   0.00000000e+00,   1.00000000e+00])
    >>> PP
    array([[  3.06161700e-17,   0.00000000e+00,   5.00000000e-01],
           [ -3.06161700e-17,   0.00000000e+00,  -5.00000000e-01],
           [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
           [  5.00000000e-01,  -5.00000000e-01,  -5.00000000e-01]])
    >>> VV
    array([[ -5.00000000e-01,   5.00000000e-01,   1.00000000e+00],
           [ -5.00000000e-01,   5.00000000e-01,   3.06161700e-17],
           [ -5.00000000e-01,   5.00000000e-01,   5.00000000e-01],
           [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00]])

    """
    # First translate the positions to center of mass reference frame.
    if P is not None:
        P = P - CoM

    # is the L vector pointing in the Z direction?
    if np.inner(L[:2], L[:2]) == 0.0:
        # the reason why we need to explicitly check for the above
        # is that formula is used in denominator
        # this just checks if we need to flip the z axis or not
        if L[2] < 0.0:
            # this is just a simple flip in direction of the z axis
            if P is not None:
                P = -P
            if V is not None:
                V = -V

        # return the values
        if V is None and P is not None:
            return L, P
        elif P is None and V is not None:
            return L, V
        else:
            return L, P, V

    # Normal vector is not aligned with simulation Z axis
    # Therefore we are going to have to apply a rotation
    # Now find the angle between modified L and the x-axis.
    LL = L.copy()
    LL[2] = 0.0
    theta = np.arccos(np.inner(LL, [1.0, 0.0, 0.0]) / np.inner(LL, LL) ** 0.5)
    if L[1] < 0.0:
        theta = -theta
    # Now rotate all the position, velocity, and L vectors by this much around
    # the z axis.
    if P is not None:
        P = rotate_vector_3D(P, 2, theta)
    if V is not None:
        V = rotate_vector_3D(V, 2, theta)
    L = rotate_vector_3D(L, 2, theta)
    # Now find the angle between L and the z-axis.
    theta = np.arccos(np.inner(L, [0.0, 0.0, 1.0]) / np.inner(L, L) ** 0.5)
    # This time we rotate around the y axis.
    if P is not None:
        P = rotate_vector_3D(P, 1, theta)
    if V is not None:
        V = rotate_vector_3D(V, 1, theta)
    L = rotate_vector_3D(L, 1, theta)

    # return the values
    if V is None and P is not None:
        return L, P
    elif P is None and V is not None:
        return L, V
    else:
        return L, P, V

def compute_rotational_velocity(CoM, L, P, V):
    r"""Computes the rotational velocity for some data around an axis.

    This is primarily for halo computations.
    Given some data, this computes the circular rotational velocity of each
    point (particle) in reference to the axis defined by the angular momentum
    vector.
    This is accomplished by converting the reference frame of the center of
    mass of the halo.

    Parameters
    ----------
    CoM : array
        The center of mass in 3D.

    L : array
        The angular momentum vector.

    P : array
        The positions of the data to be modified (i.e. particle or grid cell
        postions). The array should be Nx3.

    V : array
        The velocities of the data to be modified (i.e. particle or grid cell
        velocities). The array should be Nx3.

    Returns
    -------
    v : array
        An array N elements long that gives the circular rotational velocity
        for each datum (particle).

    Examples
    --------
    >>> CoM = np.array([0, 0, 0])
    >>> L = np.array([0, 0, 1])
    >>> P = np.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = np.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> circV = compute_rotational_velocity(CoM, L, P, V)
    >>> circV
    array([ 1.        ,  0.        ,  0.        ,  1.41421356])

    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # Find the vector in the plane of the galaxy for each position point
    # that is perpendicular to the radial vector.
    radperp = np.cross([0, 0, 1], P)
    # Find the component of the velocity along the radperp vector.
    # Unf., I don't think there's a better way to do this.
    res = np.empty(V.shape[0], dtype='float64')
    for i, rp in enumerate(radperp):
        temp = np.dot(rp, V[i]) / np.dot(rp, rp) * rp
        res[i] = np.dot(temp, temp)**0.5
    return res

def compute_parallel_velocity(CoM, L, P, V):
    r"""Computes the parallel velocity for some data around an axis.

    This is primarily for halo computations.
    Given some data, this computes the velocity component along the angular
    momentum vector.
    This is accomplished by converting the reference frame of the center of
    mass of the halo.

    Parameters
    ----------
    CoM : array
        The center of mass in 3D.

    L : array
        The angular momentum vector.

    P : array
        The positions of the data to be modified (i.e. particle or grid cell
        postions). The array should be Nx3.

    V : array
        The velocities of the data to be modified (i.e. particle or grid cell
        velocities). The array should be Nx3.

    Returns
    -------
    v : array
        An array N elements long that gives the parallel velocity for
        each datum (particle).

    Examples
    --------
    >>> CoM = np.array([0, 0, 0])
    >>> L = np.array([0, 0, 1])
    >>> P = np.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = np.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> paraV = compute_parallel_velocity(CoM, L, P, V)
    >>> paraV
    array([10, -1,  1, -1])

    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # And return just the z-axis velocities.
    return V[:,2]

def compute_radial_velocity(CoM, L, P, V):
    r"""Computes the radial velocity for some data around an axis.

    This is primarily for halo computations.
    Given some data, this computes the radial velocity component for the data.
    This is accomplished by converting the reference frame of the center of
    mass of the halo.

    Parameters
    ----------
    CoM : array
        The center of mass in 3D.

    L : array
        The angular momentum vector.

    P : array
        The positions of the data to be modified (i.e. particle or grid cell
        postions). The array should be Nx3.

    V : array
        The velocities of the data to be modified (i.e. particle or grid cell
        velocities). The array should be Nx3.

    Returns
    -------
    v : array
        An array N elements long that gives the radial velocity for
        each datum (particle).

    Examples
    --------
    >>> CoM = np.array([0, 0, 0])
    >>> L = np.array([0, 0, 1])
    >>> P = np.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = np.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> radV = compute_radial_velocity(CoM, L, P, V)
    >>> radV
    array([ 1.        ,  1.41421356 ,  0.        ,  0.])

    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # We find the tangential velocity by dotting the velocity vector
    # with the cylindrical radial vector for this point.
    # Unf., I don't think there's a better way to do this.
    P[:,2] = 0
    res = np.empty(V.shape[0], dtype='float64')
    for i, rad in enumerate(P):
        temp = np.dot(rad, V[i]) / np.dot(rad, rad) * rad
        res[i] = np.dot(temp, temp)**0.5
    return res

def compute_cylindrical_radius(CoM, L, P, V):
    r"""Compute the radius for some data around an axis in cylindrical
    coordinates.

    This is primarily for halo computations.
    Given some data, this computes the cylindrical radius for each point.
    This is accomplished by converting the reference frame of the center of
    mass of the halo.

    Parameters
    ----------
    CoM : array
        The center of mass in 3D.

    L : array
        The angular momentum vector.

    P : array
        The positions of the data to be modified (i.e. particle or grid cell
        postions). The array should be Nx3.

    V : array
        The velocities of the data to be modified (i.e. particle or grid cell
        velocities). The array should be Nx3.

    Returns
    -------
    cyl_r : array
        An array N elements long that gives the radial velocity for
        each datum (particle).

    Examples
    --------
    >>> CoM = np.array([0, 0, 0])
    >>> L = np.array([0, 0, 1])
    >>> P = np.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = np.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> cyl_r = compute_cylindrical_radius(CoM, L, P, V)
    >>> cyl_r
    array([ 1.        ,  1.41421356,  0.        ,  1.41421356])

    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # Demote all the positions to the z=0 plane, which makes the distance
    # calculation very easy.
    P[:,2] = 0
    return np.sqrt((P * P).sum(axis=1))

def ortho_find(vec1):
    r"""Find two complementary orthonormal vectors to a given vector.

    For any given non-zero vector, there are infinite pairs of vectors
    orthonormal to it.  This function gives you one arbitrary pair from
    that set along with the normalized version of the original vector.

    Parameters
    ----------
    vec1 : array_like
           An array or list to represent a 3-vector.

    Returns
    -------
    vec1 : array
           The original 3-vector array after having been normalized.

    vec2 : array
           A 3-vector array which is orthonormal to vec1.

    vec3 : array
           A 3-vector array which is orthonormal to vec1 and vec2.

    Raises
    ------
    ValueError
           If input vector is the zero vector.

    Notes
    -----
    Our initial vector is `vec1` which consists of 3 components: `x1`, `y1`,
    and `z1`.  ortho_find determines a vector, `vec2`, which is orthonormal
    to `vec1` by finding a vector which has a zero-value dot-product with
    `vec1`.

    .. math::

       vec1 \cdot vec2 = x_1 x_2 + y_1 y_2 + z_1 z_2 = 0

    As a starting point, we arbitrarily choose `vec2` to have `x2` = 1,
    `y2` = 0:

    .. math::

       vec1 \cdot vec2 = x_1 + (z_1 z_2) = 0

       \rightarrow z_2 = -(x_1 / z_1)

    Of course, this will fail if `z1` = 0, in which case, let's say use
    `z2` = 1 and `x2` = 0:

    .. math::

       \rightarrow y_2 = -(z_1 / y_1)

    Similarly, if `y1` = 0, this case will fail, in which case we use
    `y2` = 1 and `z2` = 0:

    .. math::

       \rightarrow x_2 = -(y_1 / x_1)

    Since we don't allow `vec1` to be zero, all cases are accounted for.

    Producing `vec3`, the complementary orthonormal vector to `vec1` and `vec2`
    is accomplished by simply taking the cross product of `vec1` and `vec2`.

    Examples
    --------
    >>> a = [1.0, 2.0, 3.0]
    >>> a, b, c = ortho_find(a)
    >>> a
    array([ 0.26726124,  0.53452248,  0.80178373])
    >>> b
    array([ 0.9486833 ,  0.        , -0.31622777])
    >>> c
    array([-0.16903085,  0.84515425, -0.50709255])
    """
    vec1 = np.array(vec1, dtype=np.float64)
    # Normalize
    norm = np.sqrt(np.vdot(vec1, vec1))
    if norm == 0:
        raise ValueError("Zero vector used as input.")
    vec1 /= norm
    x1 = vec1[0]
    y1 = vec1[1]
    z1 = vec1[2]
    if z1 != 0:
        x2 = 1.0
        y2 = 0.0
        z2 = -(x1 / z1)
        norm2 = (1.0 + z2 ** 2.0) ** (0.5)
    elif y1 != 0:
        x2 = 0.0
        z2 = 1.0
        y2 = -(z1 / y1)
        norm2 = (1.0 + y2 ** 2.0) ** (0.5)
    else:
        y2 = 1.0
        z2 = 0.0
        x2 = -(y1 / x1)
        norm2 = (1.0 + z2 ** 2.0) ** (0.5)
    vec2 = np.array([x2,y2,z2])
    vec2 /= norm2
    vec3 = np.cross(vec1, vec2)
    return vec1, vec2, vec3

def quartiles(a, axis=None, out=None, overwrite_input=False):
    """
    Compute the quartile values (25% and 75%) along the specified axis
    in the same way that the numpy.median calculates the median (50%) value
    alone a specified axis.  Check numpy.median for details, as it is
    virtually the same algorithm.

    Returns an array of the quartiles of the array elements [lower quartile,
    upper quartile].

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : {None, int}, optional
        Axis along which the quartiles are computed. The default (axis=None)
        is to compute the quartiles along a flattened version of the array.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output,
        but the type (of the output) will be cast if necessary.
    overwrite_input : {False, True}, optional
       If True, then allow use of memory of input array (a) for
       calculations. The input array will be modified by the call to
       quartiles. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted. Default is
       False. Note that, if `overwrite_input` is True and the input
       is not already an ndarray, an error will be raised.

    Returns
    -------
    quartiles : ndarray
        A new 2D array holding the result (unless `out` is specified, in
        which case that array is returned instead).  If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    See Also
    --------
    numpy.median, numpy.mean, numpy.percentile

    Notes
    -----
    Given a vector V of length N, the quartiles of V are the 25% and 75% values
    of a sorted copy of V, ``V_sorted`` - i.e., ``V_sorted[(N-1)/4]`` and
    ``3*V_sorted[(N-1)/4]``, when N is odd.  When N is even, it is the average
    of the two values bounding these values of ``V_sorted``.

    Examples
    --------
    >>> a = np.arange(100).reshape(10,10)
    >>> a
    array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
           [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
           [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
           [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
           [50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
           [60, 61, 62, 63, 64, 65, 66, 67, 68, 69],
           [70, 71, 72, 73, 74, 75, 76, 77, 78, 79],
           [80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
           [90, 91, 92, 93, 94, 95, 96, 97, 98, 99]])
    >>> mu.quartiles(a)
    array([ 24.5,  74.5])
    >>> mu.quartiles(a,axis=0)
    array([[ 15.,  16.,  17.,  18.,  19.,  20.,  21.,  22.,  23.,  24.],
           [ 65.,  66.,  67.,  68.,  69.,  70.,  71.,  72.,  73.,  74.]])
    >>> mu.quartiles(a,axis=1)
    array([[  1.5,  11.5,  21.5,  31.5,  41.5,  51.5,  61.5,  71.5,  81.5,
             91.5],
           [  6.5,  16.5,  26.5,  36.5,  46.5,  56.5,  66.5,  76.5,  86.5,
             96.5]])
    """
    if overwrite_input:
        if axis is None:
            sorted = a.ravel()
            sorted.sort()
        else:
            a.sort(axis=axis)
            sorted = a
    else:
        sorted = np.sort(a, axis=axis)
    if axis is None:
        axis = 0
    indexer = [slice(None)] * sorted.ndim
    indices = [int(sorted.shape[axis]/4), int(sorted.shape[axis]*.75)]
    result = []
    for index in indices:
        if sorted.shape[axis] % 2 == 1:
            # index with slice to allow mean (below) to work
            indexer[axis] = slice(index, index+1)
        else:
            indexer[axis] = slice(index-1, index+1)
        # special cases for small arrays
        if sorted.shape[axis] == 2:
            # index with slice to allow mean (below) to work
            indexer[axis] = slice(index, index+1)
        # Use mean in odd and even case to coerce data type
        # and check, use out array.
        result.append(np.mean(sorted[indexer], axis=axis, out=out))
    return np.array(result)

def get_perspective_matrix(fovy, aspect, z_near, z_far):
    """
    Given a field of view in radians, an aspect ratio, and a near
    and far plane distance, this routine computes the transformation matrix
    corresponding to perspective projection using homogenous coordinates.

    Parameters
    ----------
    fovy : scalar
        The angle in degrees of the field of view.

    aspect : scalar
        The aspect ratio of width / height for the projection.

    z_near : scalar
        The distance of the near plane from the camera.

    z_far : scalar
        The distance of the far plane from the camera.

    Returns
    -------
    persp_matrix : ndarray
        A new 4x4 2D array. Represents a perspective transformation
        in homogeneous coordinates. Note that this matrix does not
        actually perform the projection. After multiplying a 4D
        vector of the form (x_0, y_0, z_0, 1.0), the point will be
        transformed to some (x_1, y_1, z_1, w). The final projection
        is applied by performing a divide by w, that is
        (x_1/w, y_1/w, z_1/w, w/w). The matrix uses a row-major
        ordering, rather than the column major ordering typically
        used by OpenGL.

    Notes
    -----
    The usage of 4D homogeneous coordinates is for OpenGL and GPU
    hardware that automatically performs the divide by w operation.
    See the following for more details about the OpenGL perpective matrices.

    http://www.tomdalling.com/blog/modern-opengl/explaining-homogenous-coordinates-and-projective-geometry/
    http://www.songho.ca/opengl/gl_projectionmatrix.html

    """

    tan_half_fovy = np.tan(np.radians(fovy) / 2)

    result = np.zeros( (4, 4), dtype = 'float32', order = 'C')
    #result[0][0] = 1 / (aspect * tan_half_fovy)
    #result[1][1] = 1 / tan_half_fovy
    #result[2][2] = - (z_far + z_near) / (z_far - z_near)
    #result[3][2] = -1
    #result[2][3] = -(2 * z_far * z_near) / (z_far - z_near)

    f = z_far
    n = z_near

    t = tan_half_fovy * n
    b = -t * aspect
    r = t * aspect
    l = - t *  aspect

    result[0][0] = (2 * n) / (r - l)
    result[2][0] = (r + l) / (r - l)
    result[1][1] = (2 * n) / (t - b)
    result[1][2] = (t + b) / (t - b)
    result[2][2] = -(f + n) / (f - n)
    result[2][3] = -2*f*n/(f - n)
    result[3][2] = -1

    return result

def get_orthographic_matrix(maxr, aspect, z_near, z_far):
    """
    Given a field of view in radians, an aspect ratio, and a near
    and far plane distance, this routine computes the transformation matrix
    corresponding to perspective projection using homogenous coordinates.

    Parameters
    ----------
    maxr : scalar
        should be max(|x|, |y|)

    aspect : scalar
        The aspect ratio of width / height for the projection.

    z_near : scalar
        The distance of the near plane from the camera.

    z_far : scalar
        The distance of the far plane from the camera.

    Returns
    -------
    persp_matrix : ndarray
        A new 4x4 2D array. Represents a perspective transformation
        in homogeneous coordinates. Note that this matrix does not
        actually perform the projection. After multiplying a 4D
        vector of the form (x_0, y_0, z_0, 1.0), the point will be
        transformed to some (x_1, y_1, z_1, w). The final projection
        is applied by performing a divide by w, that is
        (x_1/w, y_1/w, z_1/w, w/w). The matrix uses a row-major
        ordering, rather than the column major ordering typically
        used by OpenGL.

    Notes
    -----
    The usage of 4D homogeneous coordinates is for OpenGL and GPU
    hardware that automatically performs the divide by w operation.
    See the following for more details about the OpenGL perpective matrices.

    http://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/orthographic-projection-matrix
    http://www.tomdalling.com/blog/modern-opengl/explaining-homogenous-coordinates-and-projective-geometry/
    http://www.songho.ca/opengl/gl_projectionmatrix.html

    """

    r = maxr * aspect
    t = maxr
    l = -r
    b = -t

    result = np.zeros( (4, 4), dtype = 'float32', order = 'C')
    result[0][0] = 2.0 / (r - l)
    result[1][1] = 2.0 / (t - b)
    result[2][2] = -2.0 / (z_far - z_near)
    result[3][3] = 1

    result[3][0] = - (r+l)/(r-l)
    result[3][1] = -(t+b)/(t-b)
    result[3][2] = -(z_far + z_near) / (z_far - z_near)

    return result

def get_lookat_matrix(eye, center, up):
    """
    Given the position of a camera, the point it is looking at, and
    an up-direction. Computes the lookat matrix that moves all vectors
    such that the camera is at the origin of the coordinate system,
    looking down the z-axis.

    Parameters
    ----------
    eye : array_like
        The position of the camera. Must be 3D.

    center : array_like
        The location that the camera is looking at. Must be 3D.

    up : array_like
        The direction that is considered up for the camera. Must be
        3D.

    Returns
    -------
    lookat_matrix : ndarray
        A new 4x4 2D array in homogeneous coordinates. This matrix
        moves all vectors in the same way required to move the camera
        to the origin of the coordinate system, with it pointing down
        the negative z-axis.

    """

    eye = np.array(eye)
    center = np.array(center)
    up = np.array(up)

    f = (center - eye) / np.linalg.norm(center - eye)
    s = np.cross(f, up) / np.linalg.norm(np.cross(f, up))
    u = np.cross(s, f)

    result = np.zeros ( (4, 4), dtype = 'float32', order = 'C')

    result[0][0] = s[0]
    result[0][1] = s[1]
    result[0][2] = s[2]
    result[1][0] = u[0]
    result[1][1] = u[1]
    result[1][2] = u[2]
    result[2][0] =-f[0]
    result[2][1] =-f[1]
    result[2][2] =-f[2]
    result[0][3] =-np.dot(s, eye)
    result[1][3] =-np.dot(u, eye)
    result[2][3] = np.dot(f, eye)
    result[3][3] = 1.0
    return result


def get_translate_matrix(dx, dy, dz):
    """
    Given a movement amount for each coordinate, creates a translation
    matrix that moves the vector by each amount.

    Parameters
    ----------
    dx : scalar
        A translation amount for the x-coordinate

    dy : scalar
        A translation amount for the y-coordinate

    dz : scalar
        A translation amount for the z-coordinate

    Returns
    -------
    trans_matrix : ndarray
        A new 4x4 2D array. Represents a translation by dx, dy
        and dz in each coordinate respectively.
    """
    result = np.zeros( (4, 4), dtype = 'float32', order = 'C')

    result[0][0] = 1.0
    result[1][1] = 1.0
    result[2][2] = 1.0
    result[3][3] = 1.0

    result[0][3] = dx
    result[1][3] = dy
    result[2][3] = dz

    return result

def get_scale_matrix(dx, dy, dz):
    """
    Given a scaling factor for each coordinate, returns a matrix that
    corresponds to the given scaling amounts.

    Parameters
    ----------
    dx : scalar
        A scaling factor for the x-coordinate.

    dy : scalar
        A scaling factor for the y-coordinate.

    dz : scalar
        A scaling factor for the z-coordinate.

    Returns
    -------
    scale_matrix : ndarray
        A new 4x4 2D array. Represents a scaling by dx, dy, and dz
        in each coordinate respectively.
    """
    result = np.zeros( (4, 4), dtype = 'float32', order = 'C')

    result[0][0] = dx
    result[1][1] = dy
    result[2][2] = dz
    result[3][3] = 1

    return result

def get_rotation_matrix(theta, rot_vector):
    """
    Given an angle theta and a 3D vector rot_vector, this routine
    computes the rotation matrix corresponding to rotating theta
    radians about rot_vector.

    Parameters
    ----------
    theta : scalar
        The angle in radians.

    rot_vector : array_like
        The axis of rotation.  Must be 3D.

    Returns
    -------
    rot_matrix : ndarray
         A new 3x3 2D array.  This is the representation of a
         rotation of theta radians about rot_vector in the simulation
         box coordinate frame

    See Also
    --------
    ortho_find

    Examples
    --------
    >>> a = [0,1,0]
    >>> theta = 0.785398163  # pi/4
    >>> rot = mu.get_rotation_matrix(theta,a)
    >>> rot
    array([[ 0.70710678,  0.        ,  0.70710678],
           [ 0.        ,  1.        ,  0.        ],
           [-0.70710678,  0.        ,  0.70710678]])
    >>> np.dot(rot,a)
    array([ 0.,  1.,  0.])
    # since a is an eigenvector by construction
    >>> np.dot(rot,[1,0,0])
    array([ 0.70710678,  0.        , -0.70710678])
    """

    ux = rot_vector[0]
    uy = rot_vector[1]
    uz = rot_vector[2]
    cost = np.cos(theta)
    sint = np.sin(theta)

    R = np.array([[cost+ux**2*(1-cost), ux*uy*(1-cost)-uz*sint, ux*uz*(1-cost)+uy*sint],
                  [uy*ux*(1-cost)+uz*sint, cost+uy**2*(1-cost), uy*uz*(1-cost)-ux*sint],
                  [uz*ux*(1-cost)-uy*sint, uz*uy*(1-cost)+ux*sint, cost+uz**2*(1-cost)]])

    return R

def quaternion_mult(q1, q2):
    '''

    Multiply two quaternions. The inputs are 4-component numpy arrays
    in the order [w, x, y, z].

    '''
    w = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3]
    x = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2]
    y = q1[0]*q2[2] + q1[2]*q2[0] + q1[3]*q2[1] - q1[1]*q2[3]
    z = q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1]
    return np.array([w, x, y, z])

def quaternion_to_rotation_matrix(quaternion):
    """

    This converts a quaternion representation of on orientation to
    a rotation matrix. The input is a 4-component numpy array in
    the order [w, x, y, z], and the output is a 3x3 matrix stored
    as a 2D numpy array.  We follow the approach in
    "3D Math Primer for Graphics and Game Development" by
    Dunn and Parberry.

    """

    w = quaternion[0]
    x = quaternion[1]
    y = quaternion[2]
    z = quaternion[3]

    R = np.empty((3, 3), dtype=np.float64)

    R[0][0] = 1.0 - 2.0*y**2 - 2.0*z**2
    R[0][1] = 2.0*x*y + 2.0*w*z
    R[0][2] = 2.0*x*z - 2.0*w*y

    R[1][0] = 2.0*x*y - 2.0*w*z
    R[1][1] = 1.0 - 2.0*x**2 - 2.0*z**2
    R[1][2] = 2.0*y*z + 2.0*w*x

    R[2][0] = 2.0*x*z + 2.0*w*y
    R[2][1] = 2.0*y*z - 2.0*w*x
    R[2][2] = 1.0 - 2.0*x**2 - 2.0*y**2

    return R

def rotation_matrix_to_quaternion(rot_matrix):
    '''

    Convert a rotation matrix-based representation of an
    orientation to a quaternion. The input should be a
    3x3 rotation matrix, while the output will be a
    4-component numpy array. We follow the approach in
    "3D Math Primer for Graphics and Game Development" by
    Dunn and Parberry.

    '''
    m11 = rot_matrix[0][0]
    m12 = rot_matrix[0][1]
    m13 = rot_matrix[0][2]
    m21 = rot_matrix[1][0]
    m22 = rot_matrix[1][1]
    m23 = rot_matrix[1][2]
    m31 = rot_matrix[2][0]
    m32 = rot_matrix[2][1]
    m33 = rot_matrix[2][2]

    four_w_squared_minus_1 = m11 + m22 + m33
    four_x_squared_minus_1 = m11 - m22 - m33
    four_y_squared_minus_1 = m22 - m11 - m33
    four_z_squared_minus_1 = m33 - m11 - m22
    max_index = 0
    four_max_squared_minus_1 = four_w_squared_minus_1
    if (four_x_squared_minus_1 > four_max_squared_minus_1):
        four_max_squared_minus_1 = four_x_squared_minus_1
        max_index = 1
    if (four_y_squared_minus_1 > four_max_squared_minus_1):
        four_max_squared_minus_1 = four_y_squared_minus_1
        max_index = 2
    if (four_z_squared_minus_1 > four_max_squared_minus_1):
        four_max_squared_minus_1 = four_z_squared_minus_1
        max_index = 3

    max_val = 0.5*np.sqrt(four_max_squared_minus_1 + 1.0)
    mult = 0.25 / max_val

    if (max_index == 0):
        w = max_val
        x = (m23 - m32) * mult
        y = (m31 - m13) * mult
        z = (m12 - m21) * mult
    elif (max_index == 1):
        x = max_val
        w = (m23 - m32) * mult
        y = (m12 + m21) * mult
        z = (m31 + m13) * mult
    elif (max_index == 2):
        y = max_val
        w = (m31 - m13) * mult
        x = (m12 + m21) * mult
        z = (m23 + m32) * mult
    elif (max_index == 3):
        z = max_val
        w = (m12 - m21) * mult
        x = (m31 + m13) * mult
        y = (m23 + m32) * mult

    return np.array([w, x, y, z])

def get_ortho_basis(normal):
    xprime = np.cross([0.0,1.0,0.0],normal)
    if np.sum(xprime) == 0: xprime = np.array([0.0, 0.0, 1.0])
    yprime = np.cross(normal,xprime)
    zprime = normal
    return (xprime, yprime, zprime)

def get_sph_r(coords):
    # The spherical coordinates radius is simply the magnitude of the
    # coordinate vector.

    return np.sqrt(np.sum(coords**2, axis=0))

def resize_vector(vector,vector_array):
    if len(vector_array.shape) == 4:
        res_vector = np.resize(vector,(3,1,1,1))
    else:
        res_vector = np.resize(vector,(3,1))
    return res_vector
    
def normalize_vector(vector):
    # this function normalizes
    # an input vector
    
    L2 = np.atleast_1d(np.linalg.norm(vector))
    L2[L2==0] = 1.0
    vector = vector / L2
    return vector

def get_sph_theta(coords, normal):
    # The angle (theta) with respect to the normal (J), is the arccos
    # of the dot product of the normal with the normalized coordinate
    # vector.

    res_normal = resize_vector(normal, coords)

    # check if the normal vector is normalized
    # since arccos requires the vector to be normalised
    res_normal = normalize_vector(res_normal)

    tile_shape = [1] + list(coords.shape)[1:]

    J = np.tile(res_normal,tile_shape)

    JdotCoords = np.sum(J*coords,axis=0)

    return np.arccos( JdotCoords / np.sqrt(np.sum(coords**2,axis=0)) )

def get_sph_phi(coords, normal):
    # We have freedom with respect to what axis (xprime) to define
    # the disk angle. Here I've chosen to use the axis that is
    # perpendicular to the normal and the y-axis. When normal ==
    # y-hat, then set xprime = z-hat. With this definition, when
    # normal == z-hat (as is typical), then xprime == x-hat.
    #
    # The angle is then given by the arctan of the ratio of the
    # yprime-component and the xprime-component of the coordinate
    # vector.

    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, coords)
    res_yprime = resize_vector(yprime, coords)

    tile_shape = [1] + list(coords.shape)[1:]
    Jx = np.tile(res_xprime,tile_shape)
    Jy = np.tile(res_yprime,tile_shape)

    Px = np.sum(Jx*coords,axis=0)
    Py = np.sum(Jy*coords,axis=0)

    return np.arctan2(Py,Px)

def get_cyl_r(coords, normal):
    # The cross product of the normal (J) with a coordinate vector
    # gives a vector of magnitude equal to the cylindrical radius.

    res_normal = resize_vector(normal, coords)
    res_normal = normalize_vector(res_normal)

    tile_shape = [1] + list(coords.shape)[1:]
    J = np.tile(res_normal, tile_shape)

    JcrossCoords = np.cross(J, coords, axisa=0, axisb=0, axisc=0)
    return np.sqrt(np.sum(JcrossCoords**2, axis=0))

def get_cyl_z(coords, normal):
    # The dot product of the normal (J) with the coordinate vector
    # gives the cylindrical height.

    res_normal = resize_vector(normal, coords)
    res_normal = normalize_vector(res_normal)
    
    tile_shape = [1] + list(coords.shape)[1:]
    J = np.tile(res_normal, tile_shape)

    return np.sum(J*coords, axis=0)

def get_cyl_theta(coords, normal):
    # This is identical to the spherical phi component

    return get_sph_phi(coords, normal)


def get_cyl_r_component(vectors, theta, normal):
    # The r of a vector is the vector dotted with rhat

    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, vectors)
    res_yprime = resize_vector(yprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]
    Jx = np.tile(res_xprime,tile_shape)
    Jy = np.tile(res_yprime,tile_shape)

    rhat = Jx*np.cos(theta) + Jy*np.sin(theta)

    return np.sum(vectors*rhat,axis=0)

def get_cyl_theta_component(vectors, theta, normal):
    # The theta component of a vector is the vector dotted with thetahat
    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, vectors)
    res_yprime = resize_vector(yprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]
    Jx = np.tile(res_xprime,tile_shape)
    Jy = np.tile(res_yprime,tile_shape)

    thetahat = -Jx*np.sin(theta) + Jy*np.cos(theta)

    return np.sum(vectors*thetahat, axis=0)

def get_cyl_z_component(vectors, normal):
    # The z component of a vector is the vector dotted with zhat
    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_zprime = resize_vector(zprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]
    zhat = np.tile(res_zprime, tile_shape)

    return np.sum(vectors*zhat, axis=0)

def get_sph_r_component(vectors, theta, phi, normal):
    # The r component of a vector is the vector dotted with rhat
    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, vectors)
    res_yprime = resize_vector(yprime, vectors)
    res_zprime = resize_vector(zprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]

    Jx, Jy, Jz = (
        YTArray(np.tile(rprime, tile_shape), "")
        for rprime in (res_xprime, res_yprime, res_zprime))

    rhat = Jx*np.sin(theta)*np.cos(phi) + \
           Jy*np.sin(theta)*np.sin(phi) + \
           Jz*np.cos(theta)

    return np.sum(vectors*rhat, axis=0)

def get_sph_phi_component(vectors, phi, normal):
    # The phi component of a vector is the vector dotted with phihat
    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, vectors)
    res_yprime = resize_vector(yprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]
    Jx = YTArray(np.tile(res_xprime,tile_shape), "")
    Jy = YTArray(np.tile(res_yprime,tile_shape), "")

    phihat = -Jx*np.sin(phi) + Jy*np.cos(phi)

    return np.sum(vectors*phihat, axis=0)

def get_sph_theta_component(vectors, theta, phi, normal):
    # The theta component of a vector is the vector dotted with thetahat
    normal = normalize_vector(normal)
    (xprime, yprime, zprime) = get_ortho_basis(normal)

    res_xprime = resize_vector(xprime, vectors)
    res_yprime = resize_vector(yprime, vectors)
    res_zprime = resize_vector(zprime, vectors)

    tile_shape = [1] + list(vectors.shape)[1:]
    Jx, Jy, Jz = (
        YTArray(np.tile(rprime, tile_shape), "")
        for rprime in (res_xprime, res_yprime, res_zprime))


    thetahat = Jx*np.cos(theta)*np.cos(phi) + \
               Jy*np.cos(theta)*np.sin(phi) - \
               Jz*np.sin(theta)

    return np.sum(vectors*thetahat, axis=0)
