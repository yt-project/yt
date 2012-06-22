"""
Commonly used mathematical functions.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD Physics/CASS
Author: Stephen Skory <s@skory.us>
Affiliation: UCSD Physics/CASS
Author: Geoffrey So <gsiisg@gmail.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as na
import math

def periodic_dist(a, b, period):
    r"""Find the Euclidian periodic distance between two points.
    
    Parameters
    ----------
    a : array or list
        An array or list of floats.
    
    b : array of list
        An array or list of floats.
    
    period : float or array or list
        If the volume is symmetrically periodic, this can be a single float,
        otherwise an array or list of floats giving the periodic size of the
        volume for each dimension.

    Examples
    --------
    >>> a = na.array([0.1, 0.1, 0.1])
    >>> b = na.array([0.9, 0,9, 0.9])
    >>> period = 1.
    >>> dist = periodic_dist(a, b, 1.)
    >>> dist
    0.3464102
    """
    a = na.array(a)
    b = na.array(b)
    if a.size != b.size: RunTimeError("Arrays must be the same shape.")
    c = na.empty((2, a.size), dtype="float64")
    c[0,:] = abs(a - b)
    c[1,:] = period - abs(a - b)
    d = na.amin(c, axis=0)**2
    return math.sqrt(d.sum())

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
    >>> a = na.array([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], [3, 4, 5]])
    >>> b = rotate_vector_3D(a, 2, na.pi/2)
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
        a = na.array([a])
    if a.shape[1] !=3:
        raise SyntaxError("The second dimension of the array a must be == 3!")
    if dim == 0:
        R = na.array([[1, 0,0],
            [0, na.cos(angle), na.sin(angle)],
            [0, -na.sin(angle), na.cos(angle)]])
    elif dim == 1:
        R = na.array([[na.cos(angle), 0, -na.sin(angle)],
            [0, 1, 0],
            [na.sin(angle), 0, na.cos(angle)]])
    elif dim == 2:
        R = na.array([[na.cos(angle), na.sin(angle), 0],
            [-na.sin(angle), na.cos(angle), 0],
            [0, 0, 1]])
    else:
        raise SyntaxError("dim must be 0, 1, or 2!")
    if mod:
        return na.dot(R, a.T).T[0]
    else:
        return na.dot(R, a.T).T
    

def modify_reference_frame(CoM, L, P, V):
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
        The modified positional data.
    
    V : array
        The modified velocity data.
    
    Examples
    --------
    >>> CoM = na.array([0.5, 0.5, 0.5])
    >>> L = na.array([1, 0, 0])
    >>> P = na.array([[1, 0.5, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0.5], [0, 0, 0]])
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
    if (L == na.array([0, 0, 1.])).all():
        # Whew! Nothing to do!
        return L, P, V
    # First translate the positions to center of mass reference frame.
    P = P - CoM
    # Now find the angle between modified L and the x-axis.
    LL = L.copy()
    LL[2] = 0.
    theta = na.arccos(na.inner(LL, [1.,0,0])/na.inner(LL,LL)**.5)
    if L[1] < 0:
        theta = -theta
    # Now rotate all the position, velocity, and L vectors by this much around
    # the z axis.
    P = rotate_vector_3D(P, 2, theta)
    V = rotate_vector_3D(V, 2, theta)
    L = rotate_vector_3D(L, 2, theta)
    # Now find the angle between L and the z-axis.
    theta = na.arccos(na.inner(L, [0,0,1])/na.inner(L,L)**.5)
    # This time we rotate around the y axis.
    P = rotate_vector_3D(P, 1, theta)
    V = rotate_vector_3D(V, 1, theta)
    L = rotate_vector_3D(L, 1, theta)
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
    >>> CoM = na.array([0, 0, 0])
    >>> L = na.array([0, 0, 1])
    >>> P = na.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = na.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> circV = compute_rotational_velocity(CoM, L, P, V)
    >>> circV
    array([ 1.        ,  0.        ,  0.        ,  1.41421356])

    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # Find the vector in the plane of the galaxy for each position point
    # that is perpendicular to the radial vector.
    radperp = na.cross([0, 0, 1], P)
    # Find the component of the velocity along the radperp vector.
    # Unf., I don't think there's a better way to do this.
    res = na.empty(V.shape[0], dtype='float64')
    for i, rp in enumerate(radperp):
        temp = na.dot(rp, V[i]) / na.dot(rp, rp) * rp
        res[i] = na.dot(temp, temp)**0.5
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
    >>> CoM = na.array([0, 0, 0])
    >>> L = na.array([0, 0, 1])
    >>> P = na.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = na.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
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
    >>> CoM = na.array([0, 0, 0])
    >>> L = na.array([0, 0, 1])
    >>> P = na.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = na.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
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
    res = na.empty(V.shape[0], dtype='float64')
    for i, rad in enumerate(P):
        temp = na.dot(rad, V[i]) / na.dot(rad, rad) * rad
        res[i] = na.dot(temp, temp)**0.5
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
    >>> CoM = na.array([0, 0, 0])
    >>> L = na.array([0, 0, 1])
    >>> P = na.array([[1, 0, 0], [1, 1, 1], [0, 0, 1], [1, 1, 0]])
    >>> V = na.array([[0, 1, 10], [-1, -1, -1], [1, 1, 1], [1, -1, -1]])
    >>> cyl_r = compute_cylindrical_radius(CoM, L, P, V)
    >>> cyl_r
    array([ 1.        ,  1.41421356,  0.        ,  1.41421356])
    """
    # First we translate into the simple coordinates.
    L, P, V = modify_reference_frame(CoM, L, P, V)
    # Demote all the positions to the z=0 plane, which makes the distance
    # calculation very easy.
    P[:,2] = 0
    return na.sqrt((P * P).sum(axis=1))
    
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
    vec1 = na.array(vec1, dtype=na.float64)
    # Normalize
    norm = na.sqrt(na.vdot(vec1, vec1))
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
    vec2 = na.array([x2,y2,z2])
    vec2 /= norm2
    vec3 = na.cross(vec1, vec2)
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
    >>> a = na.arange(100).reshape(10,10)
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
        sorted = na.sort(a, axis=axis)
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
        result.append(na.mean(sorted[indexer], axis=axis, out=out))
    return na.array(result)

def get_rotation_matrix(self, theta, rot_vector):
    ux = rot_vector[0]
    uy = rot_vector[1]
    uz = rot_vector[2]
    cost = na.cos(theta)
    sint = na.sin(theta)
    
    R = na.array([[cost+ux**2*(1-cost), ux*uy*(1-cost)-uz*sint, ux*uz*(1-cost)+uy*sint],
                  [uy*ux*(1-cost)+uz*sint, cost+uy**2*(1-cost), uy*uz*(1-cost)-ux*sint],
                  [uz*ux*(1-cost)-uy*sint, uz*uy*(1-cost)+ux*sint, cost+uz**2*(1-cost)]])
    
    return R

def RX(ax):
    """
    Returns
    -------
    Gives the rotation matrix about the x-axis as an array

    Example
    -------
    >>>>from yt.mods import *
    >>>>from yt.utilities.math_utils import RX
    >>>>RX(na.pi)
    >>>>array([[  1.00000000e+00,   0.00000000e+00,   0.00000000e+00],
               [  0.00000000e+00,  -1.00000000e+00,   1.22464680e-16],
               [  0.00000000e+00,  -1.22464680e-16,  -1.00000000e+00]])
    """
    rot_matrix = na.array([[1, 0, 0], \
                           [0, na.cos(ax), na.sin(ax)], \
                           [0,-na.sin(ax), na.cos(ax)]])
    return rot_matrix
def RY(ay):
    """
    Returns
    -------
    Gives the rotation matrix about the y-axis as an array

    Example
    -------
    >>>>from yt.mods import *
    >>>>from yt.utilities.math_utils import RY
    >>>>RY(na.pi)
    >>>>array([[ -1.00000000e+00,   0.00000000e+00,  -1.22464680e-16],
               [  0.00000000e+00,   1.00000000e+00,   0.00000000e+00],
               [  1.22464680e-16,   0.00000000e+00,  -1.00000000e+00]])
    """
    rot_matrix = na.array([[na.cos(ay), 0,-na.sin(ay)], \
                           [0, 1, 0], \
                           [na.sin(ay), 0, na.cos(ay)]])
    return rot_matrix
def RZ(az):
    """
    Returns
    -------
    Gives the rotation matrix about the z-axis as an array

    Example
    -------
    >>>>from yt.mods import *
    >>>>from yt.utilities.math_utils import RZ
    >>>>RZ(na.pi)
    >>>>array([[ -1.00000000e+00,   1.22464680e-16,   0.00000000e+00],
               [ -1.22464680e-16,  -1.00000000e+00,   0.00000000e+00],
               [  0.00000000e+00,   0.00000000e+00,   1.00000000e+00]])
    """
    rot_matrix = na.array([[na.cos(az), na.sin(az), 0], \
                           [-na.sin(az), na.cos(az), 0], \
                           [0, 0, 1]])
    return rot_matrix
