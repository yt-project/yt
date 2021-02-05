# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Geometry selection routines.




"""


import numpy as np

cimport cython
cimport numpy as np
cimport oct_visitors
from cython cimport floating
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc

from yt.utilities.lib.bitarray cimport ba_get_value, ba_set_value
from yt.utilities.lib.fnv_hash cimport c_fnv_hash as fnv_hash
from yt.utilities.lib.fp_utils cimport fclip, fmax, fmin, iclip, imax, imin
from yt.utilities.lib.geometry_utils cimport (
    bounded_morton_dds,
    decode_morton_64bit,
    encode_morton_64bit,
    morton_neighbors_coarse,
    morton_neighbors_refined,
)
from yt.utilities.lib.grid_traversal cimport sampler_function, walk_volume
from yt.utilities.lib.volume_container cimport VolumeContainer

from .oct_container cimport Oct, OctreeContainer
from .oct_visitors cimport cind


cdef extern from "math.h":
    double exp(double x) nogil
    float expf(float x) nogil
    long double expl(long double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil
    double fmod(double x, double y) nogil
    double log2(double x) nogil
    long int lrint(double x) nogil
    double fabs(double x) nogil

# use this as an epsilon test for grids aligned with selector
# define here to avoid the gil later
cdef np.float64_t grid_eps = np.finfo(np.float64).eps
grid_eps = 0.0

cdef inline np.float64_t dot(np.float64_t* v1,
                             np.float64_t* v2) nogil:
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

cdef inline np.float64_t norm(np.float64_t* v) nogil:
    return sqrt(dot(v, v))

# These routines are separated into a couple different categories:
#
#   * Routines for identifying intersections of an object with a bounding box
#   * Routines for identifying cells/points inside a bounding box that
#     intersect with an object
#   * Routines that speed up some type of geometric calculation

# First, bounding box / object intersection routines.
# These all respect the interface "dobj" and a set of left_edges, right_edges,
# sometimes also accepting level and mask information.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def convert_mask_to_indices(np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
            int count, int transpose = 0):
    cdef int i, j, k, cpos
    cdef np.ndarray[np.int64_t, ndim=2] indices
    indices = np.zeros((count, 3), dtype='int64')
    cpos = 0
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i, j, k] == 1:
                    if transpose == 1:
                        indices[cpos, 0] = k
                        indices[cpos, 1] = j
                        indices[cpos, 2] = i
                    else:
                        indices[cpos, 0] = i
                        indices[cpos, 1] = j
                        indices[cpos, 2] = k
                    cpos += 1
    return indices


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _mask_fill(np.ndarray[np.float64_t, ndim=1] out,
                np.int64_t offset,
                np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
                np.ndarray[floating, ndim=3] vals):
    cdef np.int64_t count = 0
    cdef int i, j, k
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i, j, k] == 1:
                    out[offset + count] = vals[i, j, k]
                    count += 1
    return count

def mask_fill(np.ndarray[np.float64_t, ndim=1] out,
              np.int64_t offset,
              np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
              np.ndarray vals):
    if vals.dtype == np.float32:
        return _mask_fill[np.float32_t](out, offset, mask, vals)
    elif vals.dtype == np.float64:
        return _mask_fill[np.float64_t](out, offset, mask, vals)
    else:
        raise RuntimeError

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def points_in_cells(
        np.float64_t[:] cx,
        np.float64_t[:] cy,
        np.float64_t[:] cz,
        np.float64_t[:] dx,
        np.float64_t[:] dy,
        np.float64_t[:] dz,
        np.float64_t[:] px,
        np.float64_t[:] py,
        np.float64_t[:] pz):
    # Take a list of cells and particles and calculate which particles
    # are enclosed within one of the cells.  This is used for querying
    # particle fields on clump/contour objects.
    # We use brute force since the cells are a relatively unordered collection.

    cdef int p, c, n_p, n_c
    cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask

    n_p = px.size
    n_c = cx.size
    mask = np.zeros(n_p, dtype="bool")

    for p in range(n_p):
        for c in range(n_c):
            if (fabs(px[p] - cx[c]) <= 0.5 * dx[c] and
                fabs(py[p] - cy[c]) <= 0.5 * dy[c] and
                fabs(pz[p] - cz[c]) <= 0.5 * dz[c]):
                mask[p] = True
                break

    return mask

include "_selection_routines/selector_object.pxi"
include "_selection_routines/point_selector.pxi"
include "_selection_routines/sphere_selector.pxi"
include "_selection_routines/region_selector.pxi"
include "_selection_routines/cut_region_selector.pxi"
include "_selection_routines/disk_selector.pxi"
include "_selection_routines/cutting_plane_selector.pxi"
include "_selection_routines/slice_selector.pxi"
include "_selection_routines/ortho_ray_selector.pxi"
include "_selection_routines/ray_selector.pxi"
include "_selection_routines/data_collection_selector.pxi"
include "_selection_routines/ellipsoid_selector.pxi"
include "_selection_routines/grid_selector.pxi"
include "_selection_routines/octree_subset_selector.pxi"
include "_selection_routines/indexed_octree_subset_selector.pxi"
include "_selection_routines/always_selector.pxi"
include "_selection_routines/compose_selector.pxi"
include "_selection_routines/halo_particles_selector.pxi"
include "_selection_routines/boolean_selectors.pxi"
