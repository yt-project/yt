"""
Simple integrators for the radiative transfer equation



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from yt.utilities.lib.fp_utils cimport fclip, i64clip
from libc.math cimport copysign, fabs
from yt.utilities.exceptions import YTDomainOverflow
from yt.utilities.lib.vec3_ops cimport subtract, cross, dot, L2_norm

cdef extern from "math.h":
    double exp(double x) nogil
    float expf(float x) nogil
    long double expl(long double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil
    double fmod(double x, double y) nogil
    double fabs(double x) nogil

cdef extern from "platform_dep.h":
    double log2(double x) nogil
    long int lrint(double x) nogil

# Finally, miscellaneous routines.

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def find_values_at_point(np.ndarray[np.float64_t, ndim=1] point,
                         np.ndarray[np.float64_t, ndim=2] left_edges,
                         np.ndarray[np.float64_t, ndim=2] right_edges,
                         np.ndarray[np.int32_t, ndim=2] dimensions,
                         field_names, grid_objects):
    # This iterates in order, first to last, and then returns with the first
    # one in which the point is located; this means if you order from highest
    # level to lowest, you will find the correct grid without consulting child
    # masking.  Note also that we will do a few relatively slow operations on
    # strings and whatnot, but they should not be terribly slow.
    cdef int ind[3]
    cdef int gi, fi, nf = len(field_names)
    cdef np.float64_t dds
    cdef np.ndarray[np.float64_t, ndim=3] field
    cdef np.ndarray[np.float64_t, ndim=1] rv = np.zeros(nf, dtype='float64')
    for gi in range(left_edges.shape[0]):
        if not ((left_edges[gi,0] < point[0] < right_edges[gi,0])
            and (left_edges[gi,1] < point[1] < right_edges[gi,1])
            and (left_edges[gi,2] < point[2] < right_edges[gi,2])):
            continue
        # We found our grid!
        for fi in range(3):
            dds = ((right_edges[gi,fi] - left_edges[gi,fi])/
                   (<np.float64_t> dimensions[gi,fi]))
            ind[fi] = <int> ((point[fi] - left_edges[gi,fi])/dds)
        grid = grid_objects[gi]
        for fi in range(nf):
            field = grid[field_names[fi]]
            rv[fi] = field[ind[0], ind[1], ind[2]]
        return rv
    raise KeyError

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def obtain_rvec(data, ftype="gas"):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] xf
    cdef np.ndarray[np.float64_t, ndim=1] yf
    cdef np.ndarray[np.float64_t, ndim=1] zf
    cdef np.ndarray[np.float64_t, ndim=2] rf
    cdef np.ndarray[np.float64_t, ndim=3] xg
    cdef np.ndarray[np.float64_t, ndim=3] yg
    cdef np.ndarray[np.float64_t, ndim=3] zg
    cdef np.ndarray[np.float64_t, ndim=4] rg
    cdef np.float64_t c[3]
    cdef int i, j, k
    center = data.get_field_parameter("center")
    c[0] = center[0]; c[1] = center[1]; c[2] = center[2]
    if len(data[ftype, 'x'].shape) == 1:
        # One dimensional data
        xf = data[ftype, 'x']
        yf = data[ftype, 'y']
        zf = data[ftype, 'z']
        rf = np.empty((3, xf.shape[0]), 'float64')
        for i in range(xf.shape[0]):
            rf[0, i] = xf[i] - c[0]
            rf[1, i] = yf[i] - c[1]
            rf[2, i] = zf[i] - c[2]
        return rf
    else:
        # Three dimensional data
        xg = data[ftype, 'x']
        yg = data[ftype, 'y']
        zg = data[ftype, 'z']
        rg = np.empty((3, xg.shape[0], xg.shape[1], xg.shape[2]), 'float64')
        for i in range(xg.shape[0]):
            for j in range(xg.shape[1]):
                for k in range(xg.shape[2]):
                    rg[0,i,j,k] = xg[i,j,k] - c[0]
                    rg[1,i,j,k] = yg[i,j,k] - c[1]
                    rg[2,i,j,k] = zg[i,j,k] - c[2]
        return rg

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t graycode(np.int64_t x):
    return x^(x>>1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t igraycode(np.int64_t x):
    cdef np.int64_t i, j
    if x == 0:
        return x
    m = <np.int64_t> ceil(log2(x)) + 1
    i, j = x, 1
    while j < m:
        i = i ^ (x>>j)
        j += 1
    return i

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t direction(np.int64_t x, np.int64_t n):
    #assert x < 2**n
    if x == 0:
        return 0
    elif x%2 == 0:
        return tsb(x-1, n)%n
    else:
        return tsb(x, n)%n

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t tsb(np.int64_t x, np.int64_t width):
    #assert x < 2**width
    cdef np.int64_t i = 0
    while x&1 and i <= width:
        x = x >> 1
        i += 1
    return i

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t bitrange(np.int64_t x, np.int64_t width,
                         np.int64_t start, np.int64_t end):
    return x >> (width-end) & ((2**(end-start))-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t rrot(np.int64_t x, np.int64_t i, np.int64_t width):
    i = i%width
    x = (x>>i) | (x<<width-i)
    return x&(2**width-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t lrot(np.int64_t x, np.int64_t i, np.int64_t width):
    i = i%width
    x = (x<<i) | (x>>width-i)
    return x&(2**width-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t transform(np.int64_t entry, np.int64_t direction,
                          np.int64_t width, np.int64_t x):
    return rrot((x^entry), direction + 1, width)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t entry(np.int64_t x):
    if x == 0: return 0
    return graycode(2*((x-1)/2))

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t setbit(np.int64_t x, np.int64_t w, np.int64_t i, np.int64_t b):
    if b == 1:
        return x | 2**(w-i-1)
    elif b == 0:
        return x & ~2**(w-i-1)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t point_to_hilbert(int order, np.int64_t p[3]):
    cdef np.int64_t h, e, d, l, b, w, i, x
    h = e = d = 0
    for i in range(order):
        l = 0
        for x in range(3):
            b = bitrange(p[3-x-1], order, i, i+1)
            l |= (b<<x)
        l = transform(e, d, 3, l)
        w = igraycode(l)
        e = e ^ lrot(entry(w), d+1, 3)
        d = (d + direction(w, 3) + 1)%3
        h = (h<<3)|w
    return h

#def hilbert_point(dimension, order, h):
#    """
#        Convert an index on the Hilbert curve of the specified dimension and
#        order to a set of point coordinates.
#    """
#    #    The bit widths in this function are:
#    #        p[*]  - order
#    #        h     - order*dimension
#    #        l     - dimension
#    #        e     - dimension
#    hwidth = order*dimension
#    e, d = 0, 0
#    p = [0]*dimension
#    for i in range(order):
#        w = utils.bitrange(h, hwidth, i*dimension, i*dimension+dimension)
#        l = utils.graycode(w)
#        l = itransform(e, d, dimension, l)
#        for j in range(dimension):
#            b = utils.bitrange(l, dimension, j, j+1)
#            p[j] = utils.setbit(p[j], order, i, b)
#        e = e ^ utils.lrot(entry(w), d+1, dimension)
#        d = (d + direction(w, dimension) + 1)%dimension
#    return p

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void hilbert_to_point(int order, np.int64_t h, np.int64_t *p):
    cdef np.int64_t hwidth, e, d, w, l, b
    cdef int i, j
    hwidth = 3 * order
    e = d = p[0] = p[1] = p[2] = 0
    for i in range(order):
        w = bitrange(h, hwidth, i*3, i*3+3)
        l = graycode(w)
        l = lrot(l, d +1, 3)^e
        for j in range(3):
            b = bitrange(l, 3, j, j+1)
            p[j] = setbit(p[j], order, i, b)
        e = e ^ lrot(entry(w), d+1, 3)
        d = (d + direction(w, 3) + 1)%3

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_hilbert_indices(int order, np.ndarray[np.int64_t, ndim=2] left_index):
    # This is inspired by the scurve package by user cortesi on GH.
    cdef int i
    cdef np.int64_t p[3]
    cdef np.ndarray[np.int64_t, ndim=1] hilbert_indices
    hilbert_indices = np.zeros(left_index.shape[0], 'int64')
    for i in range(left_index.shape[0]):
        p[0] = left_index[i, 0]
        p[1] = left_index[i, 1]
        p[2] = left_index[i, 2]
        hilbert_indices[i] = point_to_hilbert(order, p)
    return hilbert_indices

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_hilbert_points(int order, np.ndarray[np.int64_t, ndim=1] indices):
    # This is inspired by the scurve package by user cortesi on GH.
    cdef int i, j
    cdef np.int64_t p[3]
    cdef np.ndarray[np.int64_t, ndim=2] positions
    positions = np.zeros((indices.shape[0], 3), 'int64')
    for i in range(indices.shape[0]):
        hilbert_to_point(order, indices[i], p)
        for j in range(3):
            positions[i, j] = p[j]
    return positions

# yt did not invent these! :)
cdef np.uint64_t _const20 = 0x000001FFC00003FF
cdef np.uint64_t _const10 = 0x0007E007C00F801F
cdef np.uint64_t _const04 = 0x00786070C0E181C3
cdef np.uint64_t _const2a = 0x0199219243248649
cdef np.uint64_t _const2b = 0x0649249249249249
cdef np.uint64_t _const2c = 0x1249249249249249

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.uint64_t spread_bits(np.uint64_t x):
    # This magic comes from http://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints
    x=(x|(x<<20))&_const20
    x=(x|(x<<10))&_const10
    x=(x|(x<<4))&_const04
    x=(x|(x<<2))&_const2a
    x=(x|(x<<2))&_const2b
    x=(x|(x<<2))&_const2c
    return x

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_morton_indices(np.ndarray[np.uint64_t, ndim=2] left_index):
    cdef np.int64_t i, mi
    cdef np.ndarray[np.uint64_t, ndim=1] morton_indices
    morton_indices = np.zeros(left_index.shape[0], 'uint64')
    for i in range(left_index.shape[0]):
        mi = 0
        mi |= spread_bits(left_index[i,2])<<0
        mi |= spread_bits(left_index[i,1])<<1
        mi |= spread_bits(left_index[i,0])<<2
        morton_indices[i] = mi
    return morton_indices

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_morton_indices_unravel(np.ndarray[np.uint64_t, ndim=1] left_x,
                               np.ndarray[np.uint64_t, ndim=1] left_y,
                               np.ndarray[np.uint64_t, ndim=1] left_z,):
    cdef np.int64_t i, mi
    cdef np.ndarray[np.uint64_t, ndim=1] morton_indices
    morton_indices = np.zeros(left_x.shape[0], 'uint64')
    for i in range(left_x.shape[0]):
        mi = 0
        mi |= spread_bits(left_z[i])<<0
        mi |= spread_bits(left_y[i])<<1
        mi |= spread_bits(left_x[i])<<2
        morton_indices[i] = mi
    return morton_indices

ctypedef fused anyfloat:
    np.float32_t
    np.float64_t

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.int64_t position_to_morton(np.ndarray[anyfloat, ndim=1] pos_x,
                        np.ndarray[anyfloat, ndim=1] pos_y,
                        np.ndarray[anyfloat, ndim=1] pos_z,
                        np.float64_t dds[3], np.float64_t DLE[3],
                        np.float64_t DRE[3],
                        np.ndarray[np.uint64_t, ndim=1] ind,
                        int filter):
    cdef np.uint64_t mi
    cdef np.uint64_t ii[3]
    cdef np.float64_t p[3]
    cdef np.int64_t i, j, use
    cdef np.uint64_t DD[3]
    cdef np.uint64_t FLAG = ~(<np.uint64_t>0)
    for i in range(3):
        DD[i] = <np.uint64_t> ((DRE[i] - DLE[i]) / dds[i])
    for i in range(pos_x.shape[0]):
        use = 1
        p[0] = <np.float64_t> pos_x[i]
        p[1] = <np.float64_t> pos_y[i]
        p[2] = <np.float64_t> pos_z[i]
        for j in range(3):
            if p[j] < DLE[j] or p[j] > DRE[j]:
                if filter == 1:
                    # We only allow 20 levels, so this is inaccessible
                    use = 0
                    break
                return i
            ii[j] = <np.uint64_t> ((p[j] - DLE[j])/dds[j])
            ii[j] = i64clip(ii[j], 0, DD[j] - 1)
        if use == 0:
            ind[i] = FLAG
            continue
        mi = 0
        mi |= spread_bits(ii[2])<<0
        mi |= spread_bits(ii[1])<<1
        mi |= spread_bits(ii[0])<<2
        ind[i] = mi
    return pos_x.shape[0]

DEF ORDER_MAX=20

def compute_morton(np.ndarray pos_x, np.ndarray pos_y, np.ndarray pos_z,
                   domain_left_edge, domain_right_edge, filter_bbox = False):
    cdef int i
    cdef int filter
    if filter_bbox:
        filter = 1
    else:
        filter = 0
    cdef np.float64_t dds[3]
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]
    for i in range(3):
        DLE[i] = domain_left_edge[i]
        DRE[i] = domain_right_edge[i]
        dds[i] = (DRE[i] - DLE[i]) / (1 << ORDER_MAX)
    cdef np.ndarray[np.uint64_t, ndim=1] ind
    ind = np.zeros(pos_x.shape[0], dtype="uint64")
    cdef np.int64_t rv
    if pos_x.dtype == np.float32:
        rv = position_to_morton[np.float32_t](
                pos_x, pos_y, pos_z, dds, DLE, DRE, ind,
                filter)
    elif pos_x.dtype == np.float64:
        rv = position_to_morton[np.float64_t](
                pos_x, pos_y, pos_z, dds, DLE, DRE, ind,
                filter)
    else:
        print "Could not identify dtype.", pos_x.dtype
        raise NotImplementedError
    if rv < pos_x.shape[0]:
        mis = (pos_x.min(), pos_y.min(), pos_z.min())
        mas = (pos_x.max(), pos_y.max(), pos_z.max())
        raise YTDomainOverflow(mis, mas,
                               domain_left_edge, domain_right_edge)
    return ind

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def obtain_rv_vec(data, ftype="gas"):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] vxf
    cdef np.ndarray[np.float64_t, ndim=1] vyf
    cdef np.ndarray[np.float64_t, ndim=1] vzf
    cdef np.ndarray[np.float64_t, ndim=2] rvf
    cdef np.ndarray[np.float64_t, ndim=3] vxg
    cdef np.ndarray[np.float64_t, ndim=3] vyg
    cdef np.ndarray[np.float64_t, ndim=3] vzg
    cdef np.ndarray[np.float64_t, ndim=4] rvg
    cdef np.float64_t bv[3]
    cdef int i, j, k
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = np.zeros(3)
    bv[0] = bulk_velocity[0]; bv[1] = bulk_velocity[1]; bv[2] = bulk_velocity[2]
    if len(data[ftype, 'velocity_x'].shape) == 1:
        # One dimensional data
        vxf = data[ftype, 'velocity_x'].astype("float64")
        vyf = data[ftype, 'velocity_y'].astype("float64")
        vzf = data[ftype, 'velocity_z'].astype("float64")
        rvf = np.empty((3, vxf.shape[0]), 'float64')
        for i in range(vxf.shape[0]):
            rvf[0, i] = vxf[i] - bv[0]
            rvf[1, i] = vyf[i] - bv[1]
            rvf[2, i] = vzf[i] - bv[2]
        return rvf
    else:
        # Three dimensional data
        vxg = data[ftype, 'velocity_x'].astype("float64")
        vyg = data[ftype, 'velocity_y'].astype("float64")
        vzg = data[ftype, 'velocity_z'].astype("float64")
        rvg = np.empty((3, vxg.shape[0], vxg.shape[1], vxg.shape[2]), 'float64')
        for i in range(vxg.shape[0]):
            for j in range(vxg.shape[1]):
                for k in range(vxg.shape[2]):
                    rvg[0,i,j,k] = vxg[i,j,k] - bv[0]
                    rvg[1,i,j,k] = vyg[i,j,k] - bv[1]
                    rvg[2,i,j,k] = vzg[i,j,k] - bv[2]
        return rvg

cdef struct PointSet
cdef struct PointSet:
    int count
    # First index is point index, second is xyz
    np.float64_t points[3][3]
    PointSet *next

cdef inline void get_intersection(np.float64_t p0[3], np.float64_t p1[3],
                                  int ax, np.float64_t coord, PointSet *p):
    cdef np.float64_t vec[3]
    cdef np.float64_t t
    for j in range(3):
        vec[j] = p1[j] - p0[j]
    if vec[ax] == 0.0:
        return  # bail if the line is in the plane
    t = (coord - p0[ax])/vec[ax]
    # We know that if they're on opposite sides, it has to intersect.  And we
    # won't get called otherwise.
    for j in range(3):
        p.points[p.count][j] = p0[j] + vec[j] * t
    p.count += 1

@cython.cdivision(True)
def triangle_plane_intersect(int ax, np.float64_t coord,
                             np.ndarray[np.float64_t, ndim=3] triangles):
    cdef np.float64_t p0[3]
    cdef np.float64_t p1[3]
    cdef np.float64_t p3[3]
    cdef np.float64_t E0[3]
    cdef np.float64_t E1[3]
    cdef np.float64_t tri_norm[3]
    cdef np.float64_t plane_norm[3]
    cdef np.float64_t dp
    cdef int i, j, k, count, ntri, nlines
    nlines = 0
    ntri = triangles.shape[0]
    cdef PointSet *first
    cdef PointSet *last
    cdef PointSet *points
    first = last = points = NULL
    for i in range(ntri):
        count = 0

        # skip if triangle is close to being parallel to plane
        for j in range(3):
            p0[j] = triangles[i, 0, j]
            p1[j] = triangles[i, 1, j]
            p3[j] = triangles[i, 2, j]
            plane_norm[j] = 0.0
        plane_norm[ax] = 1.0
        subtract(p1, p0, E0)
        subtract(p3, p0, E1)
        cross(E0, E1, tri_norm)
        dp = dot(tri_norm, plane_norm)
        dp /= L2_norm(tri_norm)
        if (fabs(dp) > 0.995):
            continue
        
        # Now for each line segment (01, 12, 20) we check to see how many cross
        # the coordinate.
        for j in range(3):
            p3[j] = copysign(1.0, triangles[i, j, ax] - coord)
        if p3[0] * p3[1] < 0: count += 1
        if p3[1] * p3[2] < 0: count += 1
        if p3[2] * p3[0] < 0: count += 1
        if count == 2:
            nlines += 1
        elif count == 3:
            nlines += 2
        else:
            continue
        points = <PointSet *> malloc(sizeof(PointSet))
        points.count = 0
        points.next = NULL
        if p3[0] * p3[1] < 0:
            for j in range(3):
                p0[j] = triangles[i, 0, j]
                p1[j] = triangles[i, 1, j]
            get_intersection(p0, p1, ax, coord, points)
        if p3[1] * p3[2] < 0:
            for j in range(3):
                p0[j] = triangles[i, 1, j]
                p1[j] = triangles[i, 2, j]
            get_intersection(p0, p1, ax, coord, points)
        if p3[2] * p3[0] < 0:
            for j in range(3):
                p0[j] = triangles[i, 2, j]
                p1[j] = triangles[i, 0, j]
            get_intersection(p0, p1, ax, coord, points)
        if last != NULL:
            last.next = points
        if first == NULL:
            first = points
        last = points
    points = first
    cdef np.ndarray[np.float64_t, ndim=3] line_segments
    line_segments = np.empty((nlines, 2, 3), dtype="float64")
    k = 0
    while points != NULL:
        for j in range(3):
            line_segments[k, 0, j] = points.points[0][j]
            line_segments[k, 1, j] = points.points[1][j]
        k += 1
        if points.count == 3:
            for j in range(3):
                line_segments[k, 0, j] = points.points[1][j]
                line_segments[k, 1, j] = points.points[2][j]
            k += 1
        last = points
        points = points.next
        free(last)
    return line_segments
