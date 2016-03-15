"""
Marching cubes implementation



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
cimport cython
import numpy as np
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip
from libc.stdlib cimport malloc, free, abs
from fixed_interpolator cimport *

cdef extern from "marching_cubes.h":
    int tri_table[256][16]
    int edge_table[256]

cdef struct Triangle:
    Triangle *next
    np.float64_t p[3][3]
    np.float64_t val[3] # Usually only use one value

cdef struct TriangleCollection:
    int count
    Triangle *first
    Triangle *current

cdef Triangle *AddTriangle(Triangle *self,
                    np.float64_t p0[3], np.float64_t p1[3], np.float64_t p2[3]):
    cdef Triangle *nn = <Triangle *> malloc(sizeof(Triangle))
    if self != NULL:
        self.next = nn
    cdef int i
    for i in range(3):
        nn.p[0][i] = p0[i]
    for i in range(3):
        nn.p[1][i] = p1[i]
    for i in range(3):
        nn.p[2][i] = p2[i]
    nn.next = NULL
    return nn

cdef int CountTriangles(Triangle *first):
    cdef int count = 0
    cdef Triangle *this = first
    while this != NULL:
        count += 1
        this = this.next
    return count

cdef void FillTriangleValues(np.ndarray[np.float64_t, ndim=1] values,
                             Triangle *first, int nskip = 1):
    cdef Triangle *this = first
    cdef int i = 0
    cdef int j
    while this != NULL:
        for j in range(nskip):
            values[i*nskip + j] = this.val[j]
        i += 1
        this = this.next

cdef void WipeTriangles(Triangle *first):
    cdef Triangle *this = first
    cdef Triangle *last
    while this != NULL:
        last = this
        this = this.next
        free(last)

cdef void FillAndWipeTriangles(np.ndarray[np.float64_t, ndim=2] vertices,
                               Triangle *first):
    cdef int count = 0
    cdef Triangle *this = first
    cdef Triangle *last
    cdef int i, j
    while this != NULL:
        for i in range(3):
            for j in range(3):
                vertices[count, j] = this.p[i][j]
            count += 1 # Do it at the end because it's an index
        last = this
        this = this.next
        free(last)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int march_cubes(
                 np.float64_t gv[8], np.float64_t isovalue,
                 np.float64_t dds[3],
                 np.float64_t x, np.float64_t y, np.float64_t z,
                 TriangleCollection *triangles):

    cdef np.float64_t vertlist[12][3]
    cdef int cubeindex = 0
    cdef int n
    cdef int nt = 0
    for n in range(8):
        if gv[n] < isovalue:
            cubeindex |= (1 << n)
    if edge_table[cubeindex] == 0:
        return 0
    if (edge_table[cubeindex] & 1): # 0,0,0 with 1,0,0
        vertex_interp(gv[0], gv[1], isovalue, vertlist[0],
                      dds, x, y, z, 0, 1)
    if (edge_table[cubeindex] & 2): # 1,0,0 with 1,1,0
        vertex_interp(gv[1], gv[2], isovalue, vertlist[1],
                      dds, x, y, z, 1, 2)
    if (edge_table[cubeindex] & 4): # 1,1,0 with 0,1,0
        vertex_interp(gv[2], gv[3], isovalue, vertlist[2],
                      dds, x, y, z, 2, 3)
    if (edge_table[cubeindex] & 8): # 0,1,0 with 0,0,0
        vertex_interp(gv[3], gv[0], isovalue, vertlist[3],
                      dds, x, y, z, 3, 0)
    if (edge_table[cubeindex] & 16): # 0,0,1 with 1,0,1
        vertex_interp(gv[4], gv[5], isovalue, vertlist[4],
                      dds, x, y, z, 4, 5)
    if (edge_table[cubeindex] & 32): # 1,0,1 with 1,1,1
        vertex_interp(gv[5], gv[6], isovalue, vertlist[5],
                      dds, x, y, z, 5, 6)
    if (edge_table[cubeindex] & 64): # 1,1,1 with 0,1,1
        vertex_interp(gv[6], gv[7], isovalue, vertlist[6],
                      dds, x, y, z, 6, 7)
    if (edge_table[cubeindex] & 128): # 0,1,1 with 0,0,1
        vertex_interp(gv[7], gv[4], isovalue, vertlist[7],
                      dds, x, y, z, 7, 4)
    if (edge_table[cubeindex] & 256): # 0,0,0 with 0,0,1
        vertex_interp(gv[0], gv[4], isovalue, vertlist[8],
                      dds, x, y, z, 0, 4)
    if (edge_table[cubeindex] & 512): # 1,0,0 with 1,0,1
        vertex_interp(gv[1], gv[5], isovalue, vertlist[9],
                      dds, x, y, z, 1, 5)
    if (edge_table[cubeindex] & 1024): # 1,1,0 with 1,1,1
        vertex_interp(gv[2], gv[6], isovalue, vertlist[10],
                      dds, x, y, z, 2, 6)
    if (edge_table[cubeindex] & 2048): # 0,1,0 with 0,1,1
        vertex_interp(gv[3], gv[7], isovalue, vertlist[11],
                      dds, x, y, z, 3, 7)
    n = 0
    while 1:
        triangles.current = AddTriangle(triangles.current,
                    vertlist[tri_table[cubeindex][n  ]],
                    vertlist[tri_table[cubeindex][n+1]],
                    vertlist[tri_table[cubeindex][n+2]])
        triangles.count += 1
        nt += 1
        if triangles.first == NULL:
            triangles.first = triangles.current
        n += 3
        if tri_table[cubeindex][n] == -1: break
    return nt

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def march_cubes_grid(np.float64_t isovalue,
                     np.ndarray[np.float64_t, ndim=3] values,
                     np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
                     np.ndarray[np.float64_t, ndim=1] left_edge,
                     np.ndarray[np.float64_t, ndim=1] dxs,
                     obj_sample = None, int sample_type = 1):
    cdef int dims[3]
    cdef int i, j, k, n, m, nt
    cdef int offset
    cdef np.float64_t gv[8]
    cdef np.float64_t pos[3]
    cdef np.float64_t point[3]
    cdef np.float64_t idds[3]
    cdef np.float64_t *intdata = NULL
    cdef np.float64_t *sdata = NULL
    cdef np.float64_t do_sample
    cdef np.ndarray[np.float64_t, ndim=3] sample
    cdef np.ndarray[np.float64_t, ndim=1] sampled
    cdef TriangleCollection triangles
    cdef Triangle *last
    cdef Triangle *current
    if obj_sample is not None:
        sample = obj_sample
        sdata = <np.float64_t *> sample.data
        do_sample = sample_type # 1 for face, 2 for vertex
    else:
        do_sample = 0
    for i in range(3):
        dims[i] = values.shape[i] - 1
        idds[i] = 1.0 / dxs[i]
    triangles.first = triangles.current = NULL
    last = current = NULL
    triangles.count = 0
    cdef np.float64_t *data = <np.float64_t *> values.data
    cdef np.float64_t *dds = <np.float64_t *> dxs.data
    pos[0] = left_edge[0]
    for i in range(dims[0]):
        pos[1] = left_edge[1]
        for j in range(dims[1]):
            pos[2] = left_edge[2]
            for k in range(dims[2]):
                if mask[i,j,k] == 1:
                    offset = i * (dims[1] + 1) * (dims[2] + 1) \
                           + j * (dims[2] + 1) + k
                    intdata = data + offset
                    offset_fill(dims, intdata, gv)
                    nt = march_cubes(gv, isovalue, dds, pos[0], pos[1], pos[2],
                                &triangles)
                    if nt == 0 or do_sample == 0:
                        pos[2] += dds[2]
                        continue
                    if last == NULL and triangles.first != NULL:
                        current = triangles.first
                        last = NULL
                    elif last != NULL:
                        current = last.next
                    if do_sample == 1:
                        # At each triangle's center, sample our secondary field
                        while current != NULL:
                            for n in range(3):
                                point[n] = 0.0
                            for n in range(3):
                                for m in range(3):
                                    point[m] += (current.p[n][m]-pos[m])*idds[m]
                            for n in range(3):
                                point[n] /= 3.0
                            current.val[0] = offset_interpolate(dims, point,
                                                             sdata + offset)
                            last = current
                            if current.next == NULL: break
                            current = current.next
                    elif do_sample == 2:
                        while current != NULL:
                            for n in range(3):
                                for m in range(3):
                                    point[m] = (current.p[n][m]-pos[m])*idds[m]
                                current.val[n] = offset_interpolate(dims,
                                                    point, sdata + offset)
                            last = current
                            if current.next == NULL: break
                            current = current.next
                pos[2] += dds[2]
            pos[1] += dds[1]
        pos[0] += dds[0]
    # Hallo, we are all done.
    cdef np.ndarray[np.float64_t, ndim=2] vertices
    vertices = np.zeros((triangles.count*3,3), dtype='float64')
    if do_sample == 0:
        FillAndWipeTriangles(vertices, triangles.first)
        return vertices
    cdef int nskip = 0
    if do_sample == 1:
        nskip = 1
    elif do_sample == 2:
        nskip = 3
    sampled = np.zeros(triangles.count * nskip, dtype='float64')
    FillTriangleValues(sampled, triangles.first, nskip)
    FillAndWipeTriangles(vertices, triangles.first)
    return vertices, sampled

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def march_cubes_grid_flux(
                     np.float64_t isovalue,
                     np.ndarray[np.float64_t, ndim=3] values,
                     np.ndarray[np.float64_t, ndim=3] v1,
                     np.ndarray[np.float64_t, ndim=3] v2,
                     np.ndarray[np.float64_t, ndim=3] v3,
                     np.ndarray[np.float64_t, ndim=3] flux_field,
                     np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
                     np.ndarray[np.float64_t, ndim=1] left_edge,
                     np.ndarray[np.float64_t, ndim=1] dxs):
    cdef int dims[3]
    cdef int i, j, k, n, m
    cdef int offset
    cdef np.float64_t gv[8]
    cdef np.float64_t *intdata = NULL
    cdef TriangleCollection triangles
    cdef Triangle *current = NULL
    cdef Triangle *last = NULL
    cdef np.float64_t *data = <np.float64_t *> values.data
    cdef np.float64_t *v1data = <np.float64_t *> v1.data
    cdef np.float64_t *v2data = <np.float64_t *> v2.data
    cdef np.float64_t *v3data = <np.float64_t *> v3.data
    cdef np.float64_t *fdata = <np.float64_t *> flux_field.data
    cdef np.float64_t *dds = <np.float64_t *> dxs.data
    cdef np.float64_t flux = 0.0
    cdef np.float64_t temp, area, s
    cdef np.float64_t center[3]
    cdef np.float64_t point[3]
    cdef np.float64_t cell_pos[3]
    cdef np.float64_t fv[3]
    cdef np.float64_t idds[3]
    cdef np.float64_t normal[3]
    for i in range(3):
        dims[i] = values.shape[i] - 1
        idds[i] = 1.0 / dds[i]
    triangles.first = triangles.current = NULL
    triangles.count = 0
    cell_pos[0] = left_edge[0]
    for i in range(dims[0]):
        cell_pos[1] = left_edge[1]
        for j in range(dims[1]):
            cell_pos[2] = left_edge[2]
            for k in range(dims[2]):
                if mask[i,j,k] == 1:
                    offset = i * (dims[1] + 1) * (dims[2] + 1) \
                           + j * (dims[2] + 1) + k
                    intdata = data + offset
                    offset_fill(dims, intdata, gv)
                    march_cubes(gv, isovalue, dds,
                                cell_pos[0], cell_pos[1], cell_pos[2],
                                &triangles)
                    # Now our triangles collection has a bunch.  We now
                    # calculate fluxes for each.
                    if last == NULL and triangles.first != NULL:
                        current = triangles.first
                        last = NULL
                    elif last != NULL:
                        current = last.next
                    while current != NULL:
                        # Calculate the center of the triangle
                        wval = 0.0
                        for n in range(3):
                            center[n] = 0.0
                        for n in range(3):
                            for m in range(3):
                                point[m] = (current.p[n][m]-cell_pos[m])*idds[m]
                            # Now we calculate the value at this point
                            temp = offset_interpolate(dims, point, intdata)
                            #print "something", temp, point[0], point[1], point[2]
                            wval += temp
                            for m in range(3):
                                center[m] += temp * point[m]
                        # Now we divide by our normalizing factor
                        for n in range(3):
                            center[n] /= wval
                        # We have our center point of the triangle, in 0..1
                        # coordinates.  So now we interpolate our three
                        # fields.
                        fv[0] = offset_interpolate(dims, center, v1data + offset)
                        fv[1] = offset_interpolate(dims, center, v2data + offset)
                        fv[2] = offset_interpolate(dims, center, v3data + offset)
                        # We interpolate again the actual value data
                        wval = offset_interpolate(dims, center, fdata + offset)
                        # Now we have our flux vector and our field value!
                        # We just need a normal vector with which we can
                        # dot it.  The normal should be equal to the gradient
                        # in the center of the triangle, or thereabouts.
                        eval_gradient(dims, center, intdata, normal)
                        temp = 0.0
                        for n in range(3):
                            temp += normal[n]*normal[n]
                        # Take the negative, to ensure it points inwardly
                        temp = -(temp**0.5)
                        # Dump this somewhere for now
                        temp = wval * (fv[0] * normal[0] +
                                       fv[1] * normal[1] +
                                       fv[2] * normal[2])/temp
                        # Now we need the area of the triangle.  This will take
                        # a lot of time to calculate compared to the rest.
                        # We use Heron's formula.
                        for n in range(3):
                            fv[n] = 0.0
                        for n in range(3):
                            fv[0] += (current.p[0][n] - current.p[2][n])**2.0
                            fv[1] += (current.p[1][n] - current.p[0][n])**2.0
                            fv[2] += (current.p[2][n] - current.p[1][n])**2.0
                        s = 0.0
                        for n in range(3):
                            fv[n] = fv[n]**0.5
                            s += 0.5 * fv[n]
                        area = (s*(s-fv[0])*(s-fv[1])*(s-fv[2]))
                        area = area**0.5
                        flux += temp*area
                        last = current
                        if current.next == NULL: break
                        current = current.next
                cell_pos[2] += dds[2]
            cell_pos[1] += dds[1]
        cell_pos[0] += dds[0]
    # Hallo, we are all done.
    WipeTriangles(triangles.first)
    return flux

