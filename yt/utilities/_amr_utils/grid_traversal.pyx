"""
Simple integrators for the radiative transfer equation

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

import numpy as np
cimport numpy as np
cimport cython
cimport kdtree_utils
cimport healpix_interface
from stdlib cimport malloc, free, abs
from fp_utils cimport imax, fmax, imin, fmin, iclip, fclip

from cython.parallel import prange, parallel, threadid

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

cdef struct VolumeContainer
ctypedef void sample_function(
                VolumeContainer *vc,
                np.float64_t v_pos[3],
                np.float64_t v_dir[3],
                np.float64_t enter_t,
                np.float64_t exit_t,
                int index[3],
                void *data) nogil

cdef extern from "FixedInterpolator.h":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data)
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3], np.float64_t *data)
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data)
    void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                       np.float64_t grad[3])
    void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval)
    void vertex_interp(np.float64_t v1, np.float64_t v2, np.float64_t isovalue,
                       np.float64_t vl[3], np.float64_t dds[3],
                       np.float64_t x, np.float64_t y, np.float64_t z,
                       int vind1, int vind2)

cdef struct VolumeContainer:
    int n_fields
    np.float64_t **data
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    np.float64_t dds[3]
    np.float64_t idds[3]
    int dims[3]

cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef VolumeContainer *container

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self,
                  int parent_grid_id, data,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.int64_t, ndim=1] dims):
        # The data is likely brought in via a slice, so we copy it
        cdef np.ndarray[np.float64_t, ndim=3] tdata
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        self.container = <VolumeContainer *> \
            malloc(sizeof(VolumeContainer))
        cdef VolumeContainer *c = self.container # convenience
        cdef int n_fields = len(data)
        c.n_fields = n_fields
        for i in range(3):
            c.left_edge[i] = left_edge[i]
            c.right_edge[i] = right_edge[i]
            c.dims[i] = dims[i]
            c.dds[i] = (c.right_edge[i] - c.left_edge[i])/dims[i]
            c.idds[i] = 1.0/c.dds[i]
        self.my_data = data
        c.data = <np.float64_t **> malloc(sizeof(np.float64_t*) * n_fields)
        for i in range(n_fields):
            tdata = data[i]
            c.data[i] = <np.float64_t *> tdata.data

    def __dealloc__(self):
        # The data fields are not owned by the container, they are owned by us!
        # So we don't need to deallocate them.
        free(self.container.data)
        free(self.container)

cdef struct ImageContainer:
    np.float64_t *vp_pos, *vp_dir, *center, *image,
    np.float64_t pdx, pdy, bounds[4]
    int nv[2]
    int vp_strides[3]
    int im_strides[3]
    int vd_strides[3]
    np.float64_t *x_vec, *y_vec

cdef struct ImageAccumulator:
    np.float64_t rgba[3]

cdef class ImageSampler:
    cdef ImageContainer *image
    cdef sample_function *sampler
    cdef public object avp_pos, avp_dir, acenter, aimage, ax_vec, ay_vec
    def __cinit__(self, 
                  np.ndarray[np.float64_t, ndim=3] vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec):
        self.image = <ImageContainer *> malloc(sizeof(ImageContainer))
        cdef ImageContainer *imagec = self.image
        self.sampler = NULL
        cdef int i, j
        # These assignments are so we can track the objects and prevent their
        # de-allocation from reference counts.
        self.avp_pos = vp_pos
        self.avp_dir = vp_dir
        self.acenter = center
        self.aimage = image
        self.ax_vec = x_vec
        self.ay_vec = y_vec
        imagec.vp_pos = <np.float64_t *> vp_pos.data
        imagec.vp_dir = <np.float64_t *> vp_dir.data
        imagec.center = <np.float64_t *> center.data
        imagec.image = <np.float64_t *> image.data
        imagec.x_vec = <np.float64_t *> x_vec.data
        imagec.y_vec = <np.float64_t *> y_vec.data
        imagec.nv[0] = vp_pos.shape[0]
        imagec.nv[1] = vp_pos.shape[1]
        for i in range(4): imagec.bounds[i] = bounds[i]
        imagec.pdx = (bounds[1] - bounds[0])/imagec.nv[0]
        imagec.pdy = (bounds[3] - bounds[2])/imagec.nv[1]
        for i in range(3):
            imagec.vp_strides[i] = vp_pos.strides[i] / 8
            imagec.im_strides[i] = image.strides[i] / 8
        if vp_dir.ndim > 1:
            for i in range(3):
                imagec.vd_strides[i] = vp_dir.strides[i] / 8
        else:
            imagec.vd_strides[0] = imagec.vd_strides[1] = imagec.vd_strides[2] = -1
        self.setup()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_start_stop(self, np.float64_t *ex, int *rv):
        # Extrema need to be re-centered
        cdef np.float64_t cx, cy
        cdef ImageContainer *im = self.image
        cdef int i
        cx = cy = 0.0
        for i in range(3):
            cx += im.center[i] * im.x_vec[i]
            cy += im.center[i] * im.y_vec[i]
        rv[0] = lrint((ex[0] - cx - im.bounds[0])/im.pdx)
        rv[1] = rv[0] + lrint((ex[1] - ex[0])/im.pdx)
        rv[2] = lrint((ex[2] - cy - im.bounds[2])/im.pdy)
        rv[3] = rv[2] + lrint((ex[3] - ex[2])/im.pdy)

    cdef inline void copy_into(self, np.float64_t *fv, np.float64_t *tv,
                        int i, int j, int nk, int strides[3]) nogil:
        # We know the first two dimensions of our from-vector, and our
        # to-vector is flat and 'ni' long
        cdef int k
        cdef int offset = strides[0] * i + strides[1] * j
        for k in range(nk):
            tv[k] = fv[offset + k]

    cdef inline void copy_back(self, np.float64_t *fv, np.float64_t *tv,
                        int i, int j, int nk, int strides[3]) nogil:
        cdef int k
        cdef int offset = strides[0] * i + strides[1] * j
        for k in range(nk):
            tv[offset + k] = fv[k]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void calculate_extent(self, np.float64_t extrema[4],
                               VolumeContainer *vc) nogil:
        # We do this for all eight corners
        cdef np.float64_t *edges[2], temp
        edges[0] = vc.left_edge
        edges[1] = vc.right_edge
        extrema[0] = extrema[2] = 1e300; extrema[1] = extrema[3] = -1e300
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    # This should rotate it into the vector plane
                    temp  = edges[i][0] * self.image.x_vec[0]
                    temp += edges[j][1] * self.image.x_vec[1]
                    temp += edges[k][2] * self.image.x_vec[2]
                    if temp < extrema[0]: extrema[0] = temp
                    if temp > extrema[1]: extrema[1] = temp
                    temp  = edges[i][0] * self.image.y_vec[0]
                    temp += edges[j][1] * self.image.y_vec[1]
                    temp += edges[k][2] * self.image.y_vec[2]
                    if temp < extrema[2]: extrema[2] = temp
                    if temp > extrema[3]: extrema[3] = temp

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self, PartitionedGrid pg):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, j, ni, nj, nn, offset
        cdef int iter[4]
        cdef VolumeContainer *vc = pg.container
        cdef ImageContainer *im = self.image
        if self.sampler == NULL: raise RuntimeError
        cdef np.float64_t *v_pos, v_dir[3], rgba[6], extrema[4]
        hit = 0
        self.calculate_extent(extrema, vc)
        self.get_start_stop(extrema, iter)
        iter[0] = iclip(iter[0]-1, 0, im.nv[0])
        if iter[0] == 0: return
        iter[1] = iclip(iter[1]+1, 0, im.nv[0])
        if iter[1] == 0 or iter[1] == iter[0]: return
        iter[2] = iclip(iter[2]-1, 0, im.nv[1])
        if iter[2] == 0: return
        iter[3] = iclip(iter[3]+1, 0, im.nv[1])
        if iter[3] == 0 or iter[3] == iter[2]: return
        cdef ImageAccumulator *idata
        cdef void *data
        cdef int nx = (iter[1] - iter[0])
        cdef int ny = (iter[3] - iter[2])
        cdef int size = nx * ny
        if im.vd_strides[0] == -1:
            with nogil, parallel():
                idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
                v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
                for j in prange(size, schedule="dynamic"):
                    vj = j % ny
                    vi = (j - vj) / ny + iter[0]
                    vj = vj + iter[2]
                    offset = im.vp_strides[0] * vi + im.vp_strides[1] * vj
                    for i in range(3): v_pos[i] = im.vp_pos[i + offset]
                    offset = im.im_strides[0] * vi + im.im_strides[1] * vj
                    for i in range(3): idata.rgba[i] = im.image[i + offset]
                    walk_volume(vc, v_pos, im.vp_dir, self.sampler,
                                (<void *> idata))
                    for i in range(3): im.image[i + offset] = idata.rgba[i]
                free(idata)
        else:
            # If we do not have an orthographic projection, we have to cast all
            # our rays (until we can get an extrema calculation...)
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            data = <void *> idata
            for vi in range(im.nv[0]):
                for vj in range(im.nv[1]):
                    for i in range(4): idata.rgba[i] = 0.0
                    self.copy_into(im.vp_pos, v_pos, vi, vj, 3, im.vp_strides)
                    self.copy_into(im.image, idata.rgba, vi, vj, 3, im.im_strides)
                    self.copy_into(im.vp_dir, v_dir, vi, vj, 3, im.vd_strides)
                    walk_volume(vc, v_pos, v_dir, self.sampler, data)
                    self.copy_back(idata.rgba, im.image, vi, vj, 3, im.im_strides)
        return hit

cdef void projection_sampler(
                 VolumeContainer *vc, 
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
    cdef ImageAccumulator *im = <ImageAccumulator *> data
    cdef int i
    cdef np.float64_t dl = (exit_t - enter_t)
    # We need this because by default it assumes vertex-centered data.
    for i in range(3):
        if index[i] < 0 or index[i] >= vc.dims[i]: return
    cdef int di = (index[0]*(vc.dims[1])+index[1])*vc.dims[2]+index[2]
    for i in range(imin(3, vc.n_fields)):
        im.rgba[i] += vc.data[i][di] * dl

cdef class ProjectionSampler(ImageSampler):
    def setup(self):
        self.sampler = projection_sampler

cdef class GridFace:
    cdef int direction
    cdef public np.float64_t coord
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, grid, int direction, int left):
        self.direction = direction
        if left == 1:
            self.coord = grid.LeftEdge[direction]
        else:
            self.coord = grid.RightEdge[direction]
        cdef int i
        for i in range(3):
            self.left_edge[i] = grid.LeftEdge[i]
            self.right_edge[i] = grid.RightEdge[i]
        self.left_edge[direction] = self.right_edge[direction] = self.coord

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int proj_overlap(self, np.float64_t *left_edge, np.float64_t *right_edge):
        cdef int xax, yax
        xax = (self.direction + 1) % 3
        yax = (self.direction + 2) % 3
        if left_edge[xax] >= self.right_edge[xax]: return 0
        if right_edge[xax] <= self.left_edge[xax]: return 0
        if left_edge[yax] >= self.right_edge[yax]: return 0
        if right_edge[yax] <= self.left_edge[yax]: return 0
        return 1

cdef class ProtoPrism:
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef public object subgrid_faces
    cdef public int parent_grid_id
    def __cinit__(self, int parent_grid_id,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  subgrid_faces):
        self.parent_grid_id = parent_grid_id
        cdef int i
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
        self.subgrid_faces = subgrid_faces

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def sweep(self, int direction = 0, int stack = 0):
        cdef int i
        cdef GridFace face
        cdef np.float64_t proto_split[3]
        for i in range(3): proto_split[i] = self.right_edge[i]
        for face in self.subgrid_faces[direction]:
            proto_split[direction] = face.coord
            if proto_split[direction] <= self.left_edge[direction]:
                continue
            if proto_split[direction] == self.right_edge[direction]:
                if stack == 2: return [self]
                return self.sweep((direction + 1) % 3, stack + 1)
            if face.proj_overlap(self.left_edge, proto_split) == 1:
                left, right = self.split(proto_split, direction)
                LC = left.sweep((direction + 1) % 3)
                RC = right.sweep(direction)
                return LC + RC
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef object split(self, np.float64_t *sp, int direction):
        cdef int i
        cdef np.ndarray split_left = self.LeftEdge.copy()
        cdef np.ndarray split_right = self.RightEdge.copy()

        for i in range(3): split_left[i] = self.right_edge[i]
        split_left[direction] = sp[direction]
        left = ProtoPrism(self.parent_grid_id, self.LeftEdge, split_left,
                          self.subgrid_faces)

        for i in range(3): split_right[i] = self.left_edge[i]
        split_right[direction] = sp[direction]
        right = ProtoPrism(self.parent_grid_id, split_right, self.RightEdge,
                           self.subgrid_faces)

        return (left, right)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_brick(self, np.ndarray[np.float64_t, ndim=1] grid_left_edge,
                        np.ndarray[np.float64_t, ndim=1] grid_dds,
                        child_mask):
        # We get passed in the left edge, the dds (which gives dimensions) and
        # the data, which is already vertex-centered.
        cdef PartitionedGrid PG
        cdef int li[3], ri[3], idims[3], i
        for i in range(3):
            li[i] = lrint((self.left_edge[i] - grid_left_edge[i])/grid_dds[i])
            ri[i] = lrint((self.right_edge[i] - grid_left_edge[i])/grid_dds[i])
            idims[i] = ri[i] - li[i]
        if child_mask[li[0], li[1], li[2]] == 0: return []
        cdef np.ndarray[np.int64_t, ndim=1] dims = np.empty(3, dtype='int64')
        for i in range(3):
            dims[i] = idims[i]
        #cdef np.ndarray[np.float64_t, ndim=3] new_data
        #new_data = data[li[0]:ri[0]+1,li[1]:ri[1]+1,li[2]:ri[2]+1].copy()
        #PG = PartitionedGrid(self.parent_grid_id, new_data,
        #                     self.LeftEdge, self.RightEdge, dims)
        return ((li[0], ri[0]), (li[1], ri[1]), (li[2], ri[2]), dims)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int walk_volume(VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     sample_function *sampler,
                     void *data,
                     np.float64_t *return_t = NULL,
                     np.float64_t enter_t = -1.0) nogil:
    cdef int cur_ind[3], step[3], x, y, i, n, flat_ind, hit, direction
    cdef np.float64_t intersect_t = 1.0
    cdef np.float64_t iv_dir[3]
    cdef np.float64_t intersect[3], tmax[3], tdelta[3]
    cdef np.float64_t dist, alpha, dt, exit_t
    cdef np.float64_t tr, tl, temp_x, temp_y, dv
    for i in range(3):
        if (v_dir[i] < 0):
            step[i] = -1
        elif (v_dir[i] == 0):
            step[i] = 1
            tmax[i] = 1e60
            iv_dir[i] = 1e60
            tdelta[i] = 1e-60
            continue
        else:
            step[i] = 1
        x = (i+1) % 3
        y = (i+2) % 3
        iv_dir[i] = 1.0/v_dir[i]
        tl = (vc.left_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tl*v_dir[x])
        temp_y = (v_pos[y] + tl*v_dir[y])
        if vc.left_edge[x] <= temp_x and temp_x <= vc.right_edge[x] and \
           vc.left_edge[y] <= temp_y and temp_y <= vc.right_edge[y] and \
           0.0 <= tl and tl < intersect_t:
            direction = i
            intersect_t = tl
        tr = (vc.right_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tr*v_dir[x])
        temp_y = (v_pos[y] + tr*v_dir[y])
        if vc.left_edge[x] <= temp_x and temp_x <= vc.right_edge[x] and \
           vc.left_edge[y] <= temp_y and temp_y <= vc.right_edge[y] and \
           0.0 <= tr and tr < intersect_t:
            direction = i
            intersect_t = tr
    if vc.left_edge[0] <= v_pos[0] and v_pos[0] <= vc.right_edge[0] and \
       vc.left_edge[1] <= v_pos[1] and v_pos[1] <= vc.right_edge[1] and \
       vc.left_edge[2] <= v_pos[2] and v_pos[2] <= vc.right_edge[2]:
        intersect_t = 0.0
    if enter_t >= 0.0: intersect_t = enter_t
    if not ((0.0 <= intersect_t) and (intersect_t < 1.0)): return 0
    for i in range(3):
        intersect[i] = v_pos[i] + intersect_t * v_dir[i]
        cur_ind[i] = <int> floor((intersect[i] +
                                  step[i]*1e-8*vc.dds[i] -
                                  vc.left_edge[i])*vc.idds[i])
        tmax[i] = (((cur_ind[i]+step[i])*vc.dds[i])+
                    vc.left_edge[i]-v_pos[i])*iv_dir[i]
        # This deals with the asymmetry in having our indices refer to the
        # left edge of a cell, but the right edge of the brick being one
        # extra zone out.
        if cur_ind[i] == vc.dims[i] and step[i] < 0:
            cur_ind[i] = vc.dims[i] - 1
        if cur_ind[i] < 0 or cur_ind[i] >= vc.dims[i]: return 0
        if step[i] > 0:
            tmax[i] = (((cur_ind[i]+1)*vc.dds[i])
                        +vc.left_edge[i]-v_pos[i])*iv_dir[i]
        if step[i] < 0:
            tmax[i] = (((cur_ind[i]+0)*vc.dds[i])
                        +vc.left_edge[i]-v_pos[i])*iv_dir[i]
        tdelta[i] = (vc.dds[i]*iv_dir[i])
        if tdelta[i] < 0: tdelta[i] *= -1
    # We have to jumpstart our calculation
    enter_t = intersect_t
    hit = 0
    while 1:
        # dims here is one less than the dimensions of the data,
        # but we are tracing on the grid, not on the data...
        if (not (0 <= cur_ind[0] < vc.dims[0])) or \
           (not (0 <= cur_ind[1] < vc.dims[1])) or \
           (not (0 <= cur_ind[2] < vc.dims[2])):
            break
        hit += 1
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                exit_t = fmin(tmax[0], 1.0)
                sampler(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
                cur_ind[0] += step[0]
                enter_t = tmax[0]
                tmax[0] += tdelta[0]
            else:
                exit_t = fmin(tmax[2], 1.0)
                sampler(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
                cur_ind[2] += step[2]
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
        else:
            if tmax[1] < tmax[2]:
                exit_t = fmin(tmax[1], 1.0)
                sampler(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
                cur_ind[1] += step[1]
                enter_t = tmax[1]
                tmax[1] += tdelta[1]
            else:
                exit_t = fmin(tmax[2], 1.0)
                sampler(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
                cur_ind[2] += step[2]
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
        if enter_t >= 1.0: break
    if return_t != NULL: return_t[0] = exit_t
    return hit
