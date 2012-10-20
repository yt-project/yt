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
from libc.stdlib cimport malloc, free, abs
from fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip
from field_interpolation_tables cimport \
    FieldInterpolationTable, FIT_initialize_table, FIT_eval_transfer,\
    FIT_eval_transfer_with_light
from fixed_interpolator cimport *

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
    double nearbyint(double x) nogil
    double fabs(double x) nogil
    double atan(double x) nogil
    double atan2(double y, double x) nogil
    double acos(double x) nogil
    double asin(double x) nogil
    double cos(double x) nogil
    double sin(double x) nogil
    double sqrt(double x) nogil

cdef struct VolumeContainer:
    int n_fields
    np.float64_t **data
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    np.float64_t dds[3]
    np.float64_t idds[3]
    int dims[3]

ctypedef void sample_function(
                VolumeContainer *vc,
                np.float64_t v_pos[3],
                np.float64_t v_dir[3],
                np.float64_t enter_t,
                np.float64_t exit_t,
                int index[3],
                void *data) nogil

cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef public int parent_grid_id
    cdef VolumeContainer *container
    cdef kdtree_utils.kdtree *star_list
    cdef np.float64_t star_er
    cdef np.float64_t star_sigma_num
    cdef np.float64_t star_coeff

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self,
                  int parent_grid_id, data,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.int64_t, ndim=1] dims,
		  star_kdtree_container star_tree = None):
        # The data is likely brought in via a slice, so we copy it
        cdef np.ndarray[np.float64_t, ndim=3] tdata
        self.container = NULL
        self.parent_grid_id = parent_grid_id
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
        if star_tree is None:
            self.star_list = NULL
        else:
            self.set_star_tree(star_tree)

    def set_star_tree(self, star_kdtree_container star_tree):
        self.star_list = star_tree.tree
        self.star_sigma_num = 2.0*star_tree.sigma**2.0
        self.star_er = 2.326 * star_tree.sigma
        self.star_coeff = star_tree.coeff

    def __dealloc__(self):
        # The data fields are not owned by the container, they are owned by us!
        # So we don't need to deallocate them.
        if self.container == NULL: return
        if self.container.data != NULL: free(self.container.data)
        free(self.container)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def integrate_streamline(self, pos, np.float64_t h, mag):
        cdef np.float64_t cmag[1]
        cdef np.float64_t k1[3], k2[3], k3[3], k4[3]
        cdef np.float64_t newpos[3], oldpos[3]
        for i in range(3):
            newpos[i] = oldpos[i] = pos[i]
        self.get_vector_field(newpos, k1, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k1[i]*h

        if not (self.LeftEdge[0] < newpos[0] and newpos[0] < self.RightEdge[0] and \
                self.LeftEdge[1] < newpos[1] and newpos[1] < self.RightEdge[1] and \
                self.LeftEdge[2] < newpos[2] and newpos[2] < self.RightEdge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k2, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k2[i]*h

        if not (self.LeftEdge[0] <= newpos[0] and newpos[0] <= self.RightEdge[0] and \
                self.LeftEdge[1] <= newpos[1] and newpos[1] <= self.RightEdge[1] and \
                self.LeftEdge[2] <= newpos[2] and newpos[2] <= self.RightEdge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k3, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + k3[i]*h

        if not (self.LeftEdge[0] <= newpos[0] and newpos[0] <= self.RightEdge[0] and \
                self.LeftEdge[1] <= newpos[1] and newpos[1] <= self.RightEdge[1] and \
                self.LeftEdge[2] <= newpos[2] and newpos[2] <= self.RightEdge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k4, cmag)

        for i in range(3):
            pos[i] = oldpos[i] + h*(k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0)

        if mag is not None:
            for i in range(3):
                newpos[i] = pos[i]
            self.get_vector_field(newpos, k4, cmag)
            mag[0] = cmag[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_vector_field(self, np.float64_t pos[3],
                               np.float64_t *vel, np.float64_t *vel_mag):
        cdef np.float64_t dp[3]
        cdef int ci[3]
        cdef VolumeContainer *c = self.container # convenience

        for i in range(3):
            ci[i] = (int)((pos[i]-self.LeftEdge[i])/c.dds[i])
            dp[i] = (pos[i] - ci[i]*c.dds[i] - self.LeftEdge[i])/c.dds[i]

        cdef int offset = ci[0] * (c.dims[1] + 1) * (c.dims[2] + 1) \
                          + ci[1] * (c.dims[2] + 1) + ci[2]

        vel_mag[0] = 0.0
        for i in range(3):
            vel[i] = offset_interpolate(c.dims, dp, c.data[i] + offset)
            vel_mag[0] += vel[i]*vel[i]
        vel_mag[0] = np.sqrt(vel_mag[0])
        if vel_mag[0] != 0.0:
            for i in range(3):
                vel[i] /= vel_mag[0]

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
    void *supp_data

cdef class ImageSampler:
    cdef ImageContainer *image
    cdef sample_function *sampler
    cdef public object avp_pos, avp_dir, acenter, aimage, ax_vec, ay_vec
    cdef void *supp_data
    cdef np.float64_t width[3]
    def __init__(self, 
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  *args, **kwargs):
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
        imagec.nv[0] = image.shape[0]
        imagec.nv[1] = image.shape[1]
        for i in range(4): imagec.bounds[i] = bounds[i]
        imagec.pdx = (bounds[1] - bounds[0])/imagec.nv[0]
        imagec.pdy = (bounds[3] - bounds[2])/imagec.nv[1]
        for i in range(3):
            imagec.vp_strides[i] = vp_pos.strides[i] / 8
            imagec.im_strides[i] = image.strides[i] / 8
            self.width[i] = width[i]
        if vp_dir.ndim > 1:
            for i in range(3):
                imagec.vd_strides[i] = vp_dir.strides[i] / 8
        elif vp_pos.ndim == 1:
            imagec.vd_strides[0] = imagec.vd_strides[1] = imagec.vd_strides[2] = -1
        else:
            raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_start_stop(self, np.float64_t *ex, np.int64_t *rv):
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
    def __call__(self, PartitionedGrid pg, int num_threads = 0):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, j, ni, nj, nn
        cdef np.int64_t offset, iter[4]
        cdef VolumeContainer *vc = pg.container
        cdef ImageContainer *im = self.image
        self.setup(pg)
        if self.sampler == NULL: raise RuntimeError
        cdef np.float64_t *v_pos, *v_dir, rgba[6], extrema[4]
        hit = 0
        cdef np.int64_t nx, ny, size
        if im.vd_strides[0] == -1:
            self.calculate_extent(extrema, vc)
            self.get_start_stop(extrema, iter)
            iter[0] = i64clip(iter[0]-1, 0, im.nv[0])
            iter[1] = i64clip(iter[1]+1, 0, im.nv[0])
            iter[2] = i64clip(iter[2]-1, 0, im.nv[1])
            iter[3] = i64clip(iter[3]+1, 0, im.nv[1])
            nx = (iter[1] - iter[0])
            ny = (iter[3] - iter[2])
            size = nx * ny
        else:
            nx = im.nv[0]
            ny = 1
            iter[0] = iter[1] = iter[2] = iter[3] = 0
            size = nx
        cdef ImageAccumulator *idata
        cdef np.float64_t px, py 
        cdef np.float64_t width[3] 
        for i in range(3):
            width[i] = self.width[i]
        if im.vd_strides[0] == -1:
            with nogil, parallel(num_threads = num_threads):
                idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
                idata.supp_data = self.supp_data
                v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
                for j in prange(size, schedule="static",chunksize=1):
                    vj = j % ny
                    vi = (j - vj) / ny + iter[0]
                    vj = vj + iter[2]
                    # Dynamically calculate the position
                    px = width[0] * (<np.float64_t>vi)/(<np.float64_t>im.nv[0]-1) - width[0]/2.0
                    py = width[1] * (<np.float64_t>vj)/(<np.float64_t>im.nv[1]-1) - width[1]/2.0
                    v_pos[0] = im.vp_pos[0]*px + im.vp_pos[3]*py + im.vp_pos[9]
                    v_pos[1] = im.vp_pos[1]*px + im.vp_pos[4]*py + im.vp_pos[10]
                    v_pos[2] = im.vp_pos[2]*px + im.vp_pos[5]*py + im.vp_pos[11]
                    offset = im.im_strides[0] * vi + im.im_strides[1] * vj
                    for i in range(3): idata.rgba[i] = im.image[i + offset]
                    walk_volume(vc, v_pos, im.vp_dir, self.sampler,
                                (<void *> idata))
                    for i in range(3): im.image[i + offset] = idata.rgba[i]
                free(idata)
                free(v_pos)
        else:
            with nogil, parallel():
                idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
                idata.supp_data = self.supp_data
                v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
                v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
                # If we do not have a simple image plane, we have to cast all
                # our rays 
                for j in prange(size, schedule="dynamic", chunksize=100):
                    offset = j * 3
                    for i in range(3): v_pos[i] = im.vp_pos[i + offset]
                    for i in range(3): v_dir[i] = im.vp_dir[i + offset]
                    for i in range(3): idata.rgba[i] = im.image[i + offset]
                    walk_volume(vc, v_pos, v_dir, self.sampler, 
                                (<void *> idata))
                    for i in range(3): im.image[i + offset] = idata.rgba[i]
                free(v_dir)
                free(idata)
                free(v_pos)
        return hit

    cdef void setup(self, PartitionedGrid pg):
        return

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
    cdef int di = (index[0]*vc.dims[1]+index[1])*vc.dims[2]+index[2] 
    for i in range(imin(3, vc.n_fields)):
        im.rgba[i] += vc.data[i][di] * dl

cdef struct VolumeRenderAccumulator:
    int n_fits
    int n_samples
    FieldInterpolationTable *fits
    int field_table_ids[6]
    np.float64_t star_coeff
    np.float64_t star_er
    np.float64_t star_sigma_num
    kdtree_utils.kdtree *star_list
    np.float64_t *light_dir
    np.float64_t *light_rgba
    int grey_opacity


cdef class ProjectionSampler(ImageSampler):
    cdef void setup(self, PartitionedGrid pg):
        self.sampler = projection_sampler

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void interpolated_projection_sampler(
                 VolumeContainer *vc, 
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
    cdef ImageAccumulator *im = <ImageAccumulator *> data
    cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
            im.supp_data
    # we assume this has vertex-centered data.
    cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                    + index[1] * (vc.dims[2] + 1) + index[2]
    cdef np.float64_t slopes[6], dp[3], ds[3]
    cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
    cdef np.float64_t dvs[6]
    for i in range(3):
        dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
        dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
        dp[i] *= vc.idds[i]
        ds[i] = v_dir[i] * vc.idds[i] * dt
    for i in range(vri.n_samples):
        for j in range(vc.n_fields):
            dvs[j] = offset_interpolate(vc.dims, dp,
                    vc.data[j] + offset)
        for j in range(imin(3, vc.n_fields)):
            im.rgba[j] += dvs[j] * dt
        for j in range(3):
            dp[j] += ds[j]

cdef class InterpolatedProjectionSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
    def __cinit__(self, 
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  n_samples = 10):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width)
        cdef int i
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.n_samples = n_samples
        self.supp_data = <void *> self.vra

    cdef void setup(self, PartitionedGrid pg):
        self.sampler = interpolated_projection_sampler

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void volume_render_sampler(
                 VolumeContainer *vc, 
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
    cdef ImageAccumulator *im = <ImageAccumulator *> data
    cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
            im.supp_data
    # we assume this has vertex-centered data.
    cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                    + index[1] * (vc.dims[2] + 1) + index[2]
    cdef np.float64_t slopes[6], dp[3], ds[3]
    cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
    cdef np.float64_t dvs[6]
    for i in range(3):
        dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
        dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
        dp[i] *= vc.idds[i]
        ds[i] = v_dir[i] * vc.idds[i] * dt
    for i in range(vri.n_samples):
        for j in range(vc.n_fields):
            dvs[j] = offset_interpolate(vc.dims, dp,
                    vc.data[j] + offset)
        FIT_eval_transfer(dt, dvs, im.rgba, vri.n_fits, 
                vri.fits, vri.field_table_ids, vri.grey_opacity)
        for j in range(3):
            dp[j] += ds[j]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void volume_render_gradient_sampler(
                 VolumeContainer *vc, 
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
    cdef ImageAccumulator *im = <ImageAccumulator *> data
    cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
            im.supp_data
    # we assume this has vertex-centered data.
    cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                    + index[1] * (vc.dims[2] + 1) + index[2]
    cdef np.float64_t slopes[6], dp[3], ds[3]
    cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
    cdef np.float64_t dvs[6]
    cdef np.float64_t *grad
    grad = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
    for i in range(3):
        dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
        dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
        dp[i] *= vc.idds[i]
        ds[i] = v_dir[i] * vc.idds[i] * dt
    for i in range(vri.n_samples):
        for j in range(vc.n_fields):
            dvs[j] = offset_interpolate(vc.dims, dp,
                    vc.data[j] + offset)
        eval_gradient(vc.dims, dp, vc.data[0] + offset, grad)
        FIT_eval_transfer_with_light(dt, dvs, grad, 
                vri.light_dir, vri.light_rgba,
                im.rgba, vri.n_fits, 
                vri.fits, vri.field_table_ids, vri.grey_opacity)
        for j in range(3):
            dp[j] += ds[j]
    free(grad)

cdef class star_kdtree_container:
    cdef kdtree_utils.kdtree *tree
    cdef public np.float64_t sigma
    cdef public np.float64_t coeff

    def __init__(self):
        self.tree = kdtree_utils.kd_create(3)

    def add_points(self,
                   np.ndarray[np.float64_t, ndim=1] pos_x,
                   np.ndarray[np.float64_t, ndim=1] pos_y,
                   np.ndarray[np.float64_t, ndim=1] pos_z,
                   np.ndarray[np.float64_t, ndim=2] star_colors):
        cdef int i, n
        cdef np.float64_t *pointer = <np.float64_t *> star_colors.data
        for i in range(pos_x.shape[0]):
            kdtree_utils.kd_insert3(self.tree,
                pos_x[i], pos_y[i], pos_z[i], <void *> (pointer + i*3))

    def __dealloc__(self):
        kdtree_utils.kd_free(self.tree)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void volume_render_stars_sampler(
                 VolumeContainer *vc, 
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
    cdef ImageAccumulator *im = <ImageAccumulator *> data
    cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
            im.supp_data
    cdef kdtree_utils.kdres *ballq = NULL
    # we assume this has vertex-centered data.
    cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                    + index[1] * (vc.dims[2] + 1) + index[2]
    cdef np.float64_t slopes[6], dp[3], ds[3]
    cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
    cdef np.float64_t dvs[6], cell_left[3], local_dds[3], pos[3]
    cdef int nstars, dti, i, j
    cdef np.float64_t *colors = NULL, gexp, gaussian, px, py, pz
    for i in range(3):
        dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
        dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
        dp[i] *= vc.idds[i]
        ds[i] = v_dir[i] * vc.idds[i] * dt
    for i in range(vc.n_fields):
        slopes[i] = offset_interpolate(vc.dims, dp,
                        vc.data[i] + offset)
    cdef np.float64_t temp
    # Now we get the ball-tree result for the stars near our cell center.
    for i in range(3):
        cell_left[i] = index[i] * vc.dds[i] + vc.left_edge[i]
        pos[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
        local_dds[i] = v_dir[i] * dt
    ballq = kdtree_utils.kd_nearest_range3(
        vri.star_list, cell_left[0] + vc.dds[0]*0.5,
                        cell_left[1] + vc.dds[1]*0.5,
                        cell_left[2] + vc.dds[2]*0.5,
                        vri.star_er + 0.9*vc.dds[0])
                                    # ~0.866 + a bit

    nstars = kdtree_utils.kd_res_size(ballq)
    for i in range(vc.n_fields):
        temp = slopes[i]
        slopes[i] -= offset_interpolate(vc.dims, dp,
                         vc.data[i] + offset)
        slopes[i] *= -1.0/vri.n_samples
        dvs[i] = temp
    for dti in range(vri.n_samples): 
        # Now we add the contribution from stars
        kdtree_utils.kd_res_rewind(ballq)
        for i in range(nstars):
            kdtree_utils.kd_res_item3(ballq, &px, &py, &pz)
            colors = <np.float64_t *> kdtree_utils.kd_res_item_data(ballq)
            kdtree_utils.kd_res_next(ballq)
            gexp = (px - pos[0])*(px - pos[0]) \
                 + (py - pos[1])*(py - pos[1]) \
                 + (pz - pos[2])*(pz - pos[2])
            gaussian = vri.star_coeff * exp(-gexp/vri.star_sigma_num)
            for j in range(3): im.rgba[j] += gaussian*dt*colors[j]
        for i in range(3):
            pos[i] += local_dds[i]
        FIT_eval_transfer(dt, dvs, im.rgba, vri.n_fits, vri.fits,
                          vri.field_table_ids, vri.grey_opacity)
        for i in range(vc.n_fields):
            dvs[i] += slopes[i]
    kdtree_utils.kd_res_free(ballq)

cdef class VolumeRenderSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
    cdef kdtree_utils.kdtree **trees
    cdef object tree_containers
    def __cinit__(self, 
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  tf_obj, n_samples = 10,
                  star_list = None):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width)
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] temp
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.fits = <FieldInterpolationTable *> \
            malloc(sizeof(FieldInterpolationTable) * 6)
        self.vra.n_fits = tf_obj.n_field_tables
        assert(self.vra.n_fits <= 6)
        self.vra.grey_opacity = getattr(tf_obj, "grey_opacity", 0)
        self.vra.n_samples = n_samples
        self.my_field_tables = []
        for i in range(self.vra.n_fits):
            temp = tf_obj.tables[i].y
            FIT_initialize_table(&self.vra.fits[i],
                      temp.shape[0],
                      <np.float64_t *> temp.data,
                      tf_obj.tables[i].x_bounds[0],
                      tf_obj.tables[i].x_bounds[1],
                      tf_obj.field_ids[i], tf_obj.weight_field_ids[i],
                      tf_obj.weight_table_ids[i])
            self.my_field_tables.append((tf_obj.tables[i],
                                         tf_obj.tables[i].y))
        for i in range(6):
            self.vra.field_table_ids[i] = tf_obj.field_table_ids[i]
        self.supp_data = <void *> self.vra
        cdef star_kdtree_container skdc
        self.tree_containers = star_list
        if star_list is None:
            self.trees = NULL
        else:
            self.trees = <kdtree_utils.kdtree **> malloc(
                sizeof(kdtree_utils.kdtree*) * len(star_list))
            for i in range(len(star_list)):
                skdc = star_list[i]
                self.trees[i] = skdc.tree

    cdef void setup(self, PartitionedGrid pg):
        cdef star_kdtree_container star_tree
        if self.trees == NULL:
            self.sampler = volume_render_sampler
        else:
            star_tree = self.tree_containers[pg.parent_grid_id]
            self.vra.star_list = self.trees[pg.parent_grid_id]
            self.vra.star_sigma_num = 2.0*star_tree.sigma**2.0
            self.vra.star_er = 2.326 * star_tree.sigma
            self.vra.star_coeff = star_tree.coeff
            self.sampler = volume_render_stars_sampler

    def __dealloc__(self):
        return
        #free(self.vra.fits)
        #free(self.vra)

cdef class LightSourceRenderSampler(ImageSampler):
    cdef VolumeRenderAccumulator *vra
    cdef public object tf_obj
    cdef public object my_field_tables
    def __cinit__(self, 
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  tf_obj, n_samples = 10,
                  light_dir=[1.,1.,1.],
                  light_rgba=[1.,1.,1.,1.]):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width)
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] temp
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.fits = <FieldInterpolationTable *> \
            malloc(sizeof(FieldInterpolationTable) * 6)
        self.vra.n_fits = tf_obj.n_field_tables
        assert(self.vra.n_fits <= 6)
        self.vra.grey_opacity = getattr(tf_obj, "grey_opacity", 0)
        self.vra.n_samples = n_samples
        self.vra.light_dir = <np.float64_t *> malloc(sizeof(np.float64_t) * 3)
        self.vra.light_rgba = <np.float64_t *> malloc(sizeof(np.float64_t) * 4)
        light_dir /= np.sqrt(light_dir[0]**2 + light_dir[1]**2 + light_dir[2]**2)
        for i in range(3):
            self.vra.light_dir[i] = light_dir[i]
        for i in range(4):
            self.vra.light_rgba[i] = light_rgba[i]
        self.my_field_tables = []
        for i in range(self.vra.n_fits):
            temp = tf_obj.tables[i].y
            FIT_initialize_table(&self.vra.fits[i],
                      temp.shape[0],
                      <np.float64_t *> temp.data,
                      tf_obj.tables[i].x_bounds[0],
                      tf_obj.tables[i].x_bounds[1],
                      tf_obj.field_ids[i], tf_obj.weight_field_ids[i],
                      tf_obj.weight_table_ids[i])
            self.my_field_tables.append((tf_obj.tables[i],
                                         tf_obj.tables[i].y))
        for i in range(6):
            self.vra.field_table_ids[i] = tf_obj.field_table_ids[i]
        self.supp_data = <void *> self.vra

    cdef void setup(self, PartitionedGrid pg):
        self.sampler = volume_render_gradient_sampler

    def __dealloc__(self):
        return
        #free(self.vra.fits)
        #free(self.vra)
        #free(self.light_dir)
        #free(self.light_rgba)


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
                        np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask):
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
    cdef np.float64_t intersect_t = 1.1
    cdef np.float64_t iv_dir[3]
    cdef np.float64_t tmax[3], tdelta[3]
    cdef np.float64_t dist, alpha, dt, exit_t
    cdef np.float64_t tr, tl, temp_x, temp_y, dv
    direction = -1
    if vc.left_edge[0] <= v_pos[0] and v_pos[0] <= vc.right_edge[0] and \
       vc.left_edge[1] <= v_pos[1] and v_pos[1] <= vc.right_edge[1] and \
       vc.left_edge[2] <= v_pos[2] and v_pos[2] <= vc.right_edge[2]:
        intersect_t = 0.0
        direction = 3
    for i in range(3):
        if (v_dir[i] < 0):
            step[i] = -1
        elif (v_dir[i] == 0.0):
            step[i] = 0
            continue
        else:
            step[i] = 1
        iv_dir[i] = 1.0/v_dir[i]
        if direction == 3: continue
        x = (i+1) % 3
        y = (i+2) % 3
        if step[i] > 0:
            tl = (vc.left_edge[i] - v_pos[i])*iv_dir[i]
        else:
            tl = (vc.right_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tl*v_dir[x])
        temp_y = (v_pos[y] + tl*v_dir[y])
        if vc.left_edge[x] <= temp_x and temp_x <= vc.right_edge[x] and \
           vc.left_edge[y] <= temp_y and temp_y <= vc.right_edge[y] and \
           0.0 <= tl and tl < intersect_t:
            direction = i
            intersect_t = tl
    if enter_t >= 0.0: intersect_t = enter_t 
    if not ((0.0 <= intersect_t) and (intersect_t < 1.0)): return 0
    for i in range(3):
        # Two things have to be set inside this loop.
        # cur_ind[i], the current index of the grid cell the ray is in
        # tmax[i], the 't' until it crosses out of the grid cell
        tdelta[i] = step[i] * iv_dir[i] * vc.dds[i] 
        if i == direction and step[i] > 0:
            # Intersection with the left face in this direction
            cur_ind[i] = 0 
        elif i == direction and step[i] < 0:
            # Intersection with the right face in this direction
            cur_ind[i] = vc.dims[i] - 1
        else:
            # We are somewhere in the middle of the face
            temp_x = intersect_t * v_dir[i] + v_pos[i] # current position
            temp_y = ((temp_x - vc.left_edge[i])*vc.idds[i])
            cur_ind[i] =  <int> (floor(temp_y))
        if step[i] > 0:
            temp_y = (cur_ind[i] + 1) * vc.dds[i] + vc.left_edge[i]
        elif step[i] < 0:
            temp_y = cur_ind[i] * vc.dds[i] + vc.left_edge[i]
        tmax[i] = (temp_y - v_pos[i]) * iv_dir[i]
        if step[i] == 0:
            tmax[i] = 1e60
    # We have to jumpstart our calculation
    for i in range(3):
        if cur_ind[i] == vc.dims[i] and step[i] == 1:
            return 0
        if cur_ind[i] == -1 and step[i] == -1:
            return 0
    enter_t = intersect_t
    hit = 0
    while 1:
        hit += 1
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                i = 0
            else:
                i = 2
        else:
            if tmax[1] < tmax[2]:
                i = 1
            else:
                i = 2
        exit_t = fmin(tmax[i], 1.0)
        sampler(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
        cur_ind[i] += step[i]
        enter_t = tmax[i]
        tmax[i] += tdelta[i]
        if cur_ind[i] < 0 or cur_ind[i] >= vc.dims[i] or enter_t >= 1.0:
            break
    if return_t != NULL: return_t[0] = exit_t
    return hit

def hp_pix2vec_nest(long nside, long ipix):
    cdef double v[3]
    healpix_interface.pix2vec_nest(nside, ipix, v)
    cdef np.ndarray[np.float64_t, ndim=1] tr = np.empty((3,), dtype='float64')
    tr[0] = v[0]
    tr[1] = v[1]
    tr[2] = v[2]
    return tr

def arr_pix2vec_nest(long nside,
                     np.ndarray[np.int64_t, ndim=1] aipix):
    cdef int n = aipix.shape[0]
    cdef int i
    cdef double v[3]
    cdef long ipix
    cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 3), dtype='float64')
    for i in range(n):
        ipix = aipix[i]
        healpix_interface.pix2vec_nest(nside, ipix, v)
        tr[i,0] = v[0]
        tr[i,1] = v[1]
        tr[i,2] = v[2]
    return tr

def hp_vec2pix_nest(long nside, double x, double y, double z):
    cdef double v[3]
    v[0] = x
    v[1] = y
    v[2] = z
    cdef long ipix
    healpix_interface.vec2pix_nest(nside, v, &ipix)
    return ipix

def arr_vec2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] x,
                     np.ndarray[np.float64_t, ndim=1] y,
                     np.ndarray[np.float64_t, ndim=1] z):
    cdef int n = x.shape[0]
    cdef int i
    cdef double v[3]
    cdef long ipix
    cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    for i in range(n):
        v[0] = x[i]
        v[1] = y[i]
        v[2] = z[i]
        healpix_interface.vec2pix_nest(nside, v, &ipix)
        tr[i] = ipix
    return tr

def hp_pix2ang_nest(long nside, long ipnest):
    cdef double theta, phi
    healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
    return (theta, phi)

def arr_pix2ang_nest(long nside, np.ndarray[np.int64_t, ndim=1] aipnest):
    cdef int n = aipnest.shape[0]
    cdef int i
    cdef long ipnest
    cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 2), dtype='float64')
    cdef double theta, phi
    for i in range(n):
        ipnest = aipnest[i]
        healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
        tr[i,0] = theta
        tr[i,1] = phi
    return tr

def hp_ang2pix_nest(long nside, double theta, double phi):
    cdef long ipix
    healpix_interface.ang2pix_nest(nside, theta, phi, &ipix)
    return ipix

def arr_ang2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] atheta,
                     np.ndarray[np.float64_t, ndim=1] aphi):
    cdef int n = atheta.shape[0]
    cdef int i
    cdef long ipnest
    cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    cdef double theta, phi
    for i in range(n):
        theta = atheta[i]
        phi = aphi[i]
        healpix_interface.ang2pix_nest(nside, theta, phi, &ipnest)
        tr[i] = ipnest
    return tr

@cython.boundscheck(False)
@cython.cdivision(False)
@cython.wraparound(False)
def pixelize_healpix(long nside,
                     np.ndarray[np.float64_t, ndim=1] values,
                     long ntheta, long nphi,
                     np.ndarray[np.float64_t, ndim=2] irotation):
    # We will first to pix2vec, rotate, then calculate the angle
    cdef int i, j, thetai, phii
    cdef long ipix
    cdef double v0[3], v1[3]
    cdef double pi = 3.1415926
    cdef np.float64_t pi2 = pi/2.0
    cdef np.float64_t phi, theta
    cdef np.ndarray[np.float64_t, ndim=2] results
    cdef np.ndarray[np.int32_t, ndim=2] count
    results = np.zeros((ntheta, nphi), dtype="float64")
    count = np.zeros((ntheta, nphi), dtype="int32")

    cdef np.float64_t phi0 = 0
    cdef np.float64_t dphi = 2.0 * pi/(nphi-1)

    cdef np.float64_t theta0 = 0
    cdef np.float64_t dtheta = pi/(ntheta-1)
    # We assume these are the rotated theta and phi
    for thetai in range(ntheta):
        theta = theta0 + dtheta * thetai
        for phii in range(nphi):
            phi = phi0 + dphi * phii
            # We have our rotated vector
            v1[0] = cos(phi) * sin(theta)
            v1[1] = sin(phi) * sin(theta)
            v1[2] = cos(theta)
            # Now we rotate back
            for i in range(3):
                v0[i] = 0
                for j in range(3):
                    v0[i] += v1[j] * irotation[j,i]
            # Get the pixel this vector is inside
            healpix_interface.vec2pix_nest(nside, v0, &ipix)
            results[thetai, phii] = values[ipix]
            count[i, j] += 1
    return results, count
    #for i in range(ntheta):
    #    for j in range(nphi):
    #        if count[i,j] > 0:
    #            results[i,j] /= count[i,j]
    #return results, count

def healpix_aitoff_proj(np.ndarray[np.float64_t, ndim=1] pix_image,
                        long nside,
                        np.ndarray[np.float64_t, ndim=2] image,
                        np.ndarray[np.float64_t, ndim=2] irotation):
    cdef double pi = np.pi
    cdef int i, j, k, l
    cdef np.float64_t x, y, z, zb
    cdef np.float64_t dx, dy, inside
    cdef double v0[3], v1[3]
    dx = 2.0 / (image.shape[1] - 1)
    dy = 2.0 / (image.shape[0] - 1)
    cdef np.float64_t s2 = sqrt(2.0)
    cdef long ipix
    for i in range(image.shape[1]):
        x = (-1.0 + i*dx)*s2*2.0
        for j in range(image.shape[0]):
            y = (-1.0 + j * dy)*s2
            zb = (x*x/8.0 + y*y/2.0 - 1.0)
            if zb > 0: continue
            z = (1.0 - (x/4.0)**2.0 - (y/2.0)**2.0)
            z = z**0.5
            # Longitude
            phi = (2.0*atan(z*x/(2.0 * (2.0*z*z-1.0))) + pi)
            # Latitude
            # We shift it into co-latitude
            theta = (asin(z*y) + pi/2.0)
            # Now to account for rotation we translate into vectors
            v1[0] = cos(phi) * sin(theta)
            v1[1] = sin(phi) * sin(theta)
            v1[2] = cos(theta)
            for k in range(3):
                v0[k] = 0
                for l in range(3):
                    v0[k] += v1[l] * irotation[l,k]
            healpix_interface.vec2pix_nest(nside, v0, &ipix)
            #print "Rotated", v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], ipix, pix_image[ipix]
            image[j, i] = pix_image[ipix]

def arr_fisheye_vectors(int resolution, np.float64_t fov, int nimx=1, int
        nimy=1, int nimi=0, int nimj=0, np.float64_t off_theta=0.0, np.float64_t
        off_phi=0.0):
    # We now follow figures 4-7 of:
    # http://paulbourke.net/miscellaneous/domefisheye/fisheye/
    # ...but all in Cython.
    cdef np.ndarray[np.float64_t, ndim=3] vp
    cdef int i, j, k
    cdef np.float64_t r, phi, theta, px, py
    cdef np.float64_t pi = 3.1415926
    cdef np.float64_t fov_rad = fov * pi / 180.0
    cdef int nx = resolution/nimx
    cdef int ny = resolution/nimy
    vp = np.zeros((nx,ny, 3), dtype="float64")
    for i in range(nx):
        px = 2.0 * (nimi*nx + i) / (resolution) - 1.0
        for j in range(ny):
            py = 2.0 * (nimj*ny + j) / (resolution) - 1.0
            r = (px*px + py*py)**0.5
            if r == 0.0:
                phi = 0.0
            elif px < 0:
                phi = pi - asin(py / r)
            else:
                phi = asin(py / r)
            theta = r * fov_rad / 2.0
            theta += off_theta
            phi += off_phi
            vp[i,j,0] = sin(theta) * cos(phi)
            vp[i,j,1] = sin(theta) * sin(phi)
            vp[i,j,2] = cos(theta)
    return vp


