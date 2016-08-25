"""
Image sampler definitions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, calloc, free, abs
from .fixed_interpolator cimport offset_interpolate

cdef class PartitionedGrid:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self,
                  int parent_grid_id, data,
                  mask,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.int64_t, ndim=1] dims):
        # The data is likely brought in via a slice, so we copy it
        cdef np.ndarray[np.float64_t, ndim=3] tdata
        cdef np.ndarray[np.uint8_t, ndim=3] mask_data
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
        self.source_mask = mask
        mask_data = mask
        c.data = <np.float64_t **> malloc(sizeof(np.float64_t*) * n_fields)
        for i in range(n_fields):
            tdata = data[i]
            c.data[i] = <np.float64_t *> tdata.data
        c.mask = <np.uint8_t *> mask_data.data

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
        cdef np.float64_t k1[3]
        cdef np.float64_t k2[3]
        cdef np.float64_t k3[3]
        cdef np.float64_t k4[3]
        cdef np.float64_t newpos[3]
        cdef np.float64_t oldpos[3]
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

