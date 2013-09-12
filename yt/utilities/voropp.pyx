"""
Wrapping code for voro++



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from cython.operator cimport dereference as deref, preincrement as inc
from libc.stdlib cimport malloc, free, abs, calloc, labs
cimport libcpp

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "voro++.cc":
    cdef cppclass container:
        container(double xmin, double xmax, double ymin, double ymax,
                  double zmin, double zmax, int nx, int ny, int nz,
                  libcpp.bool xper, libcpp.bool yper, libcpp.bool zper, int alloc)
        void put(int n, double x, double y, double z)
        void store_cell_volumes(double *vols)

cdef class VoronoiVolume:
    cdef container *my_con
    cdef int npart
    def __init__(self, xi, yi, zi):
        self.my_con = new container(0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                                    xi, yi, zi, False, False, False, 8)
        self.npart = 0

    def __dealloc__(self):
        del self.my_con

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_array(self, np.ndarray[np.float64_t, ndim=1] xpos,
                        np.ndarray[np.float64_t, ndim=1] ypos,
                        np.ndarray[np.float64_t, ndim=1] zpos):
        cdef int i
        for i in range(xpos.shape[0]):
            self.my_con.put(self.npart, xpos[i], ypos[i], zpos[i])
            self.npart += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_volumes(self):
        cdef np.ndarray vol = np.zeros(self.npart, 'double')
        cdef double *vdouble = <double *> vol.data
        self.my_con.store_cell_volumes(vdouble)
        return vol
