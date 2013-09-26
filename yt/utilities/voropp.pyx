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

cdef extern from "voro++.hh" namespace "voro":
    cdef cppclass c_loop_all
    
    cdef cppclass voronoicell:
        double volume()

    cdef cppclass container:
        container(double xmin, double xmax, double ymin, double ymax,
                  double zmin, double zmax, int nx, int ny, int nz,
                  libcpp.bool xper, libcpp.bool yper, libcpp.bool zper, int alloc)
        void put(int n, double x, double y, double z)
        void store_cell_volumes(double *vols)
        int compute_cell(voronoicell c, c_loop_all vl)
        double sum_cell_volumes()
		
    cdef cppclass c_loop_all:
        c_loop_all(container &con)
        int inc()
        int start()

cdef class VoronoiVolume:
    cdef container *my_con
    cdef public int npart
    def __init__(self, xi, yi, zi, left_edge, right_edge):
        self.my_con = new container(left_edge[0], right_edge[0],
                                    left_edge[1], right_edge[1],
                                    left_edge[2], right_edge[2],
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
        #self.my_con.store_cell_volumes(vdouble)
        cdef c_loop_all *vl = new c_loop_all(deref(self.my_con))
        cdef voronoicell c
        if not vl.start(): return
        cdef int i = 0
        while 1:
            if self.my_con.compute_cell(c, deref(vl)):
                vol[i] = c.volume()
            if not vl.inc(): break
            i += 1
        del vl
        return vol
