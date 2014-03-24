"""
This calls the ORIGAMI routines



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
from libc.stdlib cimport malloc, free

cdef extern from "origami_tags.h":
    int compute_tags(int ng, double boxsize, double **r, int npart,
                     unsigned char *m)

cdef int printed_citation = 0

def run_origami(np.ndarray[np.float64_t, ndim=1] pos_x,
                np.ndarray[np.float64_t, ndim=1] pos_y,
                np.ndarray[np.float64_t, ndim=1] pos_z,
                double boxsize):
    # We assume these have been passed in in the correct order and
    # C-contiguous.
    global printed_citation
    if printed_citation == 0:
        print "ORIGAMI was developed by Bridget Falck and Mark Neyrinck."
        print "Please cite Falck, Neyrinck, & Szalay 2012, ApJ, 754, 2, 125."
        printed_citation = 1
    cdef int npart = pos_x.size
    if npart == 1:
        return np.zeros(1, dtype="uint8")
    assert(sizeof(unsigned char) == sizeof(np.uint8_t))
    assert(sizeof(double) == sizeof(np.float64_t))
    cdef int ng = np.round(npart**(1./3))
    assert(ng**3 == npart)
    cdef double **r = <double **> malloc(sizeof(double *) * 3)
    r[0] = <double *> pos_x.data
    r[1] = <double *> pos_y.data
    r[2] = <double *> pos_z.data
    cdef np.ndarray[np.uint8_t, ndim=1] tags = np.zeros(npart, dtype="uint8")
    cdef void *m = <void*> tags.data
    compute_tags(ng, boxsize, r, npart, <unsigned char*> m)
    free(r)
    return tags
