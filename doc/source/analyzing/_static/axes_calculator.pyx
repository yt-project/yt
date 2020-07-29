import numpy as np

cimport cython
cimport numpy as np
from libc.stdlib cimport free, malloc


cdef extern from "axes.h":
    ctypedef struct ParticleCollection:
            long npart
            double *xpos
            double *ypos
            double *zpos

    void calculate_axes(ParticleCollection *part,
             double *ax1, double *ax2, double *ax3)

def examine_axes(np.ndarray[np.float64_t, ndim=1] xpos,
                 np.ndarray[np.float64_t, ndim=1] ypos,
                 np.ndarray[np.float64_t, ndim=1] zpos):
    cdef double ax1[3]
    cdef double ax2[3]
    cdef double ax3[3]
    cdef ParticleCollection particles
    cdef int i

    particles.npart = len(xpos)
    particles.xpos = <double *> xpos.data
    particles.ypos = <double *> ypos.data
    particles.zpos = <double *> zpos.data

    for i in range(particles.npart):
        particles.xpos[i] = xpos[i]
        particles.ypos[i] = ypos[i]
        particles.zpos[i] = zpos[i]

    calculate_axes(&particles, ax1, ax2, ax3)

    return ( (ax1[0], ax1[1], ax1[2]),
             (ax2[0], ax2[1], ax2[2]),
             (ax3[0], ax3[1], ax3[2]) )
