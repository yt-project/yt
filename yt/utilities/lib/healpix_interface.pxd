"""
A light interface to a few HEALPix routines



"""


import numpy as np

cimport cython
cimport numpy as np
from libc.stdio cimport FILE, fclose, fopen


cdef extern from "healpix_vectors.h":
    int pix2vec_nest(long nside, long ipix, double *v)
    void vec2pix_nest(long nside, double *vec, long *ipix)
    void pix2ang_nest(long nside, long ipix, double *theta, double *phi)
    void ang2pix_nest(long nside, double theta, double phi, long *ipix)
