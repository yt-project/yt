"""
A light interface to a few HEALPix routines



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

from libc.stdio cimport fopen, fclose, FILE

cdef extern from "healpix_vectors.h":
    int pix2vec_nest(long nside, long ipix, double *v)
    void vec2pix_nest(long nside, double *vec, long *ipix)
    void pix2ang_nest(long nside, long ipix, double *theta, double *phi)
    void ang2pix_nest(long nside, double theta, double phi, long *ipix)
