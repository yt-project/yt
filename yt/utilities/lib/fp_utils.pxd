"""
Shareable definitions for common fp/int Cython utilities



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

cdef inline int imax(int i0, int i1) nogil:
    if i0 > i1: return i0
    return i1

cdef inline np.float64_t fmax(np.float64_t f0, np.float64_t f1) nogil:
    if f0 > f1: return f0
    return f1

cdef inline int imin(int i0, int i1) nogil:
    if i0 < i1: return i0
    return i1

cdef inline np.float64_t fmin(np.float64_t f0, np.float64_t f1) nogil:
    if f0 < f1: return f0
    return f1

cdef inline np.float64_t fabs(np.float64_t f0) nogil:
    if f0 < 0.0: return -f0
    return f0

cdef inline int iclip(int i, int a, int b) nogil:
    if i < a: return a
    if i > b: return b
    return i

cdef inline int i64clip(np.int64_t i, np.int64_t a, np.int64_t b) nogil:
    if i < a: return a
    if i > b: return b
    return i

cdef inline np.float64_t fclip(np.float64_t f,
                      np.float64_t a, np.float64_t b) nogil:
    return fmin(fmax(f, a), b)

cdef inline np.int64_t i64max(np.int64_t i0, np.int64_t i1) nogil:
    if i0 > i1: return i0
    return i1

cdef inline np.int64_t i64min(np.int64_t i0, np.int64_t i1) nogil:
    if i0 < i1: return i0
    return i1

