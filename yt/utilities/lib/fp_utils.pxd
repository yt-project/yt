"""
Shareable definitions for common fp/int Cython utilities



"""


cimport cython
cimport numpy as np


cdef inline np.int64_t imax(np.int64_t i0, np.int64_t i1) nogil:
    if i0 > i1: return i0
    return i1

cdef inline np.float64_t fmax(np.float64_t f0, np.float64_t f1) nogil:
    if f0 > f1: return f0
    return f1

cdef inline np.int64_t imin(np.int64_t i0, np.int64_t i1) nogil:
    if i0 < i1: return i0
    return i1

cdef inline np.float64_t fmin(np.float64_t f0, np.float64_t f1) nogil:
    if f0 < f1: return f0
    return f1

cdef inline np.float64_t fabs(np.float64_t f0) nogil:
    if f0 < 0.0: return -f0
    return f0

cdef inline np.int64_t iclip(np.int64_t i, np.int64_t a, np.int64_t b) nogil:
    if i < a: return a
    if i > b: return b
    return i

cdef inline np.int64_t i64clip(np.int64_t i, np.int64_t a, np.int64_t b) nogil:
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

cdef inline _ensure_code(arr):
    if hasattr(arr, "units"):
        if "code_length" == str(arr.units):
            return arr
        arr.convert_to_units("code_length")
    return arr
