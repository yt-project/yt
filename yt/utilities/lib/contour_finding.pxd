"""
Contour finding exports



"""



cimport cython
cimport numpy as np


cdef inline np.int64_t i64max(np.int64_t i0, np.int64_t i1):
    if i0 > i1: return i0
    return i1

cdef inline np.int64_t i64min(np.int64_t i0, np.int64_t i1):
    if i0 < i1: return i0
    return i1

cdef extern from "math.h":
    double fabs(double x)

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef struct ContourID

cdef struct ContourID:
    np.int64_t contour_id
    ContourID *parent
    ContourID *next
    ContourID *prev
    np.int64_t count

cdef struct CandidateContour

cdef struct CandidateContour:
    np.int64_t contour_id
    np.int64_t join_id
    CandidateContour *next
