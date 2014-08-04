"""
Contour finding exports



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

cdef ContourID *contour_create(np.int64_t contour_id,
                               ContourID *prev = ?)
cdef void contour_delete(ContourID *node)
cdef ContourID *contour_find(ContourID *node)
cdef void contour_union(ContourID *node1, ContourID *node2)
cdef int candidate_contains(CandidateContour *first,
                            np.int64_t contour_id,
                            np.int64_t join_id = ?)
cdef CandidateContour *candidate_add(CandidateContour *first,
                                     np.int64_t contour_id,
                                     np.int64_t join_id = ?)
