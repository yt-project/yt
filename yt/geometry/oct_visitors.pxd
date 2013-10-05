"""
Oct visitor definitions file




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np

cdef struct Oct
cdef struct Oct:
    np.int64_t file_ind     # index with respect to the order in which it was
                            # added
    np.int64_t domain_ind   # index within the global set of domains
    np.int64_t domain       # (opt) addl int index
    Oct **children          # Up to 8 long

cdef struct OctVisitorData:
    np.uint64_t index
    np.uint64_t last
    np.int64_t global_index
    np.int64_t pos[3]       # position in ints
    np.uint8_t ind[3]              # cell position
    void *array
    int dims
    np.int32_t domain
    np.int8_t level
    np.int8_t oref # This is the level of overref.  1 => 8 zones, 2 => 64, etc.
                   # To calculate nzones, 1 << (oref * 3)
    np.int32_t nz
                            
ctypedef void oct_visitor_function(Oct *, OctVisitorData *visitor,
                                   np.uint8_t selected)

cdef oct_visitor_function count_total_octs
cdef oct_visitor_function count_total_cells
cdef oct_visitor_function mark_octs
cdef oct_visitor_function mask_octs
cdef oct_visitor_function index_octs
cdef oct_visitor_function icoords_octs
cdef oct_visitor_function ires_octs
cdef oct_visitor_function fcoords_octs
cdef oct_visitor_function fwidth_octs
cdef oct_visitor_function copy_array_f64
cdef oct_visitor_function copy_array_i64
cdef oct_visitor_function identify_octs
cdef oct_visitor_function assign_domain_ind
cdef oct_visitor_function fill_file_indices_oind
cdef oct_visitor_function fill_file_indices_rind
cdef oct_visitor_function count_by_domain
cdef oct_visitor_function store_octree
cdef oct_visitor_function load_octree

cdef inline int cind(int i, int j, int k):
    # THIS ONLY WORKS FOR CHILDREN.  It is not general for zones.
    return (((i*2)+j)*2+k)

cdef inline int oind(OctVisitorData *data):
    cdef int d = (1 << data.oref)
    return (((data.ind[0]*d)+data.ind[1])*d+data.ind[2])

cdef inline int rind(OctVisitorData *data):
    cdef int d = (1 << data.oref)
    return (((data.ind[2]*d)+data.ind[1])*d+data.ind[0])
