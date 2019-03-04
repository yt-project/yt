cimport numpy as np
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool
from libc.stdint cimport uint32_t, uint64_t, int64_t, int32_t

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        pass

    cdef cppclass ostream:
        pass

# the following extern definitions adapted from
# http://stackoverflow.com/a/31009461/1382869
    
# obviously std::ios_base isn't a namespace, but this lets
# Cython generate the correct C++ code
cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    # you can define other constants as needed

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+

    cdef cppclass ifstream(istream):
        # constructors
        ifstream(const char*) except +
        ifstream(const char*, open_mode) except+


cdef extern from "c_kdtree.hpp":
    cdef cppclass Node:
        bool is_leaf
        uint32_t leafid
        uint32_t ndim
        double *left_edge
        double *right_edge
        uint64_t left_idx
        uint64_t children
        uint32_t split_dim
        double split
        Node* less
        Node* greater
        bool *periodic_left
        bool *periodic_right
        vector[vector[uint32_t]] left_neighbors
        vector[vector[uint32_t]] right_neighbors
        vector[uint32_t] all_neighbors
    cdef cppclass KDTree:
        uint64_t* all_idx
        uint64_t npts
        uint32_t ndim
        int64_t data_version
        uint32_t leafsize
        double* domain_left_edge
        double* domain_right_edge
        double* domain_width
        double* domain_mins
        double* domain_maxs
        bool* periodic
        bool any_periodic
        uint32_t num_leaves
        vector[Node*] leaves
        Node* root
        KDTree(double *pts, uint64_t *idx, uint64_t n, uint32_t m,
               uint32_t leafsize0, double *left_edge, double *right_edge,
               bool *periodic, int64_t data_version)
        KDTree(double *pts, uint64_t *idx, uint64_t n, uint32_t m,
               uint32_t leafsize0, double *left_edge, double *right_edge,
               bool *periodic, int64_t data_version,
               bool use_sliding_midpoint)
        KDTree(istream &ist)
        void serialize(ostream &os)
        double* wrap_pos(double* pos) nogil
        vector[uint32_t] get_neighbor_ids(double* pos) nogil
        Node* search(double* pos) nogil
        Node* search(double* pos, bool dont_wrap) nogil
        void consolidate_edges(double *leaves_le, double *leaves_re)


cdef class PyNode:
    cdef Node *_node
    cdef readonly np.uint32_t id
    cdef readonly np.uint64_t npts
    cdef readonly np.uint32_t ndim
    cdef readonly np.uint32_t num_leaves
    cdef readonly np.uint64_t start_idx
    cdef readonly np.uint64_t stop_idx
    cdef double *_domain_width
    cdef readonly object left_neighbors, right_neighbors
    cdef void _init_node(self, Node* node, uint32_t num_leaves,
                         double *domain_width)

cdef class PyKDTree:
    cdef KDTree *_tree
    cdef readonly uint64_t npts
    cdef readonly uint32_t ndim
    cdef readonly uint32_t num_leaves
    cdef readonly uint32_t leafsize
    cdef readonly int64_t data_version
    cdef double *_left_edge
    cdef double *_right_edge
    cdef bool *_periodic
    cdef readonly object leaves
    cdef readonly object _idx
    cdef void _init_tree(self, KDTree* tree)
    cdef void _make_tree(self, double *pts, bool use_sliding_midpoint)
    cdef void _make_leaves(self)
    cdef np.ndarray[np.uint32_t, ndim=1] _get_neighbor_ids(self, np.ndarray[double, ndim=1] pos)
    cdef np.ndarray[np.uint32_t, ndim=1] _get_neighbor_ids_3(self, np.float64_t pos[3])
    cdef PyNode _get(self, np.ndarray[double, ndim=1] pos)
