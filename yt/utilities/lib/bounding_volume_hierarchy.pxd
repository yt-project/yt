cimport cython 
import numpy as np
cimport numpy as np
from yt.utilities.lib.element_mappings cimport ElementSampler
from yt.utilities.lib.primitives cimport BBox, Ray

cdef extern from "mesh_triangulation.h":
    enum:
        MAX_NUM_TRI

    int HEX_NV
    int HEX_NT
    int TETRA_NV
    int TETRA_NT
    int WEDGE_NV
    int WEDGE_NT
    int triangulate_hex[MAX_NUM_TRI][3]
    int triangulate_tetra[MAX_NUM_TRI][3]
    int triangulate_wedge[MAX_NUM_TRI][3]
    int hex20_faces[6][8]

# node for the bounding volume hierarchy
cdef struct BVHNode:
    np.int64_t begin
    np.int64_t end
    BVHNode* left
    BVHNode* right
    BBox bbox

# pointer to function that computes primitive intersection
ctypedef np.int64_t (*intersect_func_type)(const void* primitives,
                                           const np.int64_t item,
                                           Ray* ray) nogil

# pointer to function that computes primitive centroids
ctypedef void (*centroid_func_type)(const void *primitives,
                                    const np.int64_t item,
                                    np.float64_t[3] centroid) nogil

# pointer to function that computes primitive bounding boxes
ctypedef void (*bbox_func_type)(const void *primitives,
                                const np.int64_t item,
                                BBox* bbox) nogil


cdef class BVH:
    cdef BVHNode* root
    cdef void* primitives
    cdef np.int64_t* prim_ids
    cdef np.float64_t** centroids
    cdef BBox* bboxes
    cdef np.float64_t* vertices
    cdef np.float64_t* field_data
    cdef np.int64_t num_prim_per_elem
    cdef np.int64_t num_prim
    cdef np.int64_t num_elem
    cdef np.int64_t num_verts_per_elem
    cdef np.int64_t num_field_per_elem
    cdef int[MAX_NUM_TRI][3] tri_array
    cdef ElementSampler sampler
    cdef centroid_func_type get_centroid
    cdef bbox_func_type get_bbox
    cdef intersect_func_type get_intersect
    cdef np.int64_t _partition(self, np.int64_t begin, np.int64_t end,
                               np.int64_t ax, np.float64_t split) nogil
    cdef void _set_up_triangles(self,
                                np.float64_t[:, :] vertices,
                                np.int64_t[:, :] indices) nogil
    cdef void _set_up_patches(self,
                              np.float64_t[:, :] vertices,
                              np.int64_t[:, :] indices) nogil
    cdef void intersect(self, Ray* ray) nogil
    cdef void _get_node_bbox(self, BVHNode* node, 
                             np.int64_t begin, np.int64_t end) nogil
    cdef void _recursive_intersect(self, Ray* ray, BVHNode* node) nogil
    cdef BVHNode* _recursive_build(self, np.int64_t begin, np.int64_t end) nogil
    cdef void _recursive_free(self, BVHNode* node) nogil
