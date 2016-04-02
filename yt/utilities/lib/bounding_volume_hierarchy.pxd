cimport cython 
import numpy as np
cimport numpy as np
from yt.utilities.lib.element_mappings cimport ElementSampler

# ray data structure
cdef struct Ray:
    np.float64_t origin[3]
    np.float64_t direction[3]
    np.float64_t inv_dir[3]
    np.float64_t data_val
    np.float64_t t_near
    np.float64_t t_far
    np.int64_t elem_id
    np.int64_t near_boundary

# axis-aligned bounding box
cdef struct BBox:
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]

# node for the bounding volume hierarchy
cdef struct BVHNode:
    np.int64_t begin
    np.int64_t end
    BVHNode* left
    BVHNode* right
    BBox bbox
    
# triangle data structure
cdef struct Triangle:
    np.float64_t p0[3]
    np.float64_t p1[3]
    np.float64_t p2[3]
    np.int64_t elem_id
    np.float64_t centroid[3]
    BBox bbox

cdef class BVH:
    cdef BVHNode* root
    cdef Triangle* triangles
    cdef np.float64_t* vertices
    cdef np.float64_t* field_data
    cdef np.int64_t num_tri_per_elem
    cdef np.int64_t num_tri
    cdef np.int64_t num_elem
    cdef np.int64_t num_verts_per_elem
    cdef ElementSampler sampler
    cdef np.int64_t _partition(self, np.int64_t begin, np.int64_t end,
                               np.int64_t ax, np.float64_t split) nogil
    cdef void intersect(self, Ray* ray) nogil
    cdef void _get_node_bbox(self, BVHNode* node, 
                             np.int64_t begin, np.int64_t end) nogil
    cdef void _recursive_intersect(self, Ray* ray, BVHNode* node) nogil
    cdef BVHNode* _recursive_build(self, np.int64_t begin, np.int64_t end) nogil
    cdef void _recursive_free(self, BVHNode* node) nogil
