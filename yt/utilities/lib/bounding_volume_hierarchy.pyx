cimport cython 
import numpy as np
cimport numpy as np
from libc.math cimport fabs, fmax, fmin
from libc.stdlib cimport malloc, free

cdef extern from "mesh_construction.h":
    enum:
        MAX_NUM_TRI
    int triangulate_hex[MAX_NUM_TRI][3]

# ray data structure
cdef struct Ray:
    np.float64_t origin[3]
    np.float64_t direction[3]
    np.float64_t data_val
    np.float64_t t_near
    np.float64_t t_far
    np.int64_t  elem_id

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
    np.int64_t v0, v1, v2
    np.int64_t elem_id
    np.float64_t centroid[3]
    BBox bbox

cdef class BVH:
    cdef BVHNode* root
    cdef Triangle* triangles
    cdef np.float64_t[:, ::1] vertices

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __init__(self,
                 np.float64_t[:, ::1] vertices,
                 np.int64_t[:, ::1] indices):
        
        self.vertices = vertices
        cdef np.int64_t num_elem = indices.shape[0]
        cdef np.int64_t num_tri = 12*num_elem
        self.triangles = <Triangle*> malloc(num_tri * sizeof(Triangle))
        
        cdef np.int64_t i, j, k
        cdef np.int64_t offset, tri_index
        cdef np.float64_t[:] p0
        cdef np.float64_t[:] p1
        cdef np.float64_t[:] p2
        cdef Triangle* tri
        for i in range(num_elem):
            offset = 12*i
            for j in range(12):
                tri_index = offset + j
                tri = &(self.triangles[tri_index])
                tri.elem_id = i
                tri.v0 = indices[i][triangulate_hex[j][0]]
                tri.v1 = indices[i][triangulate_hex[j][1]]
                tri.v2 = indices[i][triangulate_hex[j][2]]                
                p0 = vertices[tri.v0]
                p1 = vertices[tri.v1]
                p2 = vertices[tri.v2]
                for k in range(3):
                    tri.centroid[k] = (p0[k] + p1[k] + p2[k]) / 3.0
                    tri.bbox.left_edge[k] = fmin(fmin(p0[k], p1[k]), p2[k])
                    tri.bbox.right_edge[k] = fmax(fmax(p0[k], p1[k]), p2[k])
                    
        self.root = self.build(0, num_tri)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.int64_t partition(self, np.int64_t begin, np.int64_t end,
                              np.int64_t ax, np.float64_t split):
        cdef np.int64_t mid = begin
        while (begin != end):
            if self.triangles[mid].centroid[ax] > split:
                mid += 1
            elif self.triangles[begin].centroid[ax] > split:
                self.triangles[mid], self.triangles[begin] = \
                self.triangles[begin], self.triangles[mid]
                mid += 1
            begin += 1
        return mid
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _get_node_bbox(self, BVHNode* node, 
                             np.int64_t begin, np.int64_t end):
        cdef np.int64_t i, j
        cdef BBox box = self.triangles[begin].bbox
        for i in range(begin+1, end):
            for j in range(3):
                box.left_edge[j] = fmax(box.left_edge[j], 
                                        self.triangles[i].bbox.left_edge[j])
                box.right_edge[j] = fmin(box.right_edge[j], 
                                         self.triangles[i].bbox.right_edge[j])
        node.bbox = box       
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef BVHNode* build(self, np.int64_t begin, np.int64_t end):
        cdef BVHNode *node = <BVHNode* > malloc(sizeof(BVHNode))
        node.begin = begin
        node.end = end

        self._get_node_bbox(node, begin, end)
        
        if (end - begin) == 1:
            return node
        
        cdef np.int64_t ax = 0
        cdef np.float64_t d = fabs(node.bbox.right_edge[0] - 
                                   node.bbox.left_edge[0])
        if fabs(node.bbox.right_edge[1] - node.bbox.left_edge[1]) > d:
            ax = 1
        if fabs(node.bbox.right_edge[2] - node.bbox.left_edge[2]) > d:
            ax = 2

        cdef np.float64_t split = 0.5*(node.bbox.right_edge[ax] - 
                                       node.bbox.left_edge[ax])

        cdef np.int64_t mid = self.partition(begin, end, ax, split)
        if(mid == begin or mid == end):
            mid = begin + (end-begin)/2

        node.left = self.build(begin, mid)
        node.right = self.build(mid, end)

        return node
