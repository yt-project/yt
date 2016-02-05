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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.float64_t dot(const np.float64_t* a, 
                      const np.float64_t* b) nogil:
    cdef np.int64_t i
    cdef np.float64_t rv = 0.0
    for i in range(3):
        rv += a[i]*b[i]
    return rv


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void cross(const np.float64_t* a, 
                const np.float64_t* b,
                np.float64_t* c) nogil:
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]


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
    cdef np.int64_t _ray_triangle_intersect(self, Ray* ray, np.int64_t tri_index):
        # cribbed from 3D Math Primer for Graphics and Game Development, section A.16
        cdef np.float64_t p0[3]
        cdef np.float64_t p1[3]
        cdef np.float64_t p2[3]
        cdef np.float64_t e1[3]
        cdef np.float64_t e2[3]

        cdef Triangle* tri
        tri = &(self.triangles[tri_index])
        
        cdef int i
        for i in range(3):
            p0[i] = self.vertices[tri.v0][i]
            p1[i] = self.vertices[tri.v1][i]
            p2[i] = self.vertices[tri.v2][i]
            e1[i] = p1[i] - p0[i];
            e2[i] = p2[i] - p1[i];
        
        cdef np.float64_t N[3]
        cross(e1, e1, N);
        cdef np.float64_t dotprod = dot(N, ray.direction)

        if not dotprod < 0.0:
            # ray is travelling the wrong direction
            return False

        cdef np.float64_t d = dot(N, p0)
        cdef np.float64_t t = d - dot(N, ray.origin) 

        if not t <= 0.0:
            # ray origin is behind triangle
            return False

        t /= dotprod

        if not t >= ray.t_far:
            # closer intersection already found
            return False

        assert(t >= 0.0)
        assert(t <= ray.t_far)

        cdef np.float64_t p[3]
        for i in range(3):
            p[i] = ray.origin[i] + ray.direction[i]*t
        
        cdef np.float64_t u0, u1, u2, v0, v1, v2
        
        if fabs(N[0]) > fabs(N[1]):
            if fabs(N[0]) > fabs(N[2]):
                u0 = p[1]  - p0[1]
                u1 = p1[1] - p0[1]
                u2 = p2[1] - p0[1]

                v0 = p[2]  - p0[2]
                v1 = p1[2] - p0[2]
                v2 = p2[2] - p0[2]
    
            else:
                u0 = p[0]  - p0[0]
                u1 = p1[0] - p0[0]
                u2 = p2[0] - p0[0]

                v0 = p[1]  - p0[1]
                v1 = p1[1] - p0[1]
                v2 = p2[1] - p0[1]
        else:
            if fabs(N[1]) > fabs(N[2]):
                u0 = p[0]  - p0[0]
                u1 = p1[0] - p0[0]
                u2 = p2[0] - p0[0]

                v0 = p[2]  - p0[2]
                v1 = p1[2] - p0[2]
                v2 = p2[2] - p0[2]
            else:
                u0 = p[0]  - p0[0]
                u1 = p1[0] - p0[0]
                u2 = p2[0] - p0[0]

                v0 = p[1]  - p0[1]
                v1 = p1[1] - p0[1]
                v2 = p2[1] - p0[1]

        cdef np.float64_t temp = u1*v2 - u2*v1

        if not (temp is not  0.0):
            # denominator is invalid
            return False

        temp = 1.0/temp

        cdef np.float64_t a = (u0*v2 - u2*v0)*temp
        if not (a >= 0.0):
            # barycentric coord out of range
            return False

        cdef np.float64_t b = (u1*v0 - u0*v1)*temp
        if not (b >= 0.0):
            # barycentric coord out of range
            return False

        cdef np.float64_t c = 1.0 - a - b
        if not (c >= 0.0):
            # barycentric coord out of range
            return False

        # we have a hit, update ray
        ray.t_far = t
        ray.elem_id = tri.elem_id

        return True

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
