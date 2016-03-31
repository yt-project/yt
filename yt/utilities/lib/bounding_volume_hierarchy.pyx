cimport cython 
import numpy as np
cimport numpy as np
from libc.math cimport fabs, fmax, fmin
from libc.stdlib cimport malloc, free
from cython.parallel import parallel, prange
from vec3_ops cimport dot, subtract, cross
from yt.utilities.lib.element_mappings cimport ElementSampler, \
    Q1Sampler3D

cdef ElementSampler Q1Sampler = Q1Sampler3D()

cdef extern from "mesh_construction.h":
    enum:
        MAX_NUM_TRI
    int triangulate_hex[MAX_NUM_TRI][3]

# define some constants
cdef np.float64_t DETERMINANT_EPS = 1.0e-10
cdef np.float64_t INF = np.inf
cdef np.int64_t   LEAF_SIZE = 16


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_triangle_intersect(Ray* ray, const Triangle* tri) nogil:
# https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

    # edge vectors
    cdef np.float64_t e1[3]
    cdef np.float64_t e2[3]
    subtract(tri.p1, tri.p0, e1)
    subtract(tri.p2, tri.p0, e2)

    cdef np.float64_t P[3]
    cross(ray.direction, e2, P)

    cdef np.float64_t det, inv_det
    det = dot(e1, P)
    if(det > -DETERMINANT_EPS and det < DETERMINANT_EPS): 
        return False
    inv_det = 1.0 / det

    cdef np.float64_t T[3]
    subtract(ray.origin, tri.p0, T)

    cdef np.float64_t u = dot(T, P) * inv_det
    if(u < 0.0 or u > 1.0):
        return False

    cdef np.float64_t Q[3]
    cross(T, e1, Q)

    cdef np.float64_t v = dot(ray.direction, Q) * inv_det
    if(v < 0.0 or u + v  > 1.0):
        return False

    cdef np.float64_t t = dot(e2, Q) * inv_det

    if(t > DETERMINANT_EPS and t < ray.t_far):
        ray.t_far = t
        ray.elem_id = tri.elem_id
        return True

    return False


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_bbox_intersect(Ray* ray, const BBox bbox) nogil:
# https://tavianator.com/fast-branchless-raybounding-box-intersections/

    cdef np.float64_t tmin = -INF
    cdef np.float64_t tmax =  INF
 
    cdef np.int64_t i
    cdef np.float64_t t1, t2
    for i in range(3):
        t1 = (bbox.left_edge[i]  - ray.origin[i])*ray.inv_dir[i]
        t2 = (bbox.right_edge[i] - ray.origin[i])*ray.inv_dir[i] 
        tmin = fmax(tmin, fmin(t1, t2))
        tmax = fmin(tmax, fmax(t1, t2))
 
    return tmax >= fmax(tmin, 0.0)


cdef class BVH:
    '''

    This class implements a bounding volume hierarchy (BVH), a spatial acceleration
    structure for fast ray-tracing. A BVH is like a kd-tree, except that instead of 
    partitioning the *volume* of the parent to create the children, we partition the 
    triangles themselves into 'left' or 'right' sub-trees. The bounding volume for a
    node is then determined by computing the bounding volume of the triangles that
    belong to it. This allows us to quickly discard triangles that are not close 
    to intersecting a given ray.

    '''

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self,
                  np.float64_t[:, ::1] vertices,
                  np.int64_t[:, ::1] indices,
                  np.float64_t[:, ::1] field_data):

        self.vertices = vertices
        self.indices = indices
        self.field_data = field_data

        cdef np.int64_t num_elem = indices.shape[0]
        cdef np.int64_t num_verts_per_elem = indices.shape[1]
        cdef np.int64_t num_tri = 12*num_elem

        # fill our array of triangles
        cdef np.int64_t i, j, k
        cdef np.int64_t offset, tri_index
        cdef np.int64_t v0, v1, v2
        cdef Triangle* tri
        self.triangles = <Triangle*> malloc(num_tri * sizeof(Triangle))
        for i in range(num_elem):
            offset = 12*i
            for j in range(12):
                tri_index = offset + j
                tri = &(self.triangles[tri_index])
                tri.elem_id = i
                v0 = indices[i][triangulate_hex[j][0]]
                v1 = indices[i][triangulate_hex[j][1]]
                v2 = indices[i][triangulate_hex[j][2]]
                for k in range(3):
                    tri.p0[k] = vertices[v0][k]
                    tri.p1[k] = vertices[v1][k]
                    tri.p2[k] = vertices[v2][k]
                    tri.centroid[k] = (tri.p0[k] + tri.p1[k] + tri.p2[k]) / 3.0
                    tri.bbox.left_edge[k]  = fmin(fmin(tri.p0[k], tri.p1[k]), tri.p2[k])
                    tri.bbox.right_edge[k] = fmax(fmax(tri.p0[k], tri.p1[k]), tri.p2[k])

        self.root = self._recursive_build(0, num_tri)

    cdef void _recursive_free(self, BVHNode* node) nogil:
        if node.end - node.begin > LEAF_SIZE:
            self._recursive_free(node.left)
            self._recursive_free(node.right)
        free(node)

    def __dealloc__(self):
        self._recursive_free(self.root)
        free(self.triangles)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.int64_t _partition(self, np.int64_t begin, np.int64_t end,
                               np.int64_t ax, np.float64_t split) nogil:
        # this re-orders the triangle array so that all of the triangles 
        # to the left of mid have centroids less than or equal to "split" 
        # along the direction "ax". All the triangles to the right of mid 
        # will have centroids *greater* than "split" along "ax".
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
                             np.int64_t begin, np.int64_t end) nogil:
        cdef np.int64_t i, j
        cdef BBox box = self.triangles[begin].bbox
        for i in range(begin+1, end):
            for j in range(3):
                box.left_edge[j] = fmin(box.left_edge[j],
                                        self.triangles[i].bbox.left_edge[j])
                box.right_edge[j] = fmax(box.right_edge[j], 
                                         self.triangles[i].bbox.right_edge[j])
        node.bbox = box

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void intersect(self, Ray* ray) nogil:
        self._recursive_intersect(ray, self.root)
        
        if ray.elem_id < 0:
            return

        cdef np.float64_t[3] position
        cdef np.float64_t[3] direction
        cdef np.float64_t[8] field_data
        cdef np.int64_t[8] element_indices
        cdef np.float64_t[24] vertices

        cdef np.int64_t i
        for i in range(3):
            position[i] = ray.origin[i] + ray.t_far*ray.direction[i]
            
        for i in range(8):
            element_indices[i] = self.indices[ray.elem_id][i]
            field_data[i]      = self.field_data[ray.elem_id][i]

        for i in range(8):
            vertices[i*3]     = self.vertices[element_indices[i]][0]
            vertices[i*3 + 1] = self.vertices[element_indices[i]][1]
            vertices[i*3 + 2] = self.vertices[element_indices[i]][2]   

        cdef double mapped_coord[3]
        Q1Sampler.map_real_to_unit(mapped_coord, vertices, position)
        val = Q1Sampler.sample_at_unit_point(mapped_coord, field_data)
        ray.data_val = val

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _recursive_intersect(self, Ray* ray, BVHNode* node) nogil:

        # check for bbox intersection:
        if not ray_bbox_intersect(ray, node.bbox):
            return

        # check for leaf
        cdef np.int64_t i, hit
        cdef Triangle* tri
        if (node.end - node.begin) <= LEAF_SIZE:
            for i in range(node.begin, node.end):
                tri = &(self.triangles[i])
                hit = ray_triangle_intersect(ray, tri)
            return

        # if not leaf, intersect with left and right children
        self._recursive_intersect(ray, node.left)
        self._recursive_intersect(ray, node.right)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef BVHNode* _recursive_build(self, np.int64_t begin, np.int64_t end) nogil:
        cdef BVHNode *node = <BVHNode* > malloc(sizeof(BVHNode))
        node.begin = begin
        node.end = end

        self._get_node_bbox(node, begin, end)
        
        # check for leaf
        if (end - begin) <= LEAF_SIZE:
            return node
        
        # we use the "split in the middle of the longest axis approach"
        # see: http://www.vadimkravcenko.com/bvh-tree-building/

        # compute longest dimension
        cdef np.int64_t ax = 0
        cdef np.float64_t d = fabs(node.bbox.right_edge[0] - 
                                   node.bbox.left_edge[0])
        if fabs(node.bbox.right_edge[1] - node.bbox.left_edge[1]) > d:
            ax = 1
        if fabs(node.bbox.right_edge[2] - node.bbox.left_edge[2]) > d:
            ax = 2

        # split in half along that dimension
        cdef np.float64_t split = 0.5*(node.bbox.right_edge[ax] +
                                       node.bbox.left_edge[ax])

        # sort triangle list
        cdef np.int64_t mid = self._partition(begin, end, ax, split)

        if(mid == begin or mid == end):
            mid = begin + (end-begin)/2
        
        # recursively build sub-trees
        node.left = self._recursive_build(begin, mid)
        node.right = self._recursive_build(mid, end)

        return node


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void cast_rays(np.float64_t* image,
                    const np.float64_t* origins,
                    const np.float64_t* direction,
                    const int N, 
                    BVH bvh) nogil:

    cdef Ray* ray 
    cdef int i, j, k
    
    with nogil, parallel():
       
        ray = <Ray *> malloc(sizeof(Ray))
    
        for k in range(3):
            ray.direction[k] = direction[k]
            ray.inv_dir[k] = 1.0 / direction[k]

        for i in prange(N):
            for j in range(3):
                ray.origin[j] = origins[N*j + i]
            ray.t_far = INF
            ray.t_near = 0.0
            ray.data_val = 0
            ray.elem_id = -1
            bvh.intersect(ray)
            image[i] = ray.data_val

        free(ray)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def test_ray_trace(np.ndarray[np.float64_t, ndim=1] image, 
                   np.ndarray[np.float64_t, ndim=2] origins,
                   np.ndarray[np.float64_t, ndim=1] direction,
                   BVH bvh):
    
    cdef int N = origins.shape[0]
    cast_rays(&image[0], &origins[0, 0], &direction[0], N, bvh)
