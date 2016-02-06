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
    np.float64_t inv_dir[3]
    np.float64_t data_val
    np.float64_t t_near
    np.float64_t t_far
    np.int64_t elem_id

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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t ray_triangle_intersect(Ray* ray, const Triangle* tri):
    # guts of this are cribbed from "3D Math Primer for Graphics 
    # and Game Development", section A.16

    # edge vectors
    cdef np.float64_t e1[3]
    cdef np.float64_t e2[3]
    cdef int i
    for i in range(3):
        e1[i] = tri.p1[i] - tri.p0[i];
        e2[i] = tri.p2[i] - tri.p1[i];
        
    # normal vector
    cdef np.float64_t N[3]
    cross(e1, e1, N);

    cdef np.float64_t dotprod = dot(N, ray.direction)
    if not dotprod < 0.0:
        # ray is travelling the wrong direction
        return False

    # distance along ray
    cdef np.float64_t d = dot(N, tri.p0)
    cdef np.float64_t t = d - dot(N, ray.origin) 

    if not t <= 0.0:
        # ray origin is behind triangle
        return False

    t /= dotprod

    if not t >= ray.t_far:
        # closer intersection already found
        return False

#    assert(t >= 0.0)
#    assert(t <= ray.t_far)
    
    # 3D point of intersection
    cdef np.float64_t p[3]
    for i in range(3):
        p[i] = ray.origin[i] + ray.direction[i]*t
        
    cdef np.float64_t u0, u1, u2, v0, v1, v2
        
    if fabs(N[0]) > fabs(N[1]):
        if fabs(N[0]) > fabs(N[2]):
            u0 = p[1]      - tri.p0[1]
            u1 = tri.p1[1] - tri.p0[1]
            u2 = tri.p2[1] - tri.p0[1]

            v0 = p[2]      - tri.p0[2]
            v1 = tri.p1[2] - tri.p0[2]
            v2 = tri.p2[2] - tri.p0[2]
    
        else:
            u0 = p[0]      - tri.p0[0]
            u1 = tri.p1[0] - tri.p0[0]
            u2 = tri.p2[0] - tri.p0[0]
            
            v0 = p[1]      - tri.p0[1]
            v1 = tri.p1[1] - tri.p0[1]
            v2 = tri.p2[1] - tri.p0[1]
    else:
        if fabs(N[1]) > fabs(N[2]):
            u0 = p[0]      - tri.p0[0]
            u1 = tri.p1[0] - tri.p0[0]
            u2 = tri.p2[0] - tri.p0[0]

            v0 = p[2]      - tri.p0[2]
            v1 = tri.p1[2] - tri.p0[2]
            v2 = tri.p2[2] - tri.p0[2]
        else:
            u0 = p[0]      - tri.p0[0]
            u1 = tri.p1[0] - tri.p0[0]
            u2 = tri.p2[0] - tri.p0[0]

            v0 = p[1]      - tri.p0[1]
            v1 = tri.p1[1] - tri.p0[1]
            v2 = tri.p2[1] - tri.p0[1]

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
cdef np.int64_t ray_bbox_intersect(Ray* ray, const BBox bbox):

    cdef np.float64_t tx1 = (bbox.left_edge[0]  - 
                             ray.origin[0])*ray.inv_dir[0]
    cdef np.float64_t tx2 = (bbox.right_edge[0] -
                             ray.origin[0])*ray.inv_dir[0]
 
    cdef np.float64_t tmin = fmin(tx1, tx2)
    cdef np.float64_t tmax = fmax(tx1, tx2)

    cdef np.float64_t ty1 = (bbox.left_edge[1]  -
                             ray.origin[1])*ray.inv_dir[1]
    cdef np.float64_t ty2 = (bbox.right_edge[1] -
                             ray.origin[1])*ray.inv_dir[1]
 
    tmin = fmax(tmin, fmin(ty1, ty2))
    tmax = fmin(tmax, fmax(ty1, ty2))

    cdef np.float64_t tz1 = (bbox.left_edge[2]  -
                             ray.origin[2])*ray.inv_dir[2]
    cdef np.float64_t tz2 = (bbox.right_edge[2] -
                             ray.origin[2])*ray.inv_dir[2]
 
    tmin = fmax(tmin, fmin(tz1, tz2))
    tmax = fmin(tmax, fmax(tz1, tz2))

    if not (tmax > 0):
        return False
 
    return tmax >= tmin;

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
                    tri.bbox.left_edge[k] = fmin(fmin(tri.p0[k], tri.p1[k]), tri.p2[k])
                    tri.bbox.right_edge[k] = fmax(fmax(tri.p0[k], tri.p1[k]), tri.p2[k])

        # recursive build
        self.root = self._build(0, num_tri)

        self.intersect()

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
    cdef void intersect(self):
        cdef Ray ray
        ray.origin[0]     = -5.0
        ray.origin[1]     = 0.0
        ray.origin[2]     = 1.0
        ray.direction[0]  = 1.0
        ray.direction[1]  = 0.0
        ray.direction[2]  = 0.0
        ray.t_near        = 0.0
        ray.t_far         = 1e300
        
        cdef int i 
        for i in range(3):
            ray.inv_dir[i] = 1.0 / ray.direction[i]
        
        self._recursive_intersect(&ray, self.root)
        print ray.elem_id

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _recursive_intersect(self, Ray* ray, BVHNode* node):

        # check for bbox intersection:
        if not ray_bbox_intersect(ray, node.bbox):
            return

        # check for leaf
        cdef np.int64_t i, hit
        cdef Triangle* tri
        if (node.end - node.begin) == 1:
            print 'testing triangles'
            for i in range(node.begin, node.end):
                tri = &(self.triangles[i])
                hit = ray_triangle_intersect(ray, tri)
                if hit:
                    print hit

        # if not leaf, intersect with left and right children
        self._recursive_intersect(ray, node.left)
        self._recursive_intersect(ray, node.right)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef BVHNode* _build(self, np.int64_t begin, np.int64_t end):
        cdef BVHNode *node = <BVHNode* > malloc(sizeof(BVHNode))
        node.begin = begin
        node.end = end

        self._get_node_bbox(node, begin, end)
        
        # check for leaf
        if (end - begin) == 1:
            return node
        
        # compute longest dimension
        cdef np.int64_t ax = 0
        cdef np.float64_t d = fabs(node.bbox.right_edge[0] - 
                                   node.bbox.left_edge[0])
        if fabs(node.bbox.right_edge[1] - node.bbox.left_edge[1]) > d:
            ax = 1
        if fabs(node.bbox.right_edge[2] - node.bbox.left_edge[2]) > d:
            ax = 2

        # split in half along that dimension
        cdef np.float64_t split = 0.5*(node.bbox.right_edge[ax] - 
                                       node.bbox.left_edge[ax])

        # sort triangle list
        cdef np.int64_t mid = self.partition(begin, end, ax, split)
        if(mid == begin or mid == end):
            mid = begin + (end-begin)/2

        node.left = self._build(begin, mid)
        node.right = self._build(mid, end)

        return node
