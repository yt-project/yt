cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport fabs
from libc.stdlib cimport malloc, free
from cython.parallel import parallel, prange

from yt.utilities.lib.primitives cimport \
    Triangle, \
    ray_triangle_intersect, \
    triangle_centroid, \
    triangle_bbox
from yt.utilities.lib.element_mappings cimport \
    ElementSampler, \
    Q1Sampler3D, \
    P1Sampler3D, \
    W1Sampler3D

cdef ElementSampler Q1Sampler = Q1Sampler3D()
cdef ElementSampler P1Sampler = P1Sampler3D()
cdef ElementSampler W1Sampler = W1Sampler3D()

cdef extern from "platform_dep.h" nogil:
    double fmax(double x, double y)
    double fmin(double x, double y)

cdef extern from "mesh_construction.h":
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

# define some constants
cdef np.float64_t INF = np.inf
cdef np.int64_t   LEAF_SIZE = 16


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
                  np.float64_t[:, :] vertices,
                  np.int64_t[:, :] indices,
                  np.float64_t[:, :] field_data):

        self.num_elem = indices.shape[0]
        self.num_verts_per_elem = indices.shape[1]
        self.num_field_per_elem = field_data.shape[1]

        self.get_centroid = triangle_centroid
        self.get_bbox = triangle_bbox
        self.get_intersect = ray_triangle_intersect

        # We need to figure out what kind of elements we've been handed.
        cdef int[MAX_NUM_TRI][3] tri_array
        if self.num_verts_per_elem == 8:
            self.num_tri_per_elem = HEX_NT
            tri_array = triangulate_hex
            self.sampler = Q1Sampler
        elif self.num_verts_per_elem == 6:
            self.num_tri_per_elem = WEDGE_NT
            tri_array = triangulate_wedge
            self.sampler = W1Sampler
        elif self.num_verts_per_elem == 4:
            self.num_tri_per_elem = TETRA_NT
            tri_array = triangulate_tetra
            self.sampler = P1Sampler
        self.num_prim = self.num_tri_per_elem*self.num_elem

        # allocate storage
        cdef np.int64_t v_size = self.num_verts_per_elem * self.num_elem * 3
        self.vertices = <np.float64_t*> malloc(v_size * sizeof(np.float64_t))
        cdef np.int64_t f_size = self.num_field_per_elem * self.num_elem
        self.field_data = <np.float64_t*> malloc(f_size * sizeof(np.float64_t))

        # create data buffers
        cdef np.int64_t i, j, k
        cdef np.int64_t field_offset, vertex_offset
        for i in range(self.num_elem):
            for j in range(self.num_verts_per_elem):
                vertex_offset = i*self.num_verts_per_elem*3 + j*3
                for k in range(3):
                    self.vertices[vertex_offset + k] = vertices[indices[i,j]][k]
            field_offset = i*self.num_field_per_elem
            for j in range(self.num_field_per_elem):
                self.field_data[field_offset + j] = field_data[i][j]                

        # fill our array of primitives
        cdef np.int64_t offset, tri_index
        cdef np.int64_t v0, v1, v2
        cdef Triangle* tri
        self.primitives = malloc(self.num_prim * sizeof(Triangle))
        self.centroids = <np.float64_t**> malloc(self.num_prim * sizeof(np.float64_t*))
        for i in range(self.num_prim):
            self.centroids[i] = <np.float64_t*> malloc(3*sizeof(np.float64_t))
        self.bboxes = <BBox*> malloc(self.num_prim * sizeof(BBox))
        for i in range(self.num_elem):
            offset = self.num_tri_per_elem*i
            for j in range(self.num_tri_per_elem):
                tri_index = offset + j
                tri = &(<Triangle*> self.primitives)[tri_index]
                tri.elem_id = i
                v0 = indices[i][tri_array[j][0]]
                v1 = indices[i][tri_array[j][1]]
                v2 = indices[i][tri_array[j][2]]
                for k in range(3):
                    tri.p0[k] = vertices[v0][k]
                    tri.p1[k] = vertices[v1][k]
                    tri.p2[k] = vertices[v2][k]
                    self.get_centroid(self.primitives,
                                      tri_index,
                                      self.centroids[tri_index])
                    self.get_bbox(self.primitives,
                                  tri_index, 
                                  &(self.bboxes[tri_index]))

        self.root = self._recursive_build(0, self.num_prim)

    cdef void _recursive_free(self, BVHNode* node) nogil:
        if node.end - node.begin > LEAF_SIZE:
            self._recursive_free(node.left)
            self._recursive_free(node.right)
        free(node)

    def __dealloc__(self):
        self._recursive_free(self.root)
        free(self.primitives)
        cdef np.int64_t i
        for i in range(self.num_prim):
            free(self.centroids[i])
        free(self.centroids)
        free(self.bboxes)
        free(self.field_data)
        free(self.vertices)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.int64_t _partition(self, np.int64_t begin, np.int64_t end,
                               np.int64_t ax, np.float64_t split) nogil:
        # this re-orders the triangle array so that all of the triangles 
        # to the left of mid have centroids less than or equal to "split"
        # along the direction "ax". All the triangles to the right of mid
        # will have centroids *greater* than "split" along "ax".
        cdef Triangle * triangles = <Triangle*> self.primitives

        cdef np.int64_t mid = begin
        while (begin != end):
            if self.centroids[mid][ax] > split:
                mid += 1
            elif self.centroids[begin][ax] > split:
                triangles[mid], triangles[begin] = \
                triangles[begin], triangles[mid]
                self.centroids[mid], self.centroids[begin] = \
                self.centroids[begin], self.centroids[mid]
                self.bboxes[mid], self.bboxes[begin] = \
                self.bboxes[begin], self.bboxes[mid]
                mid += 1
            begin += 1
        return mid
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _get_node_bbox(self, BVHNode* node, 
                             np.int64_t begin, np.int64_t end) nogil:
        cdef np.int64_t i, j
        cdef BBox box = self.bboxes[begin]
        for i in range(begin+1, end):
            for j in range(3):
                box.left_edge[j] = fmin(box.left_edge[j],
                                        self.bboxes[i].left_edge[j])
                box.right_edge[j] = fmax(box.right_edge[j], 
                                         self.bboxes[i].right_edge[j])
        node.bbox = box

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void intersect(self, Ray* ray) nogil:
        self._recursive_intersect(ray, self.root)
        
        if ray.elem_id < 0:
            return

        cdef np.float64_t[3] position
        cdef np.int64_t i
        for i in range(3):
            position[i] = ray.origin[i] + ray.t_far*ray.direction[i]
            
        cdef np.float64_t* vertex_ptr
        cdef np.float64_t* field_ptr
        vertex_ptr = self.vertices + ray.elem_id*self.num_verts_per_elem*3
        field_ptr = self.field_data + ray.elem_id*self.num_field_per_elem

        cdef np.float64_t[4] mapped_coord
        self.sampler.map_real_to_unit(mapped_coord, vertex_ptr, position)
        if self.num_field_per_elem == 1:
            ray.data_val = field_ptr[0]
        else:
            ray.data_val = self.sampler.sample_at_unit_point(mapped_coord,
                                                             field_ptr)
        ray.near_boundary = self.sampler.check_mesh_lines(mapped_coord)

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
                hit = self.get_intersect(self.primitives, i, ray)
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
