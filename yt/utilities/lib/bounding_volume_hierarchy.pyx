cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport fabs
from libc.stdlib cimport malloc, free
from cython.parallel import parallel, prange
from grid_traversal cimport ImageSampler, \
    ImageContainer

from yt.utilities.lib.primitives cimport \
    BBox, \
    Ray, \
    ray_bbox_intersect, \
    Triangle, \
    ray_triangle_intersect, \
    triangle_centroid, \
    triangle_bbox, \
    Patch, \
    ray_patch_intersect, \
    patch_centroid, \
    patch_bbox
from yt.utilities.lib.element_mappings cimport \
    ElementSampler, \
    Q1Sampler3D, \
    P1Sampler3D, \
    W1Sampler3D, \
    S2Sampler3D
from yt.utilities.lib.vec3_ops cimport L2_norm

cdef ElementSampler Q1Sampler = Q1Sampler3D()
cdef ElementSampler P1Sampler = P1Sampler3D()
cdef ElementSampler W1Sampler = W1Sampler3D()
cdef ElementSampler S2Sampler = S2Sampler3D()

cdef extern from "platform_dep.h" nogil:
    double fmax(double x, double y)
    double fmin(double x, double y)

# define some constants
cdef np.float64_t INF = np.inf
cdef np.int64_t   LEAF_SIZE = 16


cdef class BVH:
    '''

    This class implements a bounding volume hierarchy (BVH), a spatial acceleration
    structure for fast ray-tracing. A BVH is like a kd-tree, except that instead of 
    partitioning the *volume* of the parent to create the children, we partition the 
    primitives themselves into 'left' or 'right' sub-trees. The bounding volume for a
    node is then determined by computing the bounding volume of the primitives that
    belong to it. This allows us to quickly discard primitives that are not close 
    to intersecting a given ray.

    This class is currently used to provide software 3D rendering support for
    finite element datasets. For 1st-order meshes, every element of the mesh is
    triangulated, and this set of triangles forms the primitives that will be used
    for the ray-trace. The BVH can then quickly determine which element is hit by
    each ray associated with the image plane, and the appropriate interpolation can
    be performed to sample the finite element solution at that hit position.

    Currently, 2nd-order meshes are only supported for 20-node hexahedral elements.
    There, the primitive type is a bi-quadratic patch instead of a triangle, and
    each intersection involves computing a Netwon-Raphson solve.

    See yt/utilities/lib/primitives.pyx for the definitions of both of these primitive
    types.

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

        # We need to figure out what kind of elements we've been handed.
        if self.num_verts_per_elem == 8:
            self.num_prim_per_elem = HEX_NT
            self.tri_array = triangulate_hex
            self.sampler = Q1Sampler
        elif self.num_verts_per_elem == 6:
            self.num_prim_per_elem = WEDGE_NT
            self.tri_array = triangulate_wedge
            self.sampler = W1Sampler
        elif self.num_verts_per_elem == 4:
            self.num_prim_per_elem = TETRA_NT
            self.tri_array = triangulate_tetra
            self.sampler = P1Sampler
        elif self.num_verts_per_elem == 20:
            self.num_prim_per_elem = 6
            self.sampler = S2Sampler
        else:
            raise NotImplementedError("Could not determine element type for "
                                      "nverts = %d. " % self.num_verts_per_elem)
        self.num_prim = self.num_prim_per_elem*self.num_elem

        # allocate storage
        cdef np.int64_t v_size = self.num_verts_per_elem * self.num_elem * 3
        self.vertices = <np.float64_t*> malloc(v_size * sizeof(np.float64_t))
        cdef np.int64_t f_size = self.num_field_per_elem * self.num_elem
        self.field_data = <np.float64_t*> malloc(f_size * sizeof(np.float64_t))
        self.prim_ids = <np.int64_t*> malloc(self.num_prim * sizeof(np.int64_t))
        self.centroids = <np.float64_t**> malloc(self.num_prim * sizeof(np.float64_t*))
        cdef np.int64_t i
        for i in range(self.num_prim):
            self.centroids[i] = <np.float64_t*> malloc(3*sizeof(np.float64_t))
        self.bboxes = <BBox*> malloc(self.num_prim * sizeof(BBox))

        # create data buffers
        cdef np.int64_t j, k
        cdef np.int64_t field_offset, vertex_offset
        for i in range(self.num_elem):
            for j in range(self.num_verts_per_elem):
                vertex_offset = i*self.num_verts_per_elem*3 + j*3
                for k in range(3):
                    self.vertices[vertex_offset + k] = vertices[indices[i,j]][k]
            field_offset = i*self.num_field_per_elem
            for j in range(self.num_field_per_elem):
                self.field_data[field_offset + j] = field_data[i][j]                

        # set up primitives
        if self.num_verts_per_elem == 20:
            self.primitives = malloc(self.num_prim * sizeof(Patch))
            self.get_centroid = patch_centroid
            self.get_bbox = patch_bbox
            self.get_intersect = ray_patch_intersect
            self._set_up_patches(vertices, indices)
        else:
            self.primitives = malloc(self.num_prim * sizeof(Triangle))
            self.get_centroid = triangle_centroid
            self.get_bbox = triangle_bbox
            self.get_intersect = ray_triangle_intersect
            self._set_up_triangles(vertices, indices)
        
        self.root = self._recursive_build(0, self.num_prim)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _set_up_patches(self, np.float64_t[:, :] vertices,
                              np.int64_t[:, :] indices) nogil:
        cdef Patch* patch
        cdef np.int64_t i, j, k, ind, idim
        cdef np.int64_t offset, prim_index
        for i in range(self.num_elem):
            offset = self.num_prim_per_elem*i
            for j in range(self.num_prim_per_elem):  # for each face
                prim_index = offset + j
                patch = &( <Patch*> self.primitives)[prim_index]
                self.prim_ids[prim_index] = prim_index
                patch.elem_id = i
                for k in range(8):  # for each vertex
                    ind = hex20_faces[j][k]
                    for idim in range(3):  # for each spatial dimension (yikes)
                        patch.v[k][idim] = vertices[indices[i, ind]][idim]
                self.get_centroid(self.primitives,
                                  prim_index,
                                  self.centroids[prim_index])
                self.get_bbox(self.primitives,
                              prim_index,
                              &(self.bboxes[prim_index]))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void _set_up_triangles(self, np.float64_t[:, :] vertices,
                                np.int64_t[:, :] indices) nogil:
        # fill our array of primitives
        cdef np.int64_t offset, tri_index
        cdef np.int64_t v0, v1, v2
        cdef Triangle* tri
        cdef np.int64_t i, j, k
        for i in range(self.num_elem):
            offset = self.num_prim_per_elem*i
            for j in range(self.num_prim_per_elem):
                tri_index = offset + j
                self.prim_ids[tri_index] = tri_index
                tri = &(<Triangle*> self.primitives)[tri_index]
                tri.elem_id = i
                v0 = indices[i][self.tri_array[j][0]]
                v1 = indices[i][self.tri_array[j][1]]
                v2 = indices[i][self.tri_array[j][2]]
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

    cdef void _recursive_free(self, BVHNode* node) nogil:
        if node.end - node.begin > LEAF_SIZE:
            self._recursive_free(node.left)
            self._recursive_free(node.right)
        free(node)

    def __dealloc__(self):
        self._recursive_free(self.root)
        free(self.primitives)
        free(self.prim_ids)
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
        # this re-orders the primitive array so that all of the primitives
        # to the left of mid have centroids less than or equal to "split"
        # along the direction "ax". All the primitives to the right of mid
        # will have centroids *greater* than "split" along "ax".
        cdef np.int64_t mid = begin
        while (begin != end):
            if self.centroids[mid][ax] > split:
                mid += 1
            elif self.centroids[begin][ax] > split:
                self.prim_ids[mid], self.prim_ids[begin] = \
                self.prim_ids[begin], self.prim_ids[mid]
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
        if (node.end - node.begin) <= LEAF_SIZE:
            for i in range(node.begin, node.end):
                hit = self.get_intersect(self.primitives, self.prim_ids[i], ray)
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

cdef class BVHMeshSampler(ImageSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self,
                 BVH bvh,
                 int num_threads = 0):
        '''

        This function is supposed to cast the rays and return the
        image.

        '''

        cdef int vi, vj, i, j
        cdef ImageContainer *im = self.image
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.int64_t nx, ny, size
        cdef np.float64_t width[3]
        for i in range(3):
            width[i] = self.width[i]
        nx = im.nv[0]
        ny = im.nv[1]
        size = nx * ny
        cdef Ray* ray
        with nogil, parallel():
            ray = <Ray *> malloc(sizeof(Ray))
            v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            for j in prange(size):
                vj = j % ny
                vi = (j - vj) / ny
                vj = vj
                self.vector_function(im, vi, vj, width, v_dir, v_pos)
                for i in range(3):
                    ray.origin[i] = v_pos[i]
                    ray.direction[i] = v_dir[i]
                    ray.inv_dir[i] = 1.0 / v_dir[i]
                ray.t_far = 1e37
                ray.t_near = 0.0
                ray.data_val = 0
                ray.elem_id = -1
                bvh.intersect(ray)
                im.image[vi, vj, 0] = ray.data_val
                im.image_used[vi, vj] = ray.elem_id
                im.mesh_lines[vi, vj] = ray.near_boundary
                im.zbuffer[vi, vj] = ray.t_far
            free(v_pos)
            free(v_dir)
            free(ray)
