import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free

from yt.utilities.exceptions import YTElementTypeNotRecognized

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

cdef struct TriNode:
    np.uint64_t key
    np.int64_t elem
    np.int64_t tri[3]
    TriNode* next_node

cdef np.int64_t triangles_are_equal(np.int64_t tri1[3], np.int64_t tri2[3]) nogil:
    cdef np.int64_t found
    for i in range(3):
        found = False
        for j in range(3):
            if tri1[i] == tri2[j]:
                found = True
        if not found:
            return 0
    return 1
    
cdef np.uint64_t triangle_hash(np.int64_t tri[3]) nogil:
    # http://stackoverflow.com/questions/1536393/good-hash-function-for-permutations
    cdef np.uint64_t h = 1
    for i in range(3):
        h *= (1779033703 + 2*tri[i])
    return h / 2

# should be enough, consider dynamic resizing in the future
cdef np.int64_t TABLE_SIZE = 2**24

cdef class TriSet:
    '''

    This is a hash table data structure for rapidly identifying the exterior
    triangles in a polygon mesh. We loop over each triangle in each element and
    update the TriSet for each one. We keep only the triangles that appear once,
    as these make up the exterior of the mesh.

    '''

    cdef TriNode **table
    cdef np.uint64_t num_items
    
    def __cinit__(self):
        self.table = <TriNode**> malloc(TABLE_SIZE * sizeof(TriNode*))
        for i in range(TABLE_SIZE):
            self.table[i] = NULL
        self.num_items = 0
        
    def __dealloc__(self):
        cdef np.int64_t i
        cdef TriNode *node
        cdef TriNode *delete_node
        for i in range(TABLE_SIZE):
            node = self.table[i]
            while (node != NULL):
                delete_node = node
                node = node.next_node
                free(delete_node)
            self.table[i] = NULL
        free(self.table)
    
    def get_exterior_tris(self):
        '''

        Returns two numpy arrays, one storing the exterior triangle
        indices and the other storing the corresponding element ids.

        '''

        cdef np.ndarray[np.int64_t, ndim=2] tri_indices
        tri_indices = np.empty((self.num_items, 3), dtype="int64")
        cdef np.int64_t *tri_indices_ptr = <np.int64_t*> tri_indices.data

        cdef np.ndarray[np.int64_t, ndim=1] element_map
        element_map = np.empty(self.num_items, dtype="int64")
        cdef np.int64_t *elem_map_ptr = <np.int64_t*> element_map.data

        cdef TriNode* node
        cdef np.int64_t counter = 0
        cdef np.int64_t i, j
        for i in range(TABLE_SIZE):
            node = self.table[i]
            while node != NULL:
                for j in range(3):
                    tri_indices_ptr[3*counter + j] = node.tri[j]
                elem_map_ptr[counter] = node.elem
                counter += 1
                node = node.next_node
                
        return tri_indices, element_map

    cdef TriNode* _allocate_new_node(self,
                                     np.int64_t tri[3],
                                     np.uint64_t key,
                                     np.int64_t elem) nogil:
        cdef TriNode* new_node = <TriNode* > malloc(sizeof(TriNode))
        new_node.key = key
        new_node.elem = elem
        new_node.tri[0] = tri[0]
        new_node.tri[1] = tri[1]
        new_node.tri[2] = tri[2]
        new_node.next_node = NULL
        self.num_items += 1
        return new_node
        
    @cython.cdivision(True)
    cdef void update(self, np.int64_t tri[3], np.int64_t elem) nogil:
        cdef np.uint64_t key = triangle_hash(tri)
        cdef np.uint64_t index = key % TABLE_SIZE
        cdef TriNode *node = self.table[index]
        
        if node == NULL:
            self.table[index] = self._allocate_new_node(tri, key, elem)
            return

        if key == node.key and triangles_are_equal(node.tri, tri):
            # this triangle is already here, delete it
            self.table[index] = node.next_node
            free(node)
            self.num_items -= 1
            return

        elif node.next_node == NULL:
            node.next_node = self._allocate_new_node(tri, key, elem)
            return
    
        # walk through node list
        cdef TriNode* prev = node
        node = node.next_node
        while node != NULL:
            if key == node.key and triangles_are_equal(node.tri, tri):
                # this triangle is already here, delete it
                prev.next_node = node.next_node
                free(node)
                self.num_items -= 1
                return
            if node.next_node == NULL:
                # we have reached the end; add new node
                node.next_node = self._allocate_new_node(tri, key, elem)
                return
            prev = node
            node = node.next_node


cdef class MeshInfoHolder:
    cdef np.int64_t num_elem
    cdef np.int64_t num_tri
    cdef np.int64_t num_verts
    cdef np.int64_t VPE  # num verts per element
    cdef np.int64_t TPE  # num tris per element
    cdef int[MAX_NUM_TRI][3] tri_array
    
    def __cinit__(self, np.ndarray[np.int64_t, ndim=2] indices):
        self.num_elem = indices.shape[0]
        self.VPE = indices.shape[1]

        if (self.VPE == 8 or self.VPE == 20 or self.VPE == 27):
            self.TPE = HEX_NT
            self.tri_array = triangulate_hex
        elif self.VPE == 4:
            self.TPE = TETRA_NT
            self.tri_array = triangulate_tetra
        elif self.VPE == 6:
            self.TPE = WEDGE_NT
            self.tri_array = triangulate_wedge
        else:
            raise YTElementTypeNotRecognized(3, self.VPE)

        self.num_tri = self.TPE * self.num_elem
        self.num_verts = self.num_tri * 3
        

def cull_interior_triangles(np.ndarray[np.int64_t, ndim=2] indices):
    '''

    This is used to remove interior triangles from the mesh before rendering
    it on the GPU.

    '''

    cdef MeshInfoHolder m = MeshInfoHolder(indices)
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data

    cdef TriSet s = TriSet()
    cdef np.int64_t i, j, k
    cdef np.int64_t tri[3]
    for i in range(m.num_elem):
        for j in range(m.TPE):
            for k in range(3):
                tri[k] = indices_ptr[i*m.VPE + m.tri_array[j][k]]
            s.update(tri, i)

    return s.get_exterior_tris()
    

def get_vertex_data(np.ndarray[np.float64_t, ndim=2] coords,
                    np.ndarray[np.float64_t, ndim=2] data,
                    np.ndarray[np.int64_t, ndim=2] indices):

    '''

    This converts the data array from the shape (num_elem, conn_length)
    to (num_verts, ).

    '''

    cdef MeshInfoHolder m = MeshInfoHolder(indices)
    cdef np.int64_t num_verts = coords.shape[0]
    
    cdef np.ndarray[np.float32_t, ndim=1] vertex_data
    vertex_data = np.zeros(num_verts, dtype="float32")
    cdef np.float32_t *vertex_data_ptr = <np.float32_t*> vertex_data.data
    
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data
    cdef np.float64_t *data_ptr = <np.float64_t*> data.data
    
    cdef np.int64_t i, j
    for i in range(m.num_elem):
        for j in range(m.VPE):
            vertex_data_ptr[indices_ptr[i*m.VPE + j]] = data_ptr[i*m.VPE + j]

    return vertex_data


@cython.boundscheck(False)
def triangulate_mesh(np.ndarray[np.float64_t, ndim=2] coords,
                     np.ndarray data,
                     np.ndarray[np.int64_t, ndim=2] indices):
    '''

    This converts a mesh into a flattened triangle array suitable for
    rendering on the GPU.

    '''
    cdef np.ndarray[np.int64_t, ndim=2] exterior_tris
    cdef np.ndarray[np.int64_t, ndim=1] element_map
    exterior_tris, element_map = cull_interior_triangles(indices)

    cdef np.ndarray[np.float32_t, ndim=1] vertex_data
    if data.ndim == 2:
        vertex_data = get_vertex_data(coords, data, indices)
    else:
        vertex_data = np.empty(1, dtype=np.float32)
    cdef np.float32_t *vertex_data_ptr = <np.float32_t*> vertex_data.data

    cdef np.float64_t *data_ptr = <np.float64_t*> data.data
    cdef np.float64_t *coords_ptr = <np.float64_t*> coords.data
    cdef np.int64_t *exterior_tris_ptr = <np.int64_t*> exterior_tris.data

    cdef np.int64_t num_tri = exterior_tris.shape[0]
    cdef np.int64_t num_verts = 3 * num_tri
    cdef np.int64_t num_coords = 3 * num_verts

    cdef np.ndarray[np.int32_t, ndim=1] tri_indices
    tri_indices = np.empty(num_verts, dtype=np.int32)
    cdef np.int32_t *tri_indices_ptr = <np.int32_t*> tri_indices.data
    
    cdef np.ndarray[np.float32_t, ndim=1] tri_data
    tri_data = np.empty(num_verts, dtype=np.float32)
    cdef np.float32_t *tri_data_ptr = <np.float32_t*> tri_data.data

    cdef np.ndarray[np.float32_t, ndim=1] tri_coords
    tri_coords = np.empty(num_coords, dtype=np.float32)
    cdef np.float32_t *tri_coords_ptr = <np.float32_t*> tri_coords.data
    
    cdef np.int64_t vert_index, coord_index
    cdef np.int64_t i, j, k
    for i in range(num_tri):
        for j in range(3):
            vert_index = i*3 + j
            if data.ndim == 1:
                tri_data_ptr[vert_index] = data_ptr[element_map[i]]  # element data
            else:
                tri_data_ptr[vert_index] = vertex_data_ptr[exterior_tris_ptr[i*3 + j]]
            tri_indices_ptr[vert_index] = vert_index
            for k in range(3):
                coord_index = 3*exterior_tris_ptr[i*3 + j] + k
                tri_coords_ptr[vert_index*3 + k] = coords_ptr[coord_index]

    return tri_coords, tri_data, tri_indices


def triangulate_indices(np.ndarray[np.int64_t, ndim=2] indices):
    '''

    This is like triangulate_mesh, except it only considers the 
    connectivity information, instead of also copying the vertex 
    coordinates and the data values.

    '''

    cdef MeshInfoHolder m = MeshInfoHolder(indices)
        
    cdef np.ndarray[np.int64_t, ndim=2] tri_indices
    tri_indices = np.empty((m.num_tri, 3), dtype="int64")
    
    cdef np.int64_t *tri_indices_ptr = <np.int64_t*> tri_indices.data
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data

    cdef np.int64_t i, j, k, offset
    for i in range(m.num_elem):
        for j in range(m.TPE):
            for k in range(3):
                offset = m.tri_array[j][k]
                tri_indices_ptr[i*m.TPE*3 + 3*j + k] = indices_ptr[i*m.VPE + offset]
    return tri_indices
