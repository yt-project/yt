import numpy as np
cimport numpy as np
cimport cython

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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def triangulate_indices(np.ndarray[np.int64_t, ndim=2] indices):
    cdef np.int64_t num_elem = indices.shape[0]
    cdef np.int64_t VPE = indices.shape[1]  # num verts per element
    cdef np.int64_t TPE  # num triangles per element
    cdef int[MAX_NUM_TRI][3] tri_array
    
    if (VPE == 8 or VPE == 20 or VPE == 27):
        TPE = HEX_NT
        tri_array = triangulate_hex
    elif VPE == 4:
        TPE = TETRA_NT
        tri_array = triangulate_tetra
    elif VPE == 6:
        TPE = WEDGE_NT
        tri_array = triangulate_wedge
    else:
        raise YTElementTypeNotRecognized(3, VPE)

    cdef np.int64_t num_tri = TPE * num_elem
        
    cdef np.ndarray[np.int64_t, ndim=2] tri_indices
    tri_indices = np.empty((num_tri, 3), dtype="int64")
    
    cdef np.int64_t *tri_indices_ptr = <np.int64_t*> tri_indices.data
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data

    cdef np.int64_t i, j, k, offset
    for i in range(num_elem):
        for j in range(TPE):
            for k in range(3):
                offset = tri_array[j][k]
                tri_indices_ptr[i*TPE*3 + 3*j + k] = indices_ptr[i*VPE + offset]
    return tri_indices

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def triangulate_vertex_data(np.ndarray[np.float64_t, ndim=2] coords,
                            np.ndarray[np.float64_t, ndim=2] data,
                            np.ndarray[np.int64_t, ndim=2] indices):

    cdef np.int64_t num_elem = indices.shape[0]
    cdef np.int64_t VPE = indices.shape[1]  # num verts per element
    cdef np.int64_t TPE  # num triangles per element
    cdef int[MAX_NUM_TRI][3] tri_array
    
    if (VPE == 8 or VPE == 20 or VPE == 27):
        TPE = HEX_NT
        tri_array = triangulate_hex
    elif VPE == 4:
        TPE = TETRA_NT
        tri_array = triangulate_tetra
    elif VPE == 6:
        TPE = WEDGE_NT
        tri_array = triangulate_wedge
    else:
        raise YTElementTypeNotRecognized(3, VPE)

    cdef np.int64_t num_tri = TPE * num_elem
    cdef np.int64_t num_verts = coords.shape[0]
    
    cdef np.ndarray[np.int32_t, ndim=1] tri_indices
    tri_indices = np.empty(num_tri*3, dtype="int32")
    cdef np.int32_t *tri_indices_ptr = <np.int32_t*> tri_indices.data
    
    cdef np.ndarray[np.float32_t, ndim=1] tri_coords
    tri_coords = np.empty(num_verts*3, dtype=np.float32)
    cdef np.float32_t *tri_coords_ptr = <np.float32_t*> tri_coords.data
    
    cdef np.ndarray[np.float32_t, ndim=1] tri_data
    tri_data = np.zeros(num_verts, dtype="float32")
    cdef np.float32_t *tri_data_ptr = <np.float32_t*> tri_data.data
    
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data
    cdef np.float64_t *coords_ptr = <np.float64_t*> coords.data
    cdef np.float64_t *data_ptr = <np.float64_t*> data.data
    
    cdef np.int64_t i, j, k, offset
    for i in range(num_elem):
        for j in range(TPE):
            for k in range(3):
                offset = tri_array[j][k]
                tri_indices_ptr[i*TPE*3 + 3*j + k] = indices_ptr[i*VPE + offset]
        for j in range(VPE):
            tri_data_ptr[indices_ptr[i*VPE + j]] = data_ptr[i*VPE + j]
    for i in range(num_verts):
        for j in range(3):
            tri_coords_ptr[i*3 + j] = coords_ptr[i*3 + j]

    return tri_coords, tri_data, tri_indices

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def triangulate_element_data(np.ndarray[np.float64_t, ndim=2] coords,
                             np.ndarray[np.float64_t, ndim=1] data,
                             np.ndarray[np.int64_t, ndim=2] indices):

    cdef np.int64_t num_elem = indices.shape[0]
    cdef np.int64_t VPE = indices.shape[1]  # num verts per element
    cdef np.int64_t TPE  # num triangles per element
    cdef int[MAX_NUM_TRI][3] tri_array
    
    if (VPE == 8 or VPE == 20 or VPE == 27):
        TPE = HEX_NT
        tri_array = triangulate_hex
    elif VPE == 4:
        TPE = TETRA_NT
        tri_array = triangulate_tetra
    elif VPE == 6:
        TPE = WEDGE_NT
        tri_array = triangulate_wedge
    else:
        raise YTElementTypeNotRecognized(3, VPE)

    cdef np.int64_t num_tri = TPE * num_elem
    cdef np.int64_t num_verts = num_tri*3
    
    cdef np.ndarray[np.int32_t, ndim=1] tri_indices
    tri_indices = np.empty(num_verts, dtype="int32")
    cdef np.int32_t *tri_indices_ptr = <np.int32_t*> tri_indices.data
    
    cdef np.ndarray[np.float32_t, ndim=1] tri_data
    tri_data = np.empty(num_verts, dtype=np.float32)
    cdef np.float32_t *tri_data_ptr = <np.float32_t*> tri_data.data

    cdef np.ndarray[np.float32_t, ndim=1] tri_coords
    tri_coords = np.empty(num_verts*3, dtype=np.float32)
    cdef np.float32_t *tri_coords_ptr = <np.float32_t*> tri_coords.data
    
    cdef np.float64_t *data_ptr = <np.float64_t*> data.data
    cdef np.float64_t *coords_ptr = <np.float64_t*> coords.data
    cdef np.int64_t *indices_ptr = <np.int64_t*> indices.data

    cdef np.int64_t i, j, k, l, 
    cdef np.int64_t coord_index, vert_index
    for i in range(num_elem):
        for j in range(TPE):
            for k in range(3):
                coord_index = indices_ptr[i*VPE + tri_array[j][k]]*3
                vert_index = i*TPE*3+j*3+k
                tri_data_ptr[vert_index] = data_ptr[i]
                for l in range(3):
                    tri_coords_ptr[vert_index*3 + l] = coords_ptr[coord_index + l]
                tri_indices_ptr[vert_index] = vert_index

    return tri_coords, tri_data, tri_indices
