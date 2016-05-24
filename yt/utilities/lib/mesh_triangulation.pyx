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
