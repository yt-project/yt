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
    cdef np.int64_t num_verts = indices.shape[1]
    cdef int[MAX_NUM_TRI][3] tri_array
    cdef np.int64_t tri_per_elem

    if (num_verts == 8 or num_verts == 20 or num_verts == 27):
        tri_per_elem = HEX_NT
        tri_array = triangulate_hex
    elif num_verts == 4:
        tri_per_elem = TETRA_NT
        tri_array = triangulate_tetra
    elif num_verts == 6:
        tri_per_elem = WEDGE_NT
        tri_array = triangulate_wedge
    else:
        raise YTElementTypeNotRecognized(3, num_verts)

    cdef np.ndarray[np.int64_t, ndim=2] tri_indices
    tri_indices = np.empty((tri_per_elem * num_elem, 3), dtype="int64")

    cdef np.int64_t i, j
    for i in range(num_elem):
        for j in range(tri_per_elem):
            for k in range(3):
                tri_indices[i*tri_per_elem + j][k] = indices[i][tri_array[j][k]]
    return tri_indices
