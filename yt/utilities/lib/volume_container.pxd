"""
A volume container




"""


cimport numpy as np


cdef struct VolumeContainer:
    #-----------------------------------------------------------------------------
    # Encapsulates a volume container used for volume rendering.
    #
    #    n_fields       int              : The number of fields available to the volume renderer.
    #    data           np.float64_t**   : The data within the volume container.
    #    mask           np.uint8_t*      : The mask of the volume container. It has dimensions one fewer in each direction than data.
    #    left_edge      np.float64_t[3]  : The left edge of the volume container's bounding box.
    #    right_edge     np.float64_t[3]  : The right edge of the volume container's bounding box.
    #    np.float64_t   dds[3]           : The delta dimensions, such that dds[0] = ddx, dds[1] = ddy, dds[2] = ddz.
    #    np.float64_t   idds[3]          : The inverse delta dimensions. Same as dds pattern, but the inverse. i.e. idds[0] = inv_ddx.
    #    dims           int[3]           : The dimensions of the volume container. dims[0] = x, dims[1] = y, dims[2] = z.
    #-----------------------------------------------------------------------------
    int n_fields
    np.float64_t **data
    np.uint8_t *mask
    np.float64_t left_edge[3]
    np.float64_t right_edge[3]
    np.float64_t dds[3]
    np.float64_t idds[3]
    int dims[3]

cdef inline int vc_index(VolumeContainer *vc, int i, int j, int k):
    #-----------------------------------------------------------------------------
    # vc_index(VolumeContainer *vc, int i, int j, int k)
    #    vc   VolumeContainer* : Pointer to the volume container to be indexed.
    #    i    int              : The first dimension coordinate.
    #    j    int              : The second dimension coordinate.
    #    k    int              : The third dimension coordinates.
    #
    # Returns the 3-dimensional index in the volume container given coordinates i, j, k.
    # This is used for 3-dimensional access in a flat container using C-ordering.
    # This is calculated by:
    #       vc_index = i * vc.dim[1] * vc.dims[2] + j * vc.dims[2] + k
    # and then simplified (as shown below) by combining one multiplication operation.
    #
    # 2-dimensional example:
    #       A 4 x 3 array may be represented as:
    #                                      a = [0,  1,  2,  3,
    #                                           4,  5,  6,  7,
    #                                           8,  9,  10, 11]
    #       or similarly, in a flat container in row-successive order as:
    #                          a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    #
    # To access the coordinate at (1,1) in the flat container:
    #                         i * dims[1] + j
    #                       = 1 *   3     + 1
    #                       = 4
    # The 3-dimensional case (presented below) is similar.
    #-----------------------------------------------------------------------------
    return (i*vc.dims[1]+j)*vc.dims[2]+k

cdef inline int vc_pos_index(VolumeContainer *vc, np.float64_t *spos):
    cdef int index[3]
    cdef int i
    for i in range(3):
        index[i] = <int> ((spos[i] - vc.left_edge[i]) * vc.idds[i])
    return vc_index(vc, index[0], index[1], index[2])
