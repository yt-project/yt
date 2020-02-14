cimport numpy as np

ctypedef fused integral:
    np.int8_t
    np.int32_t
    np.int64_t
    np.uint8_t
    np.uint32_t
    np.uint64_t

cdef np.uint64_t state_diagram[8][2][12]

state_diagram[0][0][:] = [ 1,  2,  0,  6, 11,  4,  5,  6, 10,  4,  7, 10]
state_diagram[0][1][:] = [ 0,  0,  0,  2,  4,  6,  4,  6,  2,  2,  4,  6]

state_diagram[1][0][:] = [ 2,  6,  9,  0, 11,  4,  7,  1,  3,  4,  2,  3]
state_diagram[1][1][:] = [ 1,  7,  3,  3,  3,  5,  7,  7,  5,  1,  5,  1]

state_diagram[2][0][:] = [ 3,  0, 10,  6,  0,  8,  5,  6,  1,  8, 11,  2]
state_diagram[2][1][:] = [ 3,  1,  7,  1,  5,  1,  3,  5,  3,  5,  7,  7]

state_diagram[3][0][:] = [ 2,  7,  9, 11,  7,  8,  3, 10,  1,  8,  2,  6]
state_diagram[3][1][:] = [ 2,  6,  4,  0,  2,  2,  0,  4,  4,  6,  6,  0]

state_diagram[4][0][:] = [ 4,  8,  1,  9,  5,  0,  1,  9, 10,  2,  7, 10]
state_diagram[4][1][:] = [ 7,  3,  1,  5,  7,  7,  5,  1,  1,  3,  3,  5]

state_diagram[5][0][:] = [ 5,  8,  1,  0,  9,  6,  1,  4,  3,  7,  5,  3]
state_diagram[5][1][:] = [ 6,  4,  2,  4,  0,  4,  6,  0,  6,  0,  2,  2]

state_diagram[6][0][:] = [ 3,  0, 11,  9,  0, 10, 11,  9,  5,  2,  8,  4]
state_diagram[6][1][:] = [ 4,  2,  6,  6,  6,  0,  2,  2,  0,  4,  0,  4]

state_diagram[7][0][:] = [ 5,  7, 11,  8,  7,  6, 11, 10,  9,  3,  5,  4]
state_diagram[7][1][:] = [ 5,  5,  5,  7,  1,  3,  1,  3,  7,  7,  1,  3]


cdef integral hilbert3d_single(integral x, integral y, integral z, int bit_length) nogil
cpdef np.ndarray hilbert3d(integral[:, :] x, int bit_length)