# This set of routines was originally written by Robert Thompson.
cimport numpy as np
import numpy as np
from libc.stdlib cimport free, malloc
from libc.math cimport sqrt
cimport cython

#@cython.boundscheck(False)
#@cython.wraparound(False)
@cython.cdivision(True)
def grid(input_fields, output_grids,
         np.ndarray[np.float64_t, ndim=1] left_edge,
         np.ndarray[np.float64_t, ndim=1] right_edge,
         np.ndarray[np.float64_t, ndim=2] ppos,
         np.ndarray[np.float64_t, ndim=1] hsml):

    cdef np.int64_t ngas
    cdef np.float64_t dds[3], pos[3], idist[3], kern
    cdef int p, i, j, k, d, ind[3], ib0[3], ib1[3], dims[3]
    cdef int nf, half_len, skip
    cdef np.float64_t dist
    cdef np.ndarray[np.float64_t, ndim=1] gas_arr
    cdef np.ndarray[np.float64_t, ndim=3] grid_arr

    ngas = input_fields[0].size
    for i in range(3):
        dims[i] = output_grids[0].shape[i]
        dds[i] = (right_edge[i] - left_edge[i])/dims[i]

    cdef np.float64_t *kernel_sum = \
        <np.float64_t *>malloc(ngas * sizeof(np.float64_t))

    nf = len(input_fields)
    assert(nf == len(output_grids))

    cdef np.float64_t **pdata = <np.float64_t**> \
            malloc(sizeof(np.float64_t *) * nf)
    cdef np.float64_t **gdata = <np.float64_t**> \
            malloc(sizeof(np.float64_t *) * nf)

    for i in range(nf):
        gas_arr = input_fields[i]
        grid_arr = output_grids[i]
        pdata[i] = <np.float64_t*> gas_arr.data
        gdata[i] = <np.float64_t*> grid_arr.data

    for p in range(ngas):
        kernel_sum[p] = 0.0
        skip = 0
        for i in range(3):
            pos[i] = ppos[p, i]
            ind[i] = <int>((pos[i] - left_edge[i]) / dds[i])
            half_len = <int>(hsml[p]/dds[i])
            ib0[i] = ind[i] - half_len
            ib1[i] = ind[i] + half_len
            if ib0[i] >= dims[i] or ib1[i] < 0:
                skip = 1
        if skip == 1: continue

        for i from ib0[0] <= i <= ib1[0]:
            idist[0] = (ind[0] - i) * (ind[0] - i) * dds[0]
            for j from ib0[1] <= j <= ib1[1]:
                idist[1] = (ind[1] - j) * (ind[1] - j) * dds[1] 
                idist[1] += idist[0]
                for k from ib0[2] <= k <= ib1[2]:
                    idist[2] = (ind[2] - k) * (ind[2] - k) * dds[2]
                    idist[2] += idist[1]
                    dist = sqrt(idist[2]) / hsml[p]
                    kernel_sum[p] += sph_kernel(dist)

        for i from ib0[0] <= i <= ib1[0]:
            idist[0] = (ind[0] - i) * (ind[0] - i) * dds[0]
            for j from ib0[1] <= j <= ib1[1]:
                idist[1] = (ind[1] - j) * (ind[1] - j) * dds[1] 
                idist[1] += idist[0]
                for k from ib0[2] <= k <= ib1[2]:
                    idist[2] = (ind[2] - k) * (ind[2] - k) * dds[2]
                    idist[2] += idist[1]
                    dist = sqrt(idist[2]) / hsml[p]
                    kern = sph_kernel(dist)
                    gi = ((i * dims[1] + j) * dims[2]) + k
                    for d in range(nf):
                        gdata[d][gi] += pdata[d][p] * kern / kernel_sum[p]

    free(kernel_sum)
    free(gdata)
    free(pdata)

##############################################
#Standard SPH kernel for use with the Grid method
cdef float sph_kernel(float x) nogil:
    cdef float kernel
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel
