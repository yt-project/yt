# This set of routines was originally written by Robert Thompson.
cimport numpy as np
import numpy as np
from libc.stdlib cimport free, malloc
from libc.math cimport sqrt
from fp_utils cimport iclip
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def smooth_particles(
        input_fields, output_grids,
        np.ndarray[np.float64_t, ndim=1] left_edge,
        np.ndarray[np.float64_t, ndim=1] right_edge,
        np.ndarray[np.float64_t, ndim=2] ppos,
        np.ndarray[np.float64_t, ndim=1] hsml):

    cdef np.int64_t ngas
    cdef np.float64_t dds[3], sdds[3], pos[3], idist[3], kern
    cdef int p, i, j, k, d, ind[3], ib0[3], ib1[3], dims[3]
    cdef int nf, half_len, skip, gi
    cdef np.float64_t dist
    cdef np.ndarray[np.float64_t, ndim=1] gas_arr
    cdef np.ndarray[np.float64_t, ndim=3] grid_arr

    ngas = input_fields[0].size
    for i in range(3):
        dims[i] = output_grids[0].shape[i]
        dds[i] = (right_edge[i] - left_edge[i])/dims[i]
        sdds[i] = dds[i] * dds[i]

    cdef np.float64_t *kernel_sum = \
        <np.float64_t *>malloc(ngas * sizeof(np.float64_t))
    cdef np.float64_t *pdist = \
        <np.float64_t *>malloc(dims[0]*dims[1]*dims[2]*
                               sizeof(np.float64_t))

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
            half_len = <int>(hsml[p]/dds[i]) + 1
            ib0[i] = ind[i] - half_len
            ib1[i] = ind[i] + half_len
            #pos[i] = ppos[p, i] - left_edge[i]
            #ind[i] = <int>(pos[i] / dds[i])
            #ib0[i] = <int>((pos[i] - hsml[i]) / dds[i]) - 1
            #ib1[i] = <int>((pos[i] + hsml[i]) / dds[i]) + 1
            if ib0[i] >= dims[i] or ib1[i] < 0:
                skip = 1
            ib0[i] = iclip(ib0[i], 0, dims[i] - 1)
            ib1[i] = iclip(ib1[i], 0, dims[i] - 1)
        if skip == 1: continue
        for i from ib0[0] <= i <= ib1[0]:
            idist[0] = (ind[0] - i) * (ind[0] - i) * sdds[0]
            for j from ib0[1] <= j <= ib1[1]:
                idist[1] = (ind[1] - j) * (ind[1] - j) * sdds[1] 
                for k from ib0[2] <= k <= ib1[2]:
                    idist[2] = (ind[2] - k) * (ind[2] - k) * sdds[2]
                    dist = idist[0] + idist[1] + idist[2]
                    dist = sqrt(dist) / hsml[p]
                    gi = ((i * dims[1] + j) * dims[2]) + k
                    pdist[gi] = sph_kernel(dist)
                    kernel_sum[p] += pdist[gi]
        for i from ib0[0] <= i <= ib1[0]:
            for j from ib0[1] <= j <= ib1[1]:
                for k from ib0[2] <= k <= ib1[2]:
                    gi = ((i * dims[1] + j) * dims[2]) + k
                    dist = pdist[gi] / kernel_sum[p]
                    for d in range(nf):
                        gdata[d][gi] += pdata[d][p] * dist
    free(kernel_sum)
    free(pdist)
    free(gdata)
    free(pdata)

##############################################
#Standard SPH kernel for use with the Grid method
cdef np.float64_t sph_kernel(np.float64_t x) nogil:
    cdef np.float64_t kernel
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel
