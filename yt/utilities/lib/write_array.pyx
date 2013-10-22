"""
Faster, cythonized file IO



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
def write_3D_array(np.ndarray[DTYPE_t, ndim=3] data, fhandle):
    assert data.dtype == DTYPE
    cdef int Nx, Ny, Nz
    Nx = data.shape[0]
    Ny = data.shape[1]
    Nz = data.shape[2]
    cdef unsigned int i, j, k

    for i in np.arange(Nz):
        for j in np.arange(Ny):
            for k in np.arange(Nx):
                fhandle.write(str(data[k, j, i]) + '\n')

@cython.boundscheck(False)
def write_3D_vector_array(np.ndarray[DTYPE_t, ndim=3] data_x, 
                          np.ndarray[DTYPE_t, ndim=3] data_y,
                          np.ndarray[DTYPE_t, ndim=3] data_z,
                          fhandle):

    assert data_x.dtype == DTYPE
    cdef int Nx, Ny, Nz
    Nx = data_x.shape[0]
    Ny = data_x.shape[1]
    Nz = data_x.shape[2]
    cdef unsigned int i, j, k

    for i in np.arange(Nz):
        for j in np.arange(Ny):
            for k in np.arange(Nx):
                fx = data_x[k, j, i]
                fy = data_y[k, j, i]
                fz = data_z[k, j, i]
                fhandle.write('{}    {}    {} \n'.format(fx, fy, fz))
