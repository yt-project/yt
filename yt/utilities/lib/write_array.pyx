"""
Faster, cythonized file IO

Author: Andrew Myers
Affiliation: UCB
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Andrew Myers.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
