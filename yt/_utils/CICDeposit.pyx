"""
Simle integrators for the radiative transfer equation

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: CASA/University of Colorado
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

def CICDeposit_3(np.ndarray[np.float64_t, ndim=1] posx,
                 np.ndarray[np.float64_t, ndim=1] posy,
                 np.ndarray[np.float64_t, ndim=1] posz,
                 np.ndarray[np.float32_t, ndim=1] mass,
                 np.int64_t npositions,
                 np.ndarray[np.float32_t, ndim=3] field,
                 np.ndarray[np.float64_t, ndim=1] leftEdge,
                 np.ndarray[np.int32_t, ndim=1] gridDimension,
                 np.float32_t cellSize):

    cdef int i1, j1, k1, n
    cdef np.float64_t xpos, ypos, zpos
    cdef np.float64_t fact, edge0, edge1, edge2
    cdef np.float32_t dx, dy, dz, dx2, dy2, dz2

    edge0 = np.float64(gridDimension[0]) - 0.5001
    edge1 = np.float64(gridDimension[1]) - 0.5001
    edge2 = np.float64(gridDimension[2]) - 0.5001
    fact = 1.0 / cellSize

    for n in xrange(npositions):

        # Compute the position of the central cell
        xpos = min([max([(posx[n] - leftEdge[0])*fact, 0.5001]), edge0])
        ypos = min([max([(posy[n] - leftEdge[1])*fact, 0.5001]), edge1])
        zpos = min([max([(posz[n] - leftEdge[2])*fact, 0.5001]), edge2])

        i1  = np.int32(xpos + 0.5)
        j1  = np.int32(ypos + 0.5)
        k1  = np.int32(zpos + 0.5)

        # Compute the weights
        dx = np.float(i1) + 0.5 - xpos
        dy = np.float(j1) + 0.5 - ypos
        dz = np.float(k1) + 0.5 - zpos
        dx2 =  1.0 - dx
        dy2 =  1.0 - dy
        dz2 =  1.0 - dz

        # Interpolate from field into sumfield
        field[i1-1,j1-1,k1-1] += mass[n] * dx  * dy  * dz
        field[i1  ,j1-1,k1-1] += mass[n] * dx2 * dy  * dz
        field[i1-1,j1  ,k1-1] += mass[n] * dx  * dy2 * dz
        field[i1  ,j1  ,k1-1] += mass[n] * dx2 * dy2 * dz
        field[i1-1,j1-1,k1  ] += mass[n] * dx  * dy  * dz2
        field[i1  ,j1-1,k1  ] += mass[n] * dx2 * dy  * dz2
        field[i1-1,j1  ,k1  ] += mass[n] * dx  * dy2 * dz2
        field[i1  ,j1  ,k1  ] += mass[n] * dx2 * dy2 * dz2
