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

@cython.boundscheck(False)
@cython.wraparound(False)
def CICDeposit_3(np.ndarray[np.float64_t, ndim=1] posx,
                 np.ndarray[np.float64_t, ndim=1] posy,
                 np.ndarray[np.float64_t, ndim=1] posz,
                 np.ndarray[np.float32_t, ndim=1] mass,
                 np.int64_t npositions,
                 np.ndarray[np.float32_t, ndim=3] field,
                 np.ndarray[np.float64_t, ndim=1] leftEdge,
                 np.ndarray[np.int32_t, ndim=1] gridDimension,
                 np.float64_t cellSize):

    cdef int i1, j1, k1, n
    cdef double xpos, ypos, zpos
    cdef double fact, edge0, edge1, edge2
    cdef double le0, le1, le2
    cdef float dx, dy, dz, dx2, dy2, dz2

    edge0 = (<float> gridDimension[0]) - 0.5001
    edge1 = (<float> gridDimension[1]) - 0.5001
    edge2 = (<float> gridDimension[2]) - 0.5001
    fact = 1.0 / cellSize

    le0 = leftEdge[0]
    le1 = leftEdge[1]
    le2 = leftEdge[2]

    for n in range(npositions):

        # Compute the position of the central cell
        xpos = fmin(fmax((posx[n] - le0)*fact, 0.5001), edge0)
        ypos = fmin(fmax((posy[n] - le1)*fact, 0.5001), edge1)
        zpos = fmin(fmax((posz[n] - le2)*fact, 0.5001), edge2)

        i1  = <int> (xpos + 0.5)
        j1  = <int> (ypos + 0.5)
        k1  = <int> (zpos + 0.5)

        # Compute the weights
        dx = (<float> i1) + 0.5 - xpos
        dy = (<float> j1) + 0.5 - ypos
        dz = (<float> k1) + 0.5 - zpos
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
