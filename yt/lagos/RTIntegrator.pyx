"""
Simle integrators for the radiative transfer equation

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
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

@cython.boundscheck(False)
def Transfer3D(np.ndarray[np.float_t, ndim=2] i_s,
               np.ndarray[np.float_t, ndim=3] o_s,
               np.ndarray[np.float_t, ndim=3] e,
               np.ndarray[np.float_t, ndim=3] a,
               int imin, int imax, int jmin, int jmax,
               int kmin, int kmax, int istride, int jstride,
               float dx):
    """
    This function accepts an incoming slab (*i_s*), a buffer
    for an outgoing set of values at every point in the grid (*o_s*),
    an emission array (*e*), an absorption array (*a*), and dimensions of
    the grid (*imin*, *imax*, *jmin*, *jmax*, *kmin*, *kmax*) as well
    as strides in the *i* and *j* directions, and a *dx* of the grid being
    integrated.
    """
    cdef int i, ii
    cdef int j, jj
    cdef int k, kk
    cdef float temp
    for i in range((imax-imin)*istride):
        ii = i + imin*istride
        for j in range((jmax-jmin)*jstride):
            jj = j + jmin*jstride
            temp = i_s[ii,jj]
            for k in range(kmax-kmin):
                o_s[i,j,k] = temp + dx*(e[i,j,k] - temp*a[i,j,k])
                temp = o_s[i,j,k]
            i_s[ii,jj] = temp

@cython.boundscheck(False)
def Transfer1D(float i_s,
               np.ndarray[np.float_t, ndim=1] o_s,
               np.ndarray[np.float_t, ndim=1] e,
               np.ndarray[np.float_t, ndim=1] a,
               np.ndarray[np.float_t, ndim=1] dx,
               int imin, int imax):
    cdef int i
    for i in range(imin, imax):
        o_s[i] = i_s + dx[i]*(e[i] - i_s*a[i])
        i_s = o_s[i]
    return i_s
