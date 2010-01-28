"""
Simple interpolators

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
def UnilinearlyInterpolate(np.ndarray[np.float64_t, ndim=1] table,
                           np.ndarray[np.float64_t, ndim=1] x_vals,
                           np.ndarray[np.float64_t, ndim=1] x_bins,
                           np.ndarray[np.int32_t, ndim=1] x_is,
                           np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef int i, x_i, y_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp = (x - x_bins[x_i]) * dx_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        output[i]  = table[x_i  ] * (xm) \
                   + table[x_i+1] * (xp)

@cython.boundscheck(False)
def BilinearlyInterpolate(np.ndarray[np.float64_t, ndim=2] table,
                          np.ndarray[np.float64_t, ndim=1] x_vals,
                          np.ndarray[np.float64_t, ndim=1] y_vals,
                          np.ndarray[np.float64_t, ndim=1] x_bins,
                          np.ndarray[np.float64_t, ndim=1] y_bins,
                          np.ndarray[np.int32_t, ndim=1] x_is,
                          np.ndarray[np.int32_t, ndim=1] y_is,
                          np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef double y, yp, ym
    cdef double dx_inv, dy_inv
    cdef int i, x_i, y_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        y_i = y_is[i]
        x = x_vals[i]
        y = y_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        dy_inv = 1.0 / (y_bins[y_i+1] - y_bins[y_i])
        xp = (x - x_bins[x_i]) * dx_inv
        yp = (y - y_bins[y_i]) * dy_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        ym = (y_bins[y_i+1] - y) * dy_inv
        output[i]  = table[x_i  , y_i  ] * (xm*ym) \
                   + table[x_i+1, y_i  ] * (xp*ym) \
                   + table[x_i  , y_i+1] * (xm*yp) \
                   + table[x_i+1, y_i+1] * (xp*yp)

@cython.boundscheck(False)
def TrilinearlyInterpolate(np.ndarray[np.float64_t, ndim=3] table,
                           np.ndarray[np.float64_t, ndim=1] x_vals,
                           np.ndarray[np.float64_t, ndim=1] y_vals,
                           np.ndarray[np.float64_t, ndim=1] z_vals,
                           np.ndarray[np.float64_t, ndim=1] x_bins,
                           np.ndarray[np.float64_t, ndim=1] y_bins,
                           np.ndarray[np.float64_t, ndim=1] z_bins,
                           np.ndarray[np.int_t, ndim=1] x_is,
                           np.ndarray[np.int_t, ndim=1] y_is,
                           np.ndarray[np.int_t, ndim=1] z_is,
                           np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef double y, yp, ym
    cdef double z, zp, zm
    cdef double dx_inv, dy_inv, dz_inv
    cdef int i, x_i, y_i, z_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        y_i = y_is[i]
        z_i = z_is[i]
        x = x_vals[i]
        y = y_vals[i]
        z = z_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        dy_inv = 1.0 / (y_bins[y_i+1] - y_bins[y_i])
        dz_inv = 1.0 / (z_bins[z_i+1] - z_bins[z_i])
        xp = (x - x_bins[x_i]) * dx_inv
        yp = (y - y_bins[y_i]) * dy_inv
        zp = (z - z_bins[z_i]) * dz_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        ym = (y_bins[y_i+1] - y) * dy_inv
        zm = (z_bins[z_i+1] - z) * dz_inv
        output[i]  = table[x_i  ,y_i  ,z_i  ] * (xm*ym*zm) \
                   + table[x_i+1,y_i  ,z_i  ] * (xp*ym*zm) \
                   + table[x_i  ,y_i+1,z_i  ] * (xm*yp*zm) \
                   + table[x_i  ,y_i  ,z_i+1] * (xm*ym*zp) \
                   + table[x_i+1,y_i  ,z_i+1] * (xp*ym*zp) \
                   + table[x_i  ,y_i+1,z_i+1] * (xm*yp*zp) \
                   + table[x_i+1,y_i+1,z_i  ] * (xp*yp*zm) \
                   + table[x_i+1,y_i+1,z_i+1] * (xp*yp*zp)
