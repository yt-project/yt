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

@cython.wraparound(False)
@cython.boundscheck(False)
def VoxelTraversal(np.ndarray[np.int_t, ndim=3] grid_mask,
                   np.ndarray[np.float64_t, ndim=3] grid_t,
                   np.ndarray[np.float64_t, ndim=1] left_edge,
                   np.ndarray[np.float64_t, ndim=1] right_edge,
                   np.ndarray[np.float64_t, ndim=1] dx,
                   np.ndarray[np.float64_t, ndim=1] u,
                   np.ndarray[np.float64_t, ndim=1] v):
    # We're roughly following Amanatides & Woo
    # Find the first place the ray hits the grid on its path
    # Do left edge then right edge in each dim
    cdef int i, x, y
    cdef np.float64_t tl, tr, intersect_t, enter_t, exit_t, dt_tolerance
    cdef np.ndarray[np.int64_t,   ndim=1] step = np.empty((3,), dtype=np.int64)
    cdef np.ndarray[np.int64_t,   ndim=1] cur_ind = np.empty((3,), dtype=np.int64)
    cdef np.ndarray[np.float64_t, ndim=1] tdelta = np.empty((3,), dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] tmax = np.empty((3,), dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] intersect = np.empty((3,), dtype=np.float64)
    intersect_t = 1
    dt_tolerance = 1e-6
    # recall p = v * t + u
    #  where p is position, v is our vector, u is the start point
    for i in range(3):
        # As long as we're iterating, set some other stuff, too
        if(v[i] < 0): step[i] = -1
        else: step[i] = 1
        x = (i+1)%3
        y = (i+2)%3
        tl = (left_edge[i] - u[i])/v[i]
        tr = (right_edge[i] - u[i])/v[i]
        if (left_edge[x] <= (u[x] + tl*v[x]) <= right_edge[x]) and \
           (left_edge[y] <= (u[y] + tl*v[y]) <= right_edge[y]) and \
           (0.0 <= tl < intersect_t):
            intersect_t = tl
        if (left_edge[x] <= (u[x] + tr*v[x]) <= right_edge[x]) and \
           (left_edge[y] <= (u[y] + tr*v[y]) <= right_edge[y]) and \
           (0.0 <= tr < intersect_t):
            intersect_t = tr
    # if fully enclosed
    if (left_edge[0] <= u[0] <= right_edge[0]) and \
       (left_edge[1] <= u[1] <= right_edge[1]) and \
       (left_edge[2] <= u[2] <= right_edge[2]):
        intersect_t = 0.0
    if not (0 <= intersect_t <= 1): return
    # Now get the indices of the intersection
    intersect = u + intersect_t * v
    cdef int ncells = 0
    for i in range(3):
        cur_ind[i] = np.floor((intersect[i] + 1e-8*dx[i] - left_edge[i])/dx[i])
        tmax[i] = (((cur_ind[i]+step[i])*dx[i])+left_edge[i]-u[i])/v[i]
        if cur_ind[i] == grid_mask.shape[i] and step[i] < 0:
            cur_ind[i] = grid_mask.shape[i] - 1
        if step[i] > 0: tmax[i] = (((cur_ind[i]+1)*dx[i])+left_edge[i]-u[i])/v[i]
        if step[i] < 0: tmax[i] = (((cur_ind[i]+0)*dx[i])+left_edge[i]-u[i])/v[i]
        tdelta[i] = np.abs((dx[i]/v[i]))
    # The variable intersect contains the point we first pierce the grid
    enter_t = intersect_t
    while 1:
        if (not (0 <= cur_ind[0] < grid_mask.shape[0])) or \
           (not (0 <= cur_ind[1] < grid_mask.shape[1])) or \
           (not (0 <= cur_ind[2] < grid_mask.shape[2])):
            break
        # Note that we are calculating t on the fly, but we get *negative* t
        # values from what they should be.
        # If we've reached t = 1, we are done.
        grid_mask[cur_ind[0], cur_ind[1], cur_ind[2]] = 1
        if (tmax[0] > 1.0) and (tmax[1] > 1.0) and (tmax[2] > 1.0):
            grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = 1.0 - enter_t
            break
        ncells += 1
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[0] - enter_t
                enter_t = tmax[0]
                tmax[0] += tdelta[0]
                cur_ind[0] += step[0]
            else:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
        else:
            if tmax[1] < tmax[2]:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[1] - enter_t
                enter_t = tmax[1]
                tmax[1] += tdelta[1]
                cur_ind[1] += step[1]
            else:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
    return
