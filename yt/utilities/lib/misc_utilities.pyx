"""
Simple utilities that don't fit anywhere else

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
cimport libc.math as math

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_profile1d(np.ndarray[np.int64_t, ndim=1] bins_x,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=1] bsource,
                  np.ndarray[np.float64_t, ndim=1] wresult,
                  np.ndarray[np.float64_t, ndim=1] bresult,
                  np.ndarray[np.float64_t, ndim=1] mresult,
                  np.ndarray[np.float64_t, ndim=1] qresult,
                  np.ndarray[np.float64_t, ndim=1] used):
    cdef int n, bin
    cdef np.float64_t wval, bval
    for n in range(bins_x.shape[0]):
        bin = bins_x[n]
        bval = bsource[n]
        wval = wsource[n]
        qresult[bin] += (wresult[bin] * wval * (bval - mresult[bin])**2) / \
            (wresult[bin] + wval)
        wresult[bin] += wval
        bresult[bin] += wval*bval
        mresult[bin] += wval * (bval - mresult[bin]) / wresult[bin]
        used[bin] = 1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_profile2d(np.ndarray[np.int64_t, ndim=1] bins_x,
                  np.ndarray[np.int64_t, ndim=1] bins_y,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=1] bsource,
                  np.ndarray[np.float64_t, ndim=2] wresult,
                  np.ndarray[np.float64_t, ndim=2] bresult,
                  np.ndarray[np.float64_t, ndim=2] mresult,
                  np.ndarray[np.float64_t, ndim=2] qresult,
                  np.ndarray[np.float64_t, ndim=2] used):
    cdef int n, bini, binj
    cdef np.int64_t bin
    cdef np.float64_t wval, bval
    for n in range(bins_x.shape[0]):
        bini = bins_x[n]
        binj = bins_y[n]
        bval = bsource[n]
        wval = wsource[n]
        qresult[bini, binj] += (wresult[bini, binj] * wval * (bval - mresult[bini, binj])**2) / \
            (wresult[bini, binj] + wval)
        wresult[bini, binj] += wval
        bresult[bini, binj] += wval*bval
        mresult[bini, binj] += wval * (bval - mresult[bini, binj]) / wresult[bini, binj]
        used[bini, binj] = 1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_profile3d(np.ndarray[np.int64_t, ndim=1] bins_x,
                  np.ndarray[np.int64_t, ndim=1] bins_y,
                  np.ndarray[np.int64_t, ndim=1] bins_z,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=1] bsource,
                  np.ndarray[np.float64_t, ndim=3] wresult,
                  np.ndarray[np.float64_t, ndim=3] bresult,
                  np.ndarray[np.float64_t, ndim=3] mresult,
                  np.ndarray[np.float64_t, ndim=3] qresult,
                  np.ndarray[np.float64_t, ndim=3] used):
    cdef int n, bini, binj, bink
    cdef np.int64_t bin
    cdef np.float64_t wval, bval
    for n in range(bins_x.shape[0]):
        bini = bins_x[n]
        binj = bins_y[n]
        bink = bins_z[n]
        bval = bsource[n]
        wval = wsource[n]
        qresult[bini, binj, bink] += (wresult[bini, binj, bink] * wval * (bval - mresult[bini, binj, bink])**2) / \
            (wresult[bini, binj, bink] + wval)
        wresult[bini, binj, bink] += wval
        bresult[bini, binj, bink] += wval*bval
        mresult[bini, binj, bink] += wval * (bval - mresult[bini, binj, bink]) / wresult[bini, binj, bink]
        used[bini, binj, bink] = 1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def get_color_bounds(np.ndarray[np.float64_t, ndim=1] px,
                     np.ndarray[np.float64_t, ndim=1] py,
                     np.ndarray[np.float64_t, ndim=1] pdx,
                     np.ndarray[np.float64_t, ndim=1] pdy,
                     np.ndarray[np.float64_t, ndim=1] value,
                     np.float64_t leftx, np.float64_t rightx,
                     np.float64_t lefty, np.float64_t righty,
                     np.float64_t mindx = -1, np.float64_t maxdx = -1):
    cdef int i
    cdef np.float64_t mi = 1e100, ma = -1e100, v
    cdef int np = px.shape[0]
    with nogil:
        for i in range(np):
            v = value[i]
            if v < mi or v > ma:
                if px[i] + pdx[i] < leftx: continue
                if px[i] - pdx[i] > rightx: continue
                if py[i] + pdy[i] < lefty: continue
                if py[i] - pdy[i] > righty: continue
                if pdx[i] < mindx or pdy[i] < mindx: continue
                if maxdx > 0 and (pdx[i] > maxdx or pdy[i] > maxdx): continue
                if v < mi: mi = v
                if v > ma: ma = v
    return (mi, ma)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def kdtree_get_choices(np.ndarray[np.float64_t, ndim=3] data,
                       np.ndarray[np.float64_t, ndim=1] l_corner,
                       np.ndarray[np.float64_t, ndim=1] r_corner):
    cdef int i, j, k, dim, n_unique, best_dim, n_best, n_grids, addit, my_split
    n_grids = data.shape[0]
    cdef np.float64_t **uniquedims, *uniques, split
    uniquedims = <np.float64_t **> alloca(3 * sizeof(np.float64_t*))
    for i in range(3):
        uniquedims[i] = <np.float64_t *> \
                alloca(2*n_grids * sizeof(np.float64_t))
    my_max = 0
    best_dim = -1
    for dim in range(3):
        n_unique = 0
        uniques = uniquedims[dim]
        for i in range(n_grids):
            # Check for disqualification
            for j in range(2):
                #print "Checking against", i,j,dim,data[i,j,dim]
                if not (l_corner[dim] < data[i, j, dim] and
                        data[i, j, dim] < r_corner[dim]):
                    #print "Skipping ", data[i,j,dim]
                    continue
                skipit = 0
                # Add our left ...
                for k in range(n_unique):
                    if uniques[k] == data[i, j, dim]:
                        skipit = 1
                        #print "Identified", uniques[k], data[i,j,dim], n_unique
                        break
                if skipit == 0:
                    uniques[n_unique] = data[i, j, dim]
                    n_unique += 1
        if n_unique > my_max:
            best_dim = dim
            my_max = n_unique
            my_split = (n_unique-1)/2
    # I recognize how lame this is.
    cdef np.ndarray[np.float64_t, ndim=1] tarr = np.empty(my_max, dtype='float64')
    for i in range(my_max):
        #print "Setting tarr: ", i, uniquedims[best_dim][i]
        tarr[i] = uniquedims[best_dim][i]
    tarr.sort()
    split = tarr[my_split]
    cdef np.ndarray[np.uint8_t, ndim=1] less_ids = np.empty(n_grids, dtype='uint8')
    cdef np.ndarray[np.uint8_t, ndim=1] greater_ids = np.empty(n_grids, dtype='uint8')
    for i in range(n_grids):
        if data[i, 0, best_dim] < split:
            less_ids[i] = 1
        else:
            less_ids[i] = 0
        if data[i, 1, best_dim] > split:
            greater_ids[i] = 1
        else:
            greater_ids[i] = 0
    # Return out unique values
    return best_dim, split, less_ids.view("bool"), greater_ids.view("bool")

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def get_box_grids_level(np.ndarray[np.float64_t, ndim=1] left_edge,
                        np.ndarray[np.float64_t, ndim=1] right_edge,
                        int level,
                        np.ndarray[np.float64_t, ndim=2] left_edges,
                        np.ndarray[np.float64_t, ndim=2] right_edges,
                        np.ndarray[np.int32_t, ndim=2] levels,
                        np.ndarray[np.int32_t, ndim=1] mask,
                        int min_index = 0):
    cdef int i, n
    cdef int nx = left_edges.shape[0]
    cdef int inside 
    for i in range(nx):
        if i < min_index or levels[i,0] != level:
            mask[i] = 0
            continue
        inside = 1
        for n in range(3):
            if left_edge[n] >= right_edges[i,n] or \
               right_edge[n] <= left_edges[i,n]:
                inside = 0
                break
        if inside == 1: mask[i] = 1
        else: mask[i] = 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def get_box_grids_below_level(
                        np.ndarray[np.float64_t, ndim=1] left_edge,
                        np.ndarray[np.float64_t, ndim=1] right_edge,
                        int level,
                        np.ndarray[np.float64_t, ndim=2] left_edges,
                        np.ndarray[np.float64_t, ndim=2] right_edges,
                        np.ndarray[np.int32_t, ndim=2] levels,
                        np.ndarray[np.int32_t, ndim=1] mask,
                        int min_level = 0):
    cdef int i, n
    cdef int nx = left_edges.shape[0]
    cdef int inside 
    for i in range(nx):
        mask[i] = 0
        if levels[i,0] <= level and levels[i,0] >= min_level:
            inside = 1
            for n in range(3):
                if left_edge[n] >= right_edges[i,n] or \
                   right_edge[n] <= left_edges[i,n]:
                    inside = 0
                    break
            if inside == 1: mask[i] = 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_values_at_point(np.ndarray[np.float64_t, ndim=1] point,
                         np.ndarray[np.float64_t, ndim=2] left_edges,
                         np.ndarray[np.float64_t, ndim=2] right_edges,
                         np.ndarray[np.int32_t, ndim=2] dimensions,
                         field_names, grid_objects):
    # This iterates in order, first to last, and then returns with the first
    # one in which the point is located; this means if you order from highest
    # level to lowest, you will find the correct grid without consulting child
    # masking.  Note also that we will do a few relatively slow operations on
    # strings and whatnot, but they should not be terribly slow.
    cdef int ind[3], gi, fi
    cdef int nf = len(field_names)
    cdef np.float64_t dds
    cdef np.ndarray[np.float64_t, ndim=3] field
    cdef np.ndarray[np.float64_t, ndim=1] rv = np.zeros(nf, dtype='float64')
    for gi in range(left_edges.shape[0]):
        if not ((left_edges[gi,0] < point[0] < right_edges[gi,0])
            and (left_edges[gi,1] < point[1] < right_edges[gi,1])
            and (left_edges[gi,2] < point[2] < right_edges[gi,2])):
            continue
        # We found our grid!
        for fi in range(3):
            dds = ((right_edges[gi,fi] - left_edges[gi,fi])/
                   (<np.float64_t> dimensions[gi,fi]))
            ind[fi] = <int> ((point[fi] - left_edges[gi,fi])/dds)
        grid = grid_objects[gi]
        for fi in range(nf):
            field = grid[field_names[fi]]
            rv[fi] = field[ind[0], ind[1], ind[2]]
        return rv
    raise KeyError

@cython.cdivision(True)
def pixelize_cylinder(np.ndarray[np.float64_t, ndim=1] radius,
                      np.ndarray[np.float64_t, ndim=1] dradius,
                      np.ndarray[np.float64_t, ndim=1] theta,
                      np.ndarray[np.float64_t, ndim=1] dtheta,
                      int buff_size,
                      np.ndarray[np.float64_t, ndim=1] field,
                      np.float64_t rmax=-1.0) :

    cdef np.ndarray[np.float64_t, ndim=2] img
    cdef np.float64_t x, y, dx, dy, r0, theta0
    cdef np.float64_t r_i, theta_i, dr_i, dtheta_i, dthetamin
    cdef int i, pi, pj
    
    if rmax < 0.0 :
        imax = radius.argmax()
        rmax = radius[imax] + dradius[imax]
          
    img = np.zeros((buff_size, buff_size))
    extents = [-rmax, rmax] * 2
    dx = (extents[1] - extents[0]) / img.shape[0]
    dy = (extents[3] - extents[2]) / img.shape[1]
      
    dthetamin = dx / rmax
      
    for i in range(radius.shape[0]):

        r0 = radius[i]
        theta0 = theta[i]
        dr_i = dradius[i]
        dtheta_i = dtheta[i]

        theta_i = theta0 - dtheta_i
        while theta_i < theta0 + dtheta_i:
            r_i = r0 - dr_i
            while r_i < r0 + dr_i:
                if rmax <= r_i:
                    r_i += 0.5*dx 
                    continue
                x = r_i * math.cos(theta_i)
                y = r_i * math.sin(theta_i)
                pi = <int>((x + rmax)/dx)
                pj = <int>((y + rmax)/dy)
                img[pi, pj] = field[i]
                r_i += 0.5*dx 
            theta_i += dthetamin

    return img
