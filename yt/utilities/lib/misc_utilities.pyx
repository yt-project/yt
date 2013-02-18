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
def lines(np.ndarray[np.float64_t, ndim=3] image, 
          np.ndarray[np.int64_t, ndim=1] xs,
          np.ndarray[np.int64_t, ndim=1] ys,
          np.ndarray[np.float64_t, ndim=2] colors,
          int points_per_color=1):
    
    cdef int nx = image.shape[0]
    cdef int ny = image.shape[1]
    cdef int nl = xs.shape[0]
    cdef np.float64_t alpha[4]
    cdef int i, j
    cdef int dx, dy, sx, sy, e2, err
    cdef np.int64_t x0, x1, y0, y1
    cdef int has_alpha = (image.shape[-1] == 4)
    for j in range(0, nl, 2):
        # From wikipedia http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
        x0 = xs[j]; y0 = ys[j]; x1 = xs[j+1]; y1 = ys[j+1]
        dx = abs(x1-x0)
        dy = abs(y1-y0)
        err = dx - dy
        if has_alpha:
            for i in range(4):
                alpha[i] = colors[j/points_per_color,i]
        else:
            for i in range(3):
                alpha[i] = colors[j/points_per_color,3]*\
                        colors[j/points_per_color,i]
        if x0 < x1: 
            sx = 1
        else:
            sx = -1
        if y0 < y1:
            sy = 1
        else:
            sy = -1
        while(1): 
            if (x0 < 0 and sx == -1): break
            elif (x0 >= nx and sx == 1): break
            elif (y0 < 0 and sy == -1): break
            elif (y0 >= nx and sy == 1): break
            if (x0 >=0 and x0 < nx and y0 >= 0 and y0 < ny):
                if has_alpha:
                    for i in range(4):
                        image[x0,y0,i] = (1.-alpha[i])*image[x0,y0,i] + alpha[i]
                else:
                    for i in range(3):
                        image[x0,y0,i] = (1.-alpha[i])*image[x0,y0,i] + alpha[i]

            if (x0 == x1 and y0 == y1):
                break
            e2 = 2*err
            if e2 > -dy:
                err = err - dy
                x0 += sx
            if e2 < dx :
                err = err + dx
                y0 += sy
    return 

def rotate_vectors(np.ndarray[np.float64_t, ndim=3] vecs,
        np.ndarray[np.float64_t, ndim=2] R):
    cdef int nx = vecs.shape[0]
    cdef int ny = vecs.shape[1]
    rotated = np.empty((nx,ny,3),dtype='float64') 
    for i in range(nx):
        for j in range(ny):
            for k in range(3):
                rotated[i,j,k] =\
                    R[k,0]*vecs[i,j,0]+R[k,1]*vecs[i,j,1]+R[k,2]*vecs[i,j,2]
    return rotated

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


def grow_flagging_field(oofield):
    cdef np.ndarray[np.uint8_t, ndim=3] ofield = oofield.astype("uint8")
    cdef np.ndarray[np.uint8_t, ndim=3] nfield
    nfield = np.zeros_like(ofield)
    cdef int i, j, k, ni, nj, nk
    cdef int oi, oj, ok
    for ni in range(ofield.shape[0]):
        for nj in range(ofield.shape[1]):
            for nk in range(ofield.shape[2]):
                for oi in range(3):
                    i = ni + (oi - 1)
                    if i < 0 or i >= ofield.shape[0]: continue
                    for oj in range(3):
                        j = nj + (oj - 1)
                        if j < 0 or j >= ofield.shape[1]: continue
                        for ok in range(3):
                            k = nk + (ok - 1)
                            if k < 0 or k >= ofield.shape[2]: continue
                            if ofield[i, j, k] == 1:
                                nfield[ni, nj, nk] = 1
    return nfield.astype("bool")
