"""
Simple utilities that don't fit anywhere else



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import get_pbar
import numpy as np
from yt.units.yt_array import YTArray
cimport numpy as np
cimport cython
cimport libc.math as math
from libc.math cimport abs, sqrt
from yt.utilities.lib.fp_utils cimport fmin, fmax, i64min, i64max
from yt.geometry.selection_routines cimport _ensure_code

from libc.stdlib cimport malloc, free
from libc.string cimport strcmp

from cython.view cimport memoryview
from cython.view cimport array as cvarray
from cpython cimport buffer


cdef extern from "platform_dep.h":
    # NOTE that size_t might not be int
    void *alloca(int)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def new_bin_profile1d(np.ndarray[np.intp_t, ndim=1] bins_x,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=2] bsource,
                  np.ndarray[np.float64_t, ndim=1] wresult,
                  np.ndarray[np.float64_t, ndim=2] bresult,
                  np.ndarray[np.float64_t, ndim=2] mresult,
                  np.ndarray[np.float64_t, ndim=2] qresult,
                  np.ndarray[np.uint8_t, ndim=1, cast=True] used):
    cdef int n, fi, bin
    cdef np.float64_t wval, bval, oldwr
    cdef int nb = bins_x.shape[0]
    cdef int nf = bsource.shape[1]
    for n in range(nb):
        bin = bins_x[n]
        wval = wsource[n]
        oldwr = wresult[bin]
        wresult[bin] += wval
        for fi in range(nf):
            bval = bsource[n,fi]
            # qresult has to have the previous wresult
            qresult[bin,fi] += (oldwr * wval * (bval - mresult[bin,fi])**2) / \
                (oldwr + wval)
            bresult[bin,fi] += wval*bval
            # mresult needs the new wresult
            mresult[bin,fi] += wval * (bval - mresult[bin,fi]) / wresult[bin]
        used[bin] = 1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def new_bin_profile2d(np.ndarray[np.intp_t, ndim=1] bins_x,
                  np.ndarray[np.intp_t, ndim=1] bins_y,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=2] bsource,
                  np.ndarray[np.float64_t, ndim=2] wresult,
                  np.ndarray[np.float64_t, ndim=3] bresult,
                  np.ndarray[np.float64_t, ndim=3] mresult,
                  np.ndarray[np.float64_t, ndim=3] qresult,
                  np.ndarray[np.uint8_t, ndim=2, cast=True] used):
    cdef int n, fi, bin_x, bin_y
    cdef np.float64_t wval, bval, oldwr
    cdef int nb = bins_x.shape[0]
    cdef int nf = bsource.shape[1]
    for n in range(nb):
        bin_x = bins_x[n]
        bin_y = bins_y[n]
        wval = wsource[n]
        oldwr = wresult[bin_x, bin_y]
        wresult[bin_x,bin_y] += wval
        for fi in range(nf):
            bval = bsource[n,fi]
            # qresult has to have the previous wresult
            qresult[bin_x,bin_y,fi] += (oldwr * wval * (bval - mresult[bin_x,bin_y,fi])**2) / \
                (oldwr + wval)
            bresult[bin_x,bin_y,fi] += wval*bval
            # mresult needs the new wresult
            mresult[bin_x,bin_y,fi] += wval * (bval - mresult[bin_x,bin_y,fi]) / wresult[bin_x,bin_y]
        used[bin_x,bin_y] = 1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def new_bin_profile3d(np.ndarray[np.intp_t, ndim=1] bins_x,
                  np.ndarray[np.intp_t, ndim=1] bins_y,
                  np.ndarray[np.intp_t, ndim=1] bins_z,
                  np.ndarray[np.float64_t, ndim=1] wsource,
                  np.ndarray[np.float64_t, ndim=2] bsource,
                  np.ndarray[np.float64_t, ndim=3] wresult,
                  np.ndarray[np.float64_t, ndim=4] bresult,
                  np.ndarray[np.float64_t, ndim=4] mresult,
                  np.ndarray[np.float64_t, ndim=4] qresult,
                  np.ndarray[np.uint8_t, ndim=3, cast=True] used):
    cdef int n, fi, bin_x, bin_y, bin_z
    cdef np.float64_t wval, bval, oldwr
    cdef int nb = bins_x.shape[0]
    cdef int nf = bsource.shape[1]
    for n in range(nb):
        bin_x = bins_x[n]
        bin_y = bins_y[n]
        bin_z = bins_z[n]
        wval = wsource[n]
        oldwr = wresult[bin_x, bin_y, bin_z]
        wresult[bin_x,bin_y,bin_z] += wval
        for fi in range(nf):
            bval = bsource[n,fi]
            # qresult has to have the previous wresult
            qresult[bin_x,bin_y,bin_z,fi] += \
                (oldwr * wval * (bval - mresult[bin_x,bin_y,bin_z,fi])**2) / \
                (oldwr + wval)
            bresult[bin_x,bin_y,bin_z,fi] += wval*bval
            # mresult needs the new wresult
            mresult[bin_x,bin_y,bin_z,fi] += wval * \
                (bval - mresult[bin_x,bin_y,bin_z,fi]) / \
                 wresult[bin_x,bin_y,bin_z]
        used[bin_x,bin_y,bin_z] = 1
    return

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
def lines(np.float64_t[:,:,:] image,
          np.int64_t[:] xs,
          np.int64_t[:] ys,
          np.float64_t[:,:] colors,
          int points_per_color=1,
          int thick=1,
	      int flip=0,
          int crop = 0):

    cdef int nx = image.shape[0]
    cdef int ny = image.shape[1]
    cdef int nl = xs.shape[0]
    cdef np.float64_t alpha[4]
    cdef np.float64_t outa
    cdef int i, j, xi, yi
    cdef int dx, dy, sx, sy, e2, err
    cdef np.int64_t x0, x1, y0, y1
    cdef int has_alpha = (image.shape[2] == 4)
    cdef int no_color = (image.shape[2] < 3)
    for j in range(0, nl, 2):
        # From wikipedia http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
        x0 = xs[j]
        y0 = ys[j]
        x1 = xs[j+1]
        y1 = ys[j+1]
        dx = abs(x1-x0)
        dy = abs(y1-y0)
        if crop == 1 and (dx > nx/2.0 or dy > ny/2.0):
            continue
        err = dx - dy

        if no_color:
            for i in range(4):
                alpha[i] = colors[j, 0]
        elif has_alpha:
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
            if (x0 < thick and sx == -1): break
            elif (x0 >= nx-thick+1 and sx == 1): break
            elif (y0 < thick and sy == -1): break
            elif (y0 >= ny-thick+1 and sy == 1): break
            if x0 >= thick and x0 < nx-thick and y0 >= thick and y0 < ny-thick:
                for xi in range(x0-thick/2, x0+(1+thick)/2):
                    for yi in range(y0-thick/2, y0+(1+thick)/2):
                        if flip:
                            yi0 = ny - yi
                        else:
                            yi0 = yi

                        if no_color:
                            image[xi, yi0, 0] = fmin(alpha[0], image[xi, yi0, 0])
                        elif has_alpha:
                            image[xi, yi0, 3] = outa = alpha[3] + image[xi, yi0, 3]*(1-alpha[3])
                            if outa != 0.0:
                                outa = 1.0/outa
                            for i in range(3):
                                image[xi, yi0, i] = \
                                        ((1.-alpha[3])*image[xi, yi0, i]*image[xi, yi0, 3]
                                         + alpha[3]*alpha[i])*outa
                        else:
                            for i in range(3):
                                image[xi, yi0, i] = \
                                        (1.-alpha[i])*image[xi,yi0,i] + alpha[i]


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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def zlines(np.ndarray[np.float64_t, ndim=3] image,
        np.ndarray[np.float64_t, ndim=2] zbuffer,
        np.ndarray[np.int64_t, ndim=1] xs,
        np.ndarray[np.int64_t, ndim=1] ys,
        np.ndarray[np.float64_t, ndim=1] zs,
        np.ndarray[np.float64_t, ndim=2] colors,
        int points_per_color=1,
        int thick=1,
        int flip=0,
        int crop = 0):

    cdef int nx = image.shape[0]
    cdef int ny = image.shape[1]
    cdef int nl = xs.shape[0]
    cdef np.float64_t[:] alpha
    cdef int i, j, c
    cdef int dx, dy, sx, sy, e2, err
    cdef np.int64_t x0, x1, y0, y1, yi0
    cdef np.float64_t z0, z1, dzx, dzy
    alpha = np.zeros(4)
    for j in range(0, nl, 2):
        # From wikipedia http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
        x0 = xs[j]
        y0 = ys[j]
        x1 = xs[j+1]
        y1 = ys[j+1]
        z0 = zs[j]
        z1 = zs[j+1]
        dx = abs(x1-x0)
        dy = abs(y1-y0)
        dzx = (z1-z0) / (dx**2 + dy**2) * dx
        dzy = (z1-z0) / (dx**2 + dy**2) * dy
        err = dx - dy
        if crop == 1 and (dx > nx/2.0 or dy > ny/2.0):
            continue

        c = j/points_per_color/2

        for i in range(3):
            alpha[i] = colors[c, i] * colors[c, 3]
        alpha[3] = colors[c, 3]

        if x0 < x1:
            sx = 1
        else:
            sx = -1
        if y0 < y1:
            sy = 1
        else:
            sy = -1
        while(1):
            if (x0 < thick and sx == -1): break
            elif (x0 >= nx-thick+1 and sx == 1): break
            elif (y0 < thick and sy == -1): break
            elif (y0 >= ny-thick+1 and sy == 1): break
            if x0 >= thick and x0 < nx-thick and y0 >= thick and y0 < ny-thick:
                for _ in range(x0-thick/2, x0+(1+thick)/2):
                    for yi in range(y0-thick/2, y0+(1+thick)/2):
                        if flip:
                            yi0 = ny - yi
                        else:
                            yi0 = yi
                        if z0 < zbuffer[x0, yi0]:
                            if alpha[3] != 1.0:
                                talpha = image[x0, yi0, 3]
                                image[x0, yi0, 3] = alpha[3] + talpha * (1 - alpha[3])
                                for i in range(3):
                                    if image[x0, yi0, 3] == 0.0:
                                        image[x0, yi0, i] = 0.0
                                    else:
                                        image[x0, yi0, i] = (alpha[3]*alpha[i] + image[x0, yi0, i]*talpha*(1.0-alpha[3]))/image[x0,yi0,3]
                            else:
                                for i in range(4):
                                    image[x0, yi0, i] = alpha[i]
                            if (1.0 - image[x0, yi0, 3] < 1.0e-4):
                                image[x0, yi0, 3] = 1.0
                                zbuffer[x0, yi0] = z0

            if (x0 == x1 and y0 == y1):
                break
            e2 = 2*err
            if e2 > -dy:
                err = err - dy
                x0 += sx
                z0 += dzx
            if e2 < dx :
                err = err + dx
                y0 += sy
                z0 += dzy
        # assert(np.abs(z0 - z1) < 1.0e-3 * (np.abs(z0) + np.abs(z1)))
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def zpoints(np.ndarray[np.float64_t, ndim=3] image,
        np.ndarray[np.float64_t, ndim=2] zbuffer,
        np.ndarray[np.int64_t, ndim=1] xs,
        np.ndarray[np.int64_t, ndim=1] ys,
        np.ndarray[np.float64_t, ndim=1] zs,
        np.ndarray[np.float64_t, ndim=2] colors,
        int points_per_color=1,
        int thick=1,
        int flip=0):

    cdef int nx = image.shape[0]
    cdef int ny = image.shape[1]
    cdef int nl = xs.shape[0]
    cdef np.float64_t[:] alpha
    cdef np.float64_t talpha
    cdef int i, j, c
    cdef np.int64_t x0, y0, yi0
    cdef np.float64_t z0
    alpha = np.zeros(4)
    for j in range(0, nl):
        x0 = xs[j]
        y0 = ys[j]
        z0 = zs[j]
        if (x0 < 0 or x0 >= nx): continue
        if (y0 < 0 or y0 >= ny): continue
        c = j/points_per_color
        for i in range(3):
            alpha[i] = colors[c, i] * colors[c, 3]
        alpha[3] = colors[c, 3]
        if flip:
            yi0 = ny - y0
        else:
            yi0 = y0

        if z0 < zbuffer[x0, yi0]:
            if alpha[3] != 1.0:
                talpha = image[x0, yi0, 3]
                image[x0, yi0, 3] = alpha[3] + talpha * (1 - alpha[3])
                for i in range(3):
                    image[x0, yi0, i] = (alpha[3]*alpha[i] + image[x0, yi0, i]*talpha*(1.0-alpha[3]))/image[x0,yi0,3]
                    if image[x0, yi0, 3] == 0.0:
                        image[x0, yi0, i] = 0.0
            else:
                for i in range(4):
                    image[x0, yi0, i] = alpha[i]
            if (1.0 - image[x0, yi0, 3] < 1.0e-4):
                zbuffer[x0, yi0] = z0
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
def kdtree_get_choices(np.ndarray[np.float64_t, ndim=3] data,
                       np.ndarray[np.float64_t, ndim=1] l_corner,
                       np.ndarray[np.float64_t, ndim=1] r_corner):
    cdef int i, j, k, dim, n_unique, best_dim, n_grids, my_split
    n_grids = data.shape[0]
    cdef np.float64_t **uniquedims
    cdef np.float64_t *uniques
    cdef np.float64_t split
    uniquedims = <np.float64_t **> alloca(3 * sizeof(np.float64_t*))
    for i in range(3):
        uniquedims[i] = <np.float64_t *> \
                alloca(2*n_grids * sizeof(np.float64_t))
    my_max = 0
    best_dim = -1
    my_split = -1
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
    if my_split < 0:
        raise RuntimeError
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
    cdef np.float64_t eps = np.finfo(np.float64).eps
    for i in range(nx):
        if i < min_index or levels[i,0] != level:
            mask[i] = 0
            continue
        inside = 1
        for n in range(3):
            if (right_edges[i,n] - left_edge[n]) <= eps or \
               (right_edge[n] - left_edges[i,n]) <= eps:
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
    cdef np.float64_t eps = np.finfo(np.float64).eps
    for i in range(nx):
        mask[i] = 0
        if levels[i,0] <= level and levels[i,0] >= min_level:
            inside = 1
            for n in range(3):
                if (right_edges[i,n] - left_edge[n]) <= eps or \
                   (right_edge[n] - left_edges[i,n]) <= eps:
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
    cdef int ind[3]
    cdef int gi, fi, nf = len(field_names)
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

#@cython.cdivision(True)
#@cython.boundscheck(False)
#@cython.wraparound(False)
def obtain_rvec(data):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] xf
    cdef np.ndarray[np.float64_t, ndim=1] yf
    cdef np.ndarray[np.float64_t, ndim=1] zf
    cdef np.ndarray[np.float64_t, ndim=2] rf
    cdef np.ndarray[np.float64_t, ndim=3] xg
    cdef np.ndarray[np.float64_t, ndim=3] yg
    cdef np.ndarray[np.float64_t, ndim=3] zg
    cdef np.ndarray[np.float64_t, ndim=4] rg
    cdef np.float64_t c[3]
    cdef int i, j, k
    center = data.get_field_parameter("center")
    c[0] = center[0]; c[1] = center[1]; c[2] = center[2]
    if len(data['x'].shape) == 1:
        # One dimensional data
        xf = data['x']
        yf = data['y']
        zf = data['z']
        rf = YTArray(np.empty((3, xf.shape[0]), 'float64'), xf.units)
        for i in range(xf.shape[0]):
            rf[0, i] = xf[i] - c[0]
            rf[1, i] = yf[i] - c[1]
            rf[2, i] = zf[i] - c[2]
        return rf
    else:
        # Three dimensional data
        xg = data['x']
        yg = data['y']
        zg = data['z']
        shape = (3, xg.shape[0], xg.shape[1], xg.shape[2])
        rg = YTArray(np.empty(shape, 'float64'), xg.units)
        #rg = YTArray(rg, xg.units)
        for i in range(xg.shape[0]):
            for j in range(xg.shape[1]):
                for k in range(xg.shape[2]):
                    rg[0,i,j,k] = xg[i,j,k] - c[0]
                    rg[1,i,j,k] = yg[i,j,k] - c[1]
                    rg[2,i,j,k] = zg[i,j,k] - c[2]
        return rg

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def obtain_rv_vec(data, field_names = ("velocity_x",
                                       "velocity_y",
                                       "velocity_z"),
                  bulk_vector = "bulk_velocity"):
    # This is just to let the pointers exist and whatnot.  We can't cdef them
    # inside conditionals.
    cdef np.ndarray[np.float64_t, ndim=1] vxf
    cdef np.ndarray[np.float64_t, ndim=1] vyf
    cdef np.ndarray[np.float64_t, ndim=1] vzf
    cdef np.ndarray[np.float64_t, ndim=2] rvf
    cdef np.ndarray[np.float64_t, ndim=3] vxg
    cdef np.ndarray[np.float64_t, ndim=3] vyg
    cdef np.ndarray[np.float64_t, ndim=3] vzg
    cdef np.ndarray[np.float64_t, ndim=4] rvg
    cdef np.float64_t bv[3]
    cdef int i, j, k
    bulk_vector = data.get_field_parameter(bulk_vector)
    if len(data[field_names[0]].shape) == 1:
        # One dimensional data
        vxf = data[field_names[0]].astype("float64")
        vyf = data[field_names[1]].astype("float64")
        vzf = data[field_names[2]].astype("float64")
        vyf.convert_to_units(vxf.units)
        vzf.convert_to_units(vxf.units)
        rvf = YTArray(np.empty((3, vxf.shape[0]), 'float64'), vxf.units)
        if bulk_vector is None:
            bv[0] = bv[1] = bv[2] = 0.0
        else:
            bulk_vector = bulk_vector.in_units(vxf.units)
            bv[0] = bulk_vector[0]
            bv[1] = bulk_vector[1]
            bv[2] = bulk_vector[2]
        for i in range(vxf.shape[0]):
            rvf[0, i] = vxf[i] - bv[0]
            rvf[1, i] = vyf[i] - bv[1]
            rvf[2, i] = vzf[i] - bv[2]
        return rvf
    else:
        # Three dimensional data
        vxg = data[field_names[0]].astype("float64")
        vyg = data[field_names[1]].astype("float64")
        vzg = data[field_names[2]].astype("float64")
        vyg.convert_to_units(vxg.units)
        vzg.convert_to_units(vxg.units)
        shape = (3, vxg.shape[0], vxg.shape[1], vxg.shape[2])
        rvg = YTArray(np.empty(shape, 'float64'), vxg.units)
        if bulk_vector is None:
            bv[0] = bv[1] = bv[2] = 0.0
        else:
            bulk_vector = bulk_vector.in_units(vxg.units)
            bv[0] = bulk_vector[0]
            bv[1] = bulk_vector[1]
            bv[2] = bulk_vector[2]
        for i in range(vxg.shape[0]):
            for j in range(vxg.shape[1]):
                for k in range(vxg.shape[2]):
                    rvg[0,i,j,k] = vxg[i,j,k] - bv[0]
                    rvg[1,i,j,k] = vyg[i,j,k] - bv[1]
                    rvg[2,i,j,k] = vzg[i,j,k] - bv[2]
        return rvg

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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def fill_region(input_fields, output_fields,
                np.int32_t output_level,
                np.ndarray[np.int64_t, ndim=1] left_index,
                np.ndarray[np.int64_t, ndim=2] ipos,
                np.ndarray[np.int64_t, ndim=1] ires,
                np.ndarray[np.int64_t, ndim=1] level_dims,
                np.int64_t refine_by = 2
                ):
    cdef int i, n
    cdef np.int64_t tot = 0, oi, oj, ok, rf
    cdef np.int64_t iind[3]
    cdef np.int64_t oind[3]
    cdef np.int64_t dim[3]
    cdef np.ndarray[np.float64_t, ndim=3] ofield
    cdef np.ndarray[np.float64_t, ndim=1] ifield
    nf = len(input_fields)
    # The variable offsets governs for each dimension and each possible
    # wrapping if we do it.  Then the wi, wj, wk indices check into each
    # [dim][wrap] inside the loops.
    cdef int wi, wj, wk
    cdef int offsets[3][3]
    cdef np.int64_t off
    for i in range(3):
        # Offsets here is a way of accounting for periodicity.  It keeps track
        # of how to offset our grid as we loop over the icoords.
        dim[i] = output_fields[0].shape[i]
        offsets[i][0] = offsets[i][2] = 0
        offsets[i][1] = 1
        if left_index[i] < 0:
            offsets[i][2] = 1
        if left_index[i] + dim[i] >= level_dims[i]:
            offsets[i][0] = 1
    for n in range(nf):
        tot = 0
        ofield = output_fields[n]
        ifield = input_fields[n]
        for i in range(ipos.shape[0]):
            rf = refine_by**(output_level - ires[i])
            for wi in range(3):
                if offsets[0][wi] == 0: continue
                off = (left_index[0] + level_dims[0]*(wi-1))
                iind[0] = ipos[i, 0] * rf - off
                # rf here is the "refinement factor", or, the number of zones
                # that this zone could potentially contribute to our filled
                # grid.
                for oi in range(rf):
                    # Now we need to apply our offset
                    oind[0] = oi + iind[0]
                    if oind[0] < 0:
                        continue
                    elif oind[0] >= dim[0]:
                        break
                    for wj in range(3):
                        if offsets[1][wj] == 0: continue
                        off = (left_index[1] + level_dims[1]*(wj-1))
                        iind[1] = ipos[i, 1] * rf - off
                        for oj in range(rf):
                            oind[1] = oj + iind[1]
                            if oind[1] < 0:
                                continue
                            elif oind[1] >= dim[1]:
                                break
                            for wk in range(3):
                                if offsets[2][wk] == 0: continue
                                off = (left_index[2] + level_dims[2]*(wk-1))
                                iind[2] = ipos[i, 2] * rf - off
                                for ok in range(rf):
                                    oind[2] = ok + iind[2]
                                    if oind[2] < 0:
                                        continue
                                    elif oind[2] >= dim[2]:
                                        break
                                    ofield[oind[0], oind[1], oind[2]] = \
                                        ifield[i]
                                    tot += 1
    return tot

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def fill_region_float(np.ndarray[np.float64_t, ndim=2] fcoords,
                      np.ndarray[np.float64_t, ndim=2] fwidth,
                      np.ndarray[np.float64_t, ndim=1] data,
                      np.ndarray[np.float64_t, ndim=1] box_left_edge,
                      np.ndarray[np.float64_t, ndim=1] box_right_edge,
                      np.ndarray[np.float64_t, ndim=3] dest,
                      int antialias = 1,
                      period = None,
                      int check_period = 1):
    cdef np.float64_t ds_period[3]
    cdef np.float64_t box_dds[3]
    cdef np.float64_t box_idds[3]
    cdef np.float64_t width[3]
    cdef np.float64_t LE[3]
    cdef np.float64_t RE[3]
    cdef np.int64_t i, j, k, p, xi, yi
    cdef np.int64_t dims[3]
    cdef np.int64_t ld[3]
    cdef np.int64_t ud[3]
    cdef np.float64_t overlap[3]
    cdef np.float64_t dsp
    cdef np.float64_t osp[3]
    cdef np.float64_t odsp[3]
    cdef np.float64_t sp[3]
    cdef np.float64_t lfd[3]
    cdef np.float64_t ufd[3]
    # These are the temp vars we get from the arrays
    # Some periodicity helpers
    cdef int diter[3][2]
    cdef np.float64_t diterv[3][2]
    if period is not None:
        for i in range(3):
            ds_period[i] = period[i]
    else:
        ds_period[0] = ds_period[1] = ds_period[2] = 0.0
    box_left_edge = _ensure_code(box_left_edge)
    box_right_edge = _ensure_code(box_right_edge)
    _ensure_code(fcoords)
    _ensure_code(fwidth)
    for i in range(3):
        LE[i] = box_left_edge[i]
        RE[i] = box_right_edge[i]
        width[i] = RE[i] - LE[i]
        dims[i] = dest.shape[i]
        box_dds[i] = width[i] / dims[i]
        box_idds[i] = 1.0/box_dds[i]
        diter[i][0] = diter[i][1] = 0
        diterv[i][0] = diterv[i][1] = 0.0
        overlap[i] = 1.0 
    with nogil:
        for p in range(fcoords.shape[0]):
            for i in range(3):
               diter[i][1] = 999
               odsp[i] = fwidth[p,i]*0.5
               osp[i] = fcoords[p,i] # already centered
               overlap[i] = 1.0
            dsp = data[p]
            if check_period == 1:
                for i in range(3):
                    if (osp[i] - odsp[i] < LE[i]):
                        diter[i][1] = +1
                        diterv[i][1] = ds_period[i]
                    elif (osp[i] + odsp[i] > RE[i]):
                        diter[i][1] = -1
                        diterv[i][1] = -ds_period[i]
            for xi in range(2):
                if diter[0][xi] == 999: continue
                sp[0] = osp[0] + diterv[0][xi]
                if (sp[0] + odsp[0] < LE[0]) or (sp[0] - odsp[0] > RE[0]): continue
                for yi in range(2):
                    if diter[1][yi] == 999: continue
                    sp[1] = osp[1] + diterv[1][yi]
                    if (sp[1] + odsp[1] < LE[1]) or (sp[1] - odsp[1] > RE[1]): continue
                    for zi in range(2):
                        if diter[2][zi] == 999: continue
                        sp[2] = osp[2] + diterv[2][yi]
                        if (sp[2] + odsp[2] < LE[2]) or (sp[2] - odsp[2] > RE[2]): continue
                        for i in range(3):
                            ld[i] = <np.int64_t> fmax(((sp[i]-odsp[i]-LE[i])*box_idds[i]),0)
                            # NOTE: This is a different way of doing it than in the C
                            # routines.  In C, we were implicitly casting the
                            # initialization to int, but *not* the conditional, which
                            # was allowed an extra value:
                            #     for(j=lc;j<rc;j++)
                            # here, when assigning lc (double) to j (int) it got
                            # truncated, but no similar truncation was done in the
                            # comparison of j to rc (double).  So give ourselves a
                            # bonus row and bonus column here.
                            ud[i] = <np.int64_t> fmin(((sp[i]+odsp[i]-LE[i])*box_idds[i] + 1), dims[i])
                        for i in range(ld[0], ud[0]):
                            if antialias == 1:
                                lfd[0] = box_dds[0] * i + LE[0]
                                ufd[0] = box_dds[0] * (i + 1) + LE[0]
                                overlap[0] = ((fmin(ufd[0], sp[0]+odsp[0])
                                           - fmax(lfd[0], (sp[0]-odsp[0])))*box_idds[0])
                            if overlap[0] < 0.0: continue
                            for j in range(ld[1], ud[1]):
                                if antialias == 1:
                                    lfd[1] = box_dds[1] * j + LE[1]
                                    ufd[1] = box_dds[1] * (j + 1) + LE[1]
                                    overlap[1] = ((fmin(ufd[1], sp[1]+odsp[1])
                                               - fmax(lfd[1], (sp[1]-odsp[1])))*box_idds[1])
                                if overlap[1] < 0.0: continue
                                for k in range(ld[2], ud[2]):
                                    if antialias == 1:
                                        lfd[2] = box_dds[2] * k + LE[2]
                                        ufd[2] = box_dds[2] * (k + 1) + LE[2]
                                        overlap[2] = ((fmin(ufd[2], sp[2]+odsp[2])
                                                   - fmax(lfd[2], (sp[2]-odsp[2])))*box_idds[2])
                                        if overlap[2] < 0.0: continue
                                        dest[i,j,k] += dsp * (overlap[0]*overlap[1]*overlap[2])
                                    else:
                                        dest[i,j,k] = dsp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def gravitational_binding_energy(
        np.float64_t[:] mass,
        np.float64_t[:] x,
        np.float64_t[:] y,
        np.float64_t[:] z,
        int truncate,
        np.float64_t kinetic):

    cdef int q_outer, q_inner, n_q, i
    cdef np.float64_t mass_o, x_o, y_o, z_o
    cdef np.float64_t mass_i, x_i, y_i, z_i
    cdef np.float64_t this_potential, total_potential
    total_potential = 0.

    i = 0
    n_q = mass.size
    pbar = get_pbar("Calculating potential for %d cells" % n_q,
                    0.5 * (n_q**2 - n_q))
    for q_outer in range(n_q - 1):
        this_potential = 0.
        mass_o = mass[q_outer]
        x_o = x[q_outer]
        y_o = y[q_outer]
        z_o = z[q_outer]
        for q_inner in range(q_outer + 1, n_q):
            mass_i = mass[q_inner]
            x_i = x[q_inner]
            y_i = y[q_inner]
            z_i = z[q_inner]
            this_potential += mass_o * mass_i / \
              sqrt((x_i - x_o) * (x_i - x_o) +
                   (y_i - y_o) * (y_i - y_o) +
                   (z_i - z_o) * (z_i - z_o))
        i += n_q - q_outer
        pbar.update(i)
        total_potential += this_potential
        if truncate and total_potential / kinetic > 1.:
            break
    pbar.finish()

    return total_potential

# The OnceIndirect code is from:
# http://stackoverflow.com/questions/10465091/assembling-a-cython-memoryview-from-numpy-arrays/12991519#12991519
# This is under the CC-BY-SA license.

cdef class OnceIndirect:
    cdef object _objects
    cdef void** buf
    cdef int ndim
    cdef int n_rows
    cdef int buf_len
    cdef Py_ssize_t* shape
    cdef Py_ssize_t* strides
    cdef Py_ssize_t* suboffsets
    cdef Py_ssize_t itemsize
    cdef bytes format
    cdef int is_readonly

    def __cinit__(self, object rows, want_writable=True, want_format=True, allow_indirect=False):
        """
        Set want_writable to False if you don't want writable data. (This may
        prevent copies.)
        Set want_format to False if your input doesn't support PyBUF_FORMAT (unlikely)
        Set allow_indirect to True if you are ok with the memoryview being indirect
        in dimensions other than the first. (This may prevent copies.)

        An example usage:

        cdef double[::cython.view.indirect, ::1] vectors =
            OnceIndirect([object.vector for object in objects])
        """
        demand = buffer.PyBUF_INDIRECT if allow_indirect else buffer.PyBUF_STRIDES
        if want_writable:
            demand |= buffer.PyBUF_WRITABLE
        if want_format:
            demand |= buffer.PyBUF_FORMAT
        self._objects = [memoryview(row, demand) for row in rows]
        self.n_rows = len(self._objects)
        self.buf_len = sizeof(void*) * self.n_rows
        self.buf = <void**>malloc(self.buf_len)
        self.ndim = 1 + self._objects[0].ndim
        self.shape = <Py_ssize_t*>malloc(sizeof(Py_ssize_t) * self.ndim)
        self.strides = <Py_ssize_t*>malloc(sizeof(Py_ssize_t) * self.ndim)
        self.suboffsets = <Py_ssize_t*>malloc(sizeof(Py_ssize_t) * self.ndim)

        cdef memoryview example_obj = self._objects[0]
        self.itemsize = example_obj.itemsize

        if want_format:
            self.format = example_obj.view.format
        else:
            self.format = None
        self.is_readonly |= example_obj.view.readonly

        for dim in range(self.ndim):
            if dim == 0:
                self.shape[dim] = self.n_rows
                self.strides[dim] = sizeof(void*)
                self.suboffsets[dim] = 0
            else:
                self.shape[dim] = example_obj.view.shape[dim - 1]
                self.strides[dim] = example_obj.view.strides[dim - 1]
                if example_obj.view.suboffsets == NULL:
                    self.suboffsets[dim] = -1
                else:
                    self.suboffsets[dim] = example_obj.suboffsets[dim - 1]

        cdef memoryview obj
        cdef int i = 0
        for obj in self._objects:
            assert_similar(example_obj, obj)
            self.buf[i] = obj.view.buf
            i += 1

    def __getbuffer__(self, Py_buffer* buff, int flags):
        if (flags & buffer.PyBUF_INDIRECT) != buffer.PyBUF_INDIRECT:
            raise Exception("don't want to copy data")
        if flags & buffer.PyBUF_WRITABLE and self.is_readonly:
            raise Exception("couldn't provide writable, you should have demanded it earlier")
        if flags & buffer.PyBUF_FORMAT:
            if self.format is None:
                raise Exception("couldn't provide format, you should have demanded it earlier")
            buff.format = self.format
        else:
            buff.format = NULL

        buff.buf = <void*>self.buf
        buff.obj = self
        buff.len = self.buf_len
        buff.readonly = self.is_readonly
        buff.ndim = self.ndim
        buff.shape = self.shape
        buff.strides = self.strides
        buff.suboffsets = self.suboffsets
        buff.itemsize = self.itemsize
        buff.internal = NULL

    def __dealloc__(self):
        free(self.buf)
        free(self.shape)
        free(self.strides)
        free(self.suboffsets)

cdef int assert_similar(memoryview left_, memoryview right_) except -1:
    cdef Py_buffer left = left_.view
    cdef Py_buffer right = right_.view
    assert left.ndim == right.ndim
    cdef int i
    for i in range(left.ndim):
        assert left.shape[i] == right.shape[i], (left_.shape, right_.shape)
        assert left.strides[i] == right.strides[i], (left_.strides, right_.strides)

    if left.suboffsets == NULL:
        assert right.suboffsets == NULL, (left_.suboffsets, right_.suboffsets)
    else:
        for i in range(left.ndim):
            assert left.suboffsets[i] == right.suboffsets[i], (left_.suboffsets, right_.suboffsets)

    if left.format == NULL:
        assert right.format == NULL, (bytes(left.format), bytes(right.format))
    else:
        #alternatively, compare as Python strings:
        #assert bytes(left.format) == bytes(right.format)
        assert strcmp(left.format, right.format) == 0, (bytes(left.format), bytes(right.format))
    return 0
