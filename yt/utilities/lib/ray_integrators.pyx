"""
Simle integrators for the radiative transfer equation



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
from libc.stdlib cimport malloc, free, abs

cdef extern from "math.h":
    double exp(double x)
    float expf(float x)

@cython.boundscheck(False)
def Transfer3D(np.ndarray[np.float64_t, ndim=3] i_s,
               np.ndarray[np.float64_t, ndim=4] o_s,
               np.ndarray[np.float64_t, ndim=4] e,
               np.ndarray[np.float64_t, ndim=4] a,
               int imin, int imax, int jmin, int jmax,
               int kmin, int kmax, int istride, int jstride,
               np.float64_t dx):
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
    cdef int k
    cdef int n, nn
    nn = o_s.shape[3] # This might be slow
    cdef np.float64_t *temp = <np.float64_t *>malloc(sizeof(np.float64_t) * nn)
    for i in range((imax-imin)*istride):
        ii = i + imin*istride
        for j in range((jmax-jmin)*jstride):
            jj = j + jmin*jstride
            # Not sure about the ordering of the loops here
            for n in range(nn):
                temp[n] = i_s[ii,jj,n]
            for k in range(kmax-kmin):
                for n in range(nn):
                    o_s[i,j,k,n] = temp[n] + dx*(e[i,j,k,n] - temp[n]*a[i,j,k,n])
                    temp[n] = o_s[i,j,k,n]
            for n in range(nn):
                i_s[ii,jj,n] = temp[n]
    free(temp)

@cython.boundscheck(False)
def TransferShells(np.ndarray[np.float64_t, ndim=3] i_s,
                   np.ndarray[np.float64_t, ndim=3] data,
                   np.ndarray[np.float64_t, ndim=2] shells):
    """
    This function accepts an incoming slab (*i_s*), a buffer of *data*,
    and a list of shells specified as [ (value, tolerance, r, g, b), ... ].
    """
    cdef int i, ii
    cdef int j, jj
    cdef int k, kk
    cdef int n, nn
    cdef np.float64_t dist
    ii = data.shape[0]
    jj = data.shape[1]
    kk = data.shape[2]
    nn = shells.shape[0]
    cdef float rgba[4]
    cdef float alpha
    for i in range(ii):
        for j in range(jj):
            # Not sure about the ordering of the loops here
            for k in range(kk):
                for n in range(nn):
                    dist = shells[n, 0] - data[i,j,k]
                    if dist < 0: dist *= -1.0
                    if dist < shells[n,1]:
                        dist = exp(-dist/8.0)
                        rgba[0] = shells[n,2]
                        rgba[1] = shells[n,3]
                        rgba[2] = shells[n,4]
                        rgba[3] = shells[n,5]
                        alpha = i_s[i,j,3]
                        dist *= dist # This might improve appearance
                        i_s[i,j,0] += (1.0 - alpha)*rgba[0]*dist*rgba[3]
                        i_s[i,j,1] += (1.0 - alpha)*rgba[1]*dist*rgba[3]
                        i_s[i,j,2] += (1.0 - alpha)*rgba[2]*dist*rgba[3]
                        i_s[i,j,3] += (1.0 - alpha)*rgba[3]*dist*rgba[3]
                        break

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
                   np.ndarray[np.float64_t, ndim=3] grid_dt,
                   np.ndarray[np.float64_t, ndim=1] left_edge,
                   np.ndarray[np.float64_t, ndim=1] right_edge,
                   np.ndarray[np.float64_t, ndim=1] dx,
                   np.ndarray[np.float64_t, ndim=1] u,
                   np.ndarray[np.float64_t, ndim=1] v):
    # We're roughly following Amanatides & Woo
    # Find the first place the ray hits the grid on its path
    # Do left edge then right edge in each dim
    cdef int i, x, y
    cdef np.float64_t tl, tr, intersect_t, enter_t
    cdef np.float64_t iv_dir[3]
    cdef np.float64_t tdelta[3]
    cdef np.float64_t tmax[3]
    cdef np.float64_t intersect[3]
    cdef np.int64_t cur_ind[3]
    cdef np.int64_t step[3]
    intersect_t = 1
    # recall p = v * t + u
    #  where p is position, v is our vector, u is the start point
    for i in range(3):
        # As long as we're iterating, set some other stuff, too
        if(v[i] < 0):
            step[i] = -1
        elif (v[i] == 0):
            step[i] = 1
            tmax[i] = 1e60
            iv_dir[i] = 1e60
            tdelta[i] = 1e-60
            continue
        else:
            step[i] = 1
        x = (i+1)%3
        y = (i+2)%3
        iv_dir[i] = 1.0/v[i]
        tl = (left_edge[i] - u[i])*iv_dir[i]
        tr = (right_edge[i] - u[i])*iv_dir[i]
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
    for i in range(3):
        intersect[i] = u[i] + intersect_t * v[i]
        cur_ind[i] = np.floor((intersect[i] + 1e-8*dx[i] - left_edge[i])/dx[i])
        tmax[i] = (((cur_ind[i]+step[i])*dx[i])+left_edge[i]-u[i])*iv_dir[i]
        if cur_ind[i] == grid_mask.shape[i] and step[i] < 0:
            cur_ind[i] = grid_mask.shape[i] - 1
        if step[i] > 0: tmax[i] = (((cur_ind[i]+1)*dx[i])+left_edge[i]-u[i])*iv_dir[i]
        if step[i] < 0: tmax[i] = (((cur_ind[i]+0)*dx[i])+left_edge[i]-u[i])*iv_dir[i]
        tdelta[i] = (dx[i]*iv_dir[i])
        if tdelta[i] < 0: tdelta[i] *= -1
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
            grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = 1.0
            grid_dt[cur_ind[0], cur_ind[1], cur_ind[2]] = 1.0 - enter_t
            break
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = enter_t
                grid_dt[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[0] - enter_t
                enter_t = tmax[0]
                tmax[0] += tdelta[0]
                cur_ind[0] += step[0]
            else:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = enter_t
                grid_dt[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
        else:
            if tmax[1] < tmax[2]:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = enter_t
                grid_dt[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[1] - enter_t
                enter_t = tmax[1]
                tmax[1] += tdelta[1]
                cur_ind[1] += step[1]
            else:
                grid_t[cur_ind[0], cur_ind[1], cur_ind[2]] = enter_t
                grid_dt[cur_ind[0], cur_ind[1], cur_ind[2]] = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
    return

@cython.wraparound(False)
@cython.boundscheck(False)
def PlaneVoxelIntegration(np.ndarray[np.float64_t, ndim=1] left_edge,
                          np.ndarray[np.float64_t, ndim=1] right_edge,
                          np.ndarray[np.float64_t, ndim=1] dx,
                          np.ndarray[np.float64_t, ndim=2] ug,
                          np.ndarray[np.float64_t, ndim=1] v,
                          np.ndarray[np.float64_t, ndim=2] image,
                          np.ndarray[np.float64_t, ndim=3] data,
                          np.ndarray[np.float64_t, ndim=2] shells):
    # We're roughly following Amanatides & Woo on a ray-by-ray basis
    # Note that for now it's just shells, but this can and should be
    # generalized to transfer functions
    cdef int i, vi
    cdef int nv = ug.shape[0]
    cdef int nshells = shells.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] u = np.empty((3,), dtype=np.float64)
    # Copy things into temporary location for passing between functions
    for vi in range(nv):
        for i in range(3): u[i] = ug[vi, i]
        integrate_ray(u, v, left_edge, right_edge, dx,
                      nshells, vi, data, shells, image)

@cython.wraparound(False)
@cython.boundscheck(False)
def integrate_ray(np.ndarray[np.float64_t, ndim=1] u,
                  np.ndarray[np.float64_t, ndim=1] v,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.float64_t, ndim=1] dx,
                  int nshells, int ind,
                  np.ndarray[np.float64_t, ndim=3] data,
                  np.ndarray[np.float64_t, ndim=2] shells,
                  np.ndarray[np.float64_t, ndim=2] image):
    cdef int x, y, i, n
    cdef int step[3]
    cdef np.float64_t intersect_t = 1
    cdef np.float64_t tl, tr, enter_t
    cdef np.int64_t cur_ind[3]
    cdef np.float64_t tdelta[3]
    cdef np.float64_t tmax[3]
    cdef np.float64_t intersect[3]
    cdef np.float64_t dv
    cdef np.float64_t dist, alpha
    cdef int dims[3]
    cdef np.float64_t temp_x, temp_y
    for i in range(3):
        # As long as we're iterating, set some other stuff, too
        dims[i] = data.shape[i]
        if(v[i] < 0): step[i] = -1
        else: step[i] = 1
        x = (i+1)%3
        y = (i+2)%3
        tl = (left_edge[i] - u[i])/v[i]
        tr = (right_edge[i] - u[i])/v[i]
        temp_x = (u[x] + tl*v[x])
        temp_y = (u[y] + tl*v[y])
        if (left_edge[x] <= temp_x) and (temp_x <= right_edge[x]) and \
           (left_edge[y] <= temp_y) and (temp_y <= right_edge[y]) and \
           (0.0 <= tl) and (tl < intersect_t):
            intersect_t = tl
        temp_x = (u[x] + tr*v[x])
        temp_y = (u[y] + tr*v[y])
        if (left_edge[x] <= temp_x) and (temp_x <= right_edge[x]) and \
           (left_edge[y] <= temp_y) and (temp_y <= right_edge[y]) and \
           (0.0 <= tr) and (tr < intersect_t):
            intersect_t = tr
    # if fully enclosed
    if (left_edge[0] <= u[0] <= right_edge[0]) and \
       (left_edge[1] <= u[1] <= right_edge[1]) and \
       (left_edge[2] <= u[2] <= right_edge[2]):
        intersect_t = 0.0
    if not (0 <= intersect_t <= 1):
        #print "Returning: intersect_t ==", intersect_t
        return
    # Now get the indices of the intersection
    for i in range(3): intersect[i] = u[i] + intersect_t * v[i]
    for i in range(3):
        cur_ind[i] = np.floor((intersect[i] + 1e-8*dx[i] - left_edge[i])/dx[i])
        tmax[i] = (((cur_ind[i]+step[i])*dx[i])+left_edge[i]-u[i])/v[i]
        if cur_ind[i] == dims[i] and step[i] < 0:
            cur_ind[i] = dims[i] - 1
        if step[i] > 0: tmax[i] = (((cur_ind[i]+1)*dx[i])+left_edge[i]-u[i])/v[i]
        if step[i] < 0: tmax[i] = (((cur_ind[i]+0)*dx[i])+left_edge[i]-u[i])/v[i]
        tdelta[i] = (dx[i]/v[i])
        if tdelta[i] < 0: tdelta[i] *= -1
    # The variable intersect contains the point we first pierce the grid
    enter_t = intersect_t
    if (not (0 <= cur_ind[0] < dims[0])) or \
       (not (0 <= cur_ind[1] < dims[1])) or \
       (not (0 <= cur_ind[2] < dims[2])):
        #print "Returning: cur_ind", cur_ind[0], cur_ind[1], cur_ind[2]
        #print "  dims:     ", dims[0], dims[1], dims[2]
        #print "  intersect:",  intersect[0], intersect[1], intersect[2]
        #print "  intersect:", intersect_t
        #print "  u        :", u[0], u[1], u[2]
        #
        return
    #print cur_ind[0], dims[0], cur_ind[1], dims[1], cur_ind[2], dims[2]
    dv = data[cur_ind[0], cur_ind[1], cur_ind[2]]
    #dt = 1e300
    while 1:
        if image[ind,3] >= 1.0: break
        if (not (0 <= cur_ind[0] < dims[0])) or \
           (not (0 <= cur_ind[1] < dims[1])) or \
           (not (0 <= cur_ind[2] < dims[2])):
            break
        # Do our transfer here
        for n in range(nshells):
            dist = shells[n, 0] - dv
            if dist < shells[n,1]:
                dist = exp(-dist/8.0)
                alpha = (1.0 - shells[n,5])*shells[n,5]#*dt
                image[ind,0] += alpha*shells[n,2]*dist
                image[ind,1] += alpha*shells[n,3]*dist
                image[ind,2] += alpha*shells[n,4]*dist
                image[ind,3] += alpha*shells[n,5]*dist
                #image[ind,i] += rgba[i]*dist*rgba[3]/dt
                #print rgba[i], image[ind,i], dist, dt
                break
        if (tmax[0] > 1.0) and (tmax[1] > 1.0) and (tmax[2] > 1.0):
            dt = 1.0 - enter_t
            break
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                dt = tmax[0] - enter_t
                enter_t = tmax[0]
                tmax[0] += tdelta[0]
                cur_ind[0] += step[0]
            else:
                dt = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
        else:
            if tmax[1] < tmax[2]:
                dt = tmax[1] - enter_t
                enter_t = tmax[1]
                tmax[1] += tdelta[1]
                cur_ind[1] += step[1]
            else:
                dt = tmax[2] - enter_t
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
                cur_ind[2] += step[2]
        dv = data[cur_ind[0], cur_ind[1], cur_ind[2]]
