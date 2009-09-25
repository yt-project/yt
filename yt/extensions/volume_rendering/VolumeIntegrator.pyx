"""
Simple integrators for the radiative transfer equation

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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
from stdlib cimport malloc, free, abs

cdef inline int imax(int i0, int i1):
    if i0 > i1: return i0
    return i1

cdef inline np.float64_t fmax(np.float64_t f0, np.float64_t f1):
    if f0 > f1: return f0
    return f1

cdef inline int imin(int i0, int i1):
    if i0 < i1: return i0
    return i1

cdef inline np.float64_t fmin(np.float64_t f0, np.float64_t f1):
    if f0 < f1: return f0
    return f1

cdef inline int iclip(int i, int a, int b):
    return imin(imax(i, a), b)

cdef inline np.float64_t fclip(np.float64_t f,
                      np.float64_t a, np.float64_t b):
    return fmin(fmax(f, a), b)

cdef extern from "math.h":
    double exp(double x)
    float expf(float x)
    double floor(double x)
    double ceil(double x)
    double fmod(double x, double y)

cdef extern from "FixedInterpolator.h":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data)
cdef class VectorPlane

cdef class TransferFunctionProxy:
    cdef np.float64_t x_bounds[2]
    cdef np.float64_t *vs[3]
    cdef int nbins
    cdef np.float64_t dbin
    cdef public object tf_obj
    cdef public object ci
    cdef VectorPlane vp
    def __cinit__(self, tf_obj):
        self.tf_obj = tf_obj
        cdef np.ndarray[np.float64_t, ndim=1] temp
        temp = tf_obj.red.y
        self.vs[0] = <np.float64_t *> temp.data
        temp = tf_obj.green.y
        self.vs[1] = <np.float64_t *> temp.data
        temp = tf_obj.blue.y
        self.vs[2] = <np.float64_t *> temp.data
        self.x_bounds[0] = tf_obj.x_bounds[0]
        self.x_bounds[1] = tf_obj.x_bounds[1]
        self.nbins = tf_obj.nbins
        self.dbin = (self.x_bounds[1] - self.x_bounds[0])/self.nbins

    cdef void eval_transfer(self,
                       np.float64_t dv0,
                       np.float64_t dv1,
                       np.float64_t dt,
                       np.float64_t rgba[4]):
        cdef int i
        cdef int b0, b1
        cdef np.float64_t tf0, tf1
        rgba[0] += dt
        return
        for i in range(3):
            # First locate our points
            b0 = iclip(<int> floor((dv0 - self.x_bounds[0]) / self.dbin),
                       0, self.nbins-1)
            b1 = b0 + 1
            tf0 = (dv0 - self.x_bounds[0] - self.dbin * b0) * self.vs[i][b0] + \
                  (self.x_bounds[0] + self.dbin * b1 - dv0) * self.vs[i][b1]

            b0 = iclip(<int> floor((dv1 - self.x_bounds[0]) / self.dbin),
                       0, self.nbins-1)
            b1 = b0 + 1
            tf1 = (dv1 - self.x_bounds[0] - self.dbin * b0) * self.vs[i][b0] + \
                  (self.x_bounds[0] + self.dbin * b1 - dv1) * self.vs[i][b1]
            # alpha blending goes here
            rgba[i] += dt*0.5*(tf0+tf1)
        # We should update alpha here

cdef class VectorPlane:
    cdef public object avp_pos, avp_dir, acenter, aimage
    cdef np.float64_t *vp_pos, *vp_dir, *center, *image,
    cdef np.float64_t pdx, pdy, bounds[4]
    cdef int nv
    cdef public object ax_vec, ay_vec
    cdef np.float64_t *x_vec, *y_vec
    cdef public object xpos, ypos, zpos

    def __cinit__(self, 
                  np.ndarray[np.float64_t, ndim=3] vp_pos,
                  np.ndarray[np.float64_t, ndim=1] vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec):
        cdef int i, j
        self.avp_pos = vp_pos
        self.avp_dir = vp_dir
        self.acenter = center
        self.aimage = image
        self.ax_vec = x_vec
        self.ay_vec = y_vec
        self.vp_pos = <np.float64_t *> vp_pos.data
        self.vp_dir = <np.float64_t *> vp_dir.data
        self.center = <np.float64_t *> center.data
        self.image = <np.float64_t *> image.data
        self.x_vec = <np.float64_t *> x_vec.data
        self.y_vec = <np.float64_t *> y_vec.data
        self.nv = vp_pos.shape[0]
        for i in range(4): self.bounds[i] = bounds[i]
        self.pdx = (self.bounds[1] - self.bounds[0])/self.nv
        self.pdy = (self.bounds[3] - self.bounds[2])/self.nv
        self.xpos = []
        self.ypos = []
        self.zpos = []
        for i in range(self.nv):
            self.xpos.append([])
            self.ypos.append([])
            self.zpos.append([])
            for j in range(self.nv):
                self.xpos[-1].append([])
                self.ypos[-1].append([])
                self.zpos[-1].append([])

    cdef void get_start_stop(self, np.float64_t *ex, int *rv):
        rv[0] = iclip(<int> floor((ex[0] - self.bounds[0])/self.pdx), 0, self.nv)
        rv[1] = iclip(<int> floor((ex[2] - self.bounds[2])/self.pdy), 0, self.nv)
        rv[2] = iclip(rv[0] + <int> ceil((ex[1] - ex[0])/self.pdx), 0, self.nv)
        rv[3] = iclip(rv[1] + <int> ceil((ex[3] - ex[2])/self.pdy), 0, self.nv)

    cdef void copy(self, np.float64_t *fv, np.float64_t *tv,
                   int nk, int i, int j):
        # We know the first two dimensions of our from-vector, and our
        # to-vector is flat and 'ni' long
        cdef int k
        for k in range(nk):
            tv[k] = fv[(i*self.nv+j)*self.nv+k]

    cdef void copy_back(self, np.float64_t *fv, np.float64_t *tv,
                       int nk, int i, int j):
        cdef int k
        for k in range(nk):
            tv[(i*self.nv+j)*self.nv+k] = fv[k] 

cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef np.float64_t *data
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef np.float64_t dds[3]
    cdef int dims[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self,
                  np.ndarray[np.float64_t, ndim=3] data,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.int64_t, ndim=1] dims):
        # The data is likely brought in via a slice, so we copy it
        cdef int i
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        cdef np.ndarray[np.float64_t, ndim=3] tdata = data.copy()
        self.my_data = tdata
        self.data = <np.float64_t *> (tdata.data)
        
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
            self.dims[i] = dims[i]
            self.dds[i] = (self.right_edge[i] - self.left_edge[i])/dims[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def cast_plane(self, TransferFunctionProxy tf, VectorPlane vp):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, ni, nj, nn
        cdef int iter[4]
        cdef np.float64_t v_pos[3], v_dir[3], rgba[4], extrema[4]
        tf.vp = vp
        self.calculate_extent(vp, extrema)
        vp.get_start_stop(extrema, iter)
        hit = 0
        for vi in range(vp.nv):#iter[0], iter[2]):
            for vj in range(vp.nv):#iter[1], iter[3]):
                #vp.copy(vp.vp_pos, v_pos, vi, vj, 3)
                #vp.copy(vp.image, rgba, vi, vj, 4)
                for i in range(3): v_pos[i] = vp.avp_pos[vi, vj, i]
                for i in range(4): rgba[i] = vp.aimage[vi, vj, i]
                tf.ci = (vi, vj)
                hit += self.integrate_ray(v_pos, vp.vp_dir, rgba, tf)
                #vp.copy_back(rgba, vp.image, vi, vj, 4)
                for i in range(4): vp.aimage[vi,vj,i] = rgba[i]
        return (tf.vp.xpos, tf.vp.ypos, tf.vp.zpos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void calculate_extent(self, VectorPlane vp,
                               np.float64_t extrema[4]):
        # We do this for all eight corners
        cdef np.float64_t *edges[2], temp
        edges[0] = self.left_edge
        edges[1] = self.right_edge
        extrema[0] = extrema[2] = 1e300; extrema[1] = extrema[3] = -1e300
        cdef int i, j, k
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    temp  = edges[i][0] * vp.x_vec[0]
                    temp += edges[j][1] * vp.x_vec[1]
                    temp += edges[k][2] * vp.x_vec[2]
                    if temp < extrema[0]: extrema[0] = temp
                    if temp > extrema[1]: extrema[1] = temp
                    temp  = edges[i][0] * vp.y_vec[0]
                    temp += edges[j][1] * vp.y_vec[1]
                    temp += edges[k][2] * vp.y_vec[2]
                    if temp < extrema[2]: extrema[2] = temp
                    if temp > extrema[3]: extrema[3] = temp
        #print extrema[0], extrema[1], extrema[2], extrema[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int integrate_ray(self, np.float64_t v_pos[3],
                                 np.float64_t v_dir[3],
                                 np.float64_t rgba[4],
                                 TransferFunctionProxy tf):
        cdef int cur_ind[3], step[3], x, y, i, n, flat_ind, hit
        cdef np.float64_t intersect_t = 1.0
        cdef np.float64_t intersect[3], tmax[3], tdelta[3]
        cdef np.float64_t enter_t, dist, alpha, dt
        cdef np.float64_t tr, tl, temp_x, temp_y, dv
        for i in range(3):
            if (v_dir[i] < 0):
                step[i] = -1
            else:
                step[i] = 1
            x = (i+1) % 3
            y = (i+2) % 3
            tl = (self.left_edge[i] - v_pos[i])/v_dir[i]
            tr = (self.right_edge[i] - v_pos[i])/v_dir[i]
            temp_x = (v_pos[x] + tl*v_dir[x])
            temp_y = (v_pos[y] + tl*v_dir[y])
            if self.left_edge[x] <= temp_x and temp_x <= self.right_edge[x] and \
               self.left_edge[y] <= temp_y and temp_y <= self.right_edge[y] and \
               0.0 <= tl and tl < intersect_t:
                intersect_t = tl
            temp_x = (v_pos[x] + tr*v_dir[x])
            temp_y = (v_pos[y] + tr*v_dir[y])
            if self.left_edge[x] <= temp_x and temp_x <= self.right_edge[x] and \
               self.left_edge[y] <= temp_y and temp_y <= self.right_edge[y] and \
               0.0 <= tr and tr < intersect_t:
                intersect_t = tr
        if self.left_edge[0] <= v_pos[0] and v_pos[0] <= self.right_edge[0] and \
           self.left_edge[1] <= v_pos[1] and v_pos[1] <= self.right_edge[1] and \
           self.left_edge[2] <= v_pos[2] and v_pos[2] <= self.right_edge[2]:
            intersect_t = 0.0
        if not ((0.0 <= intersect_t) and (intersect_t < 1.0)): return 0
        for i in range(3):
            intersect[i] = v_pos[i] + intersect_t * v_dir[i]
            cur_ind[i] = <int> floor((intersect[i] + 1e-8*self.dds[i] -
                                      self.left_edge[i])/self.dds[i])
            tmax[i] = (((cur_ind[i]+step[i])*self.dds[i])+
                        self.left_edge[i]-v_pos[i])/v_dir[i]
            if cur_ind[i] == self.dims[i] and step[i] < 0:
                cur_ind[i] = self.dims[i] - 1
            if cur_ind[i] < 0 or cur_ind[i] >= self.dims[i]: return 0
            if step[i] > 0:
                tmax[i] = (((cur_ind[i]+1)*self.dds[i])
                            +self.left_edge[i]-v_pos[i])/v_dir[i]
            if step[i] < 0:
                tmax[i] = (((cur_ind[i]+0)*self.dds[i])
                            +self.left_edge[i]-v_pos[i])/v_dir[i]
            tdelta[i] = (self.dds[i]/v_dir[i])
            if tdelta[i] < 0: tdelta[i] *= -1
        # We have to jumpstart our calculation
        enter_t = intersect_t
        dv0 = self.get_val(v_pos, v_dir, enter_t, tf)
        while 1:
            if (not (0 <= cur_ind[0] < self.dims[0])) or \
               (not (0 <= cur_ind[1] < self.dims[1])) or \
               (not (0 <= cur_ind[2] < self.dims[2])):
                break
            hit += 1
            # Do our transfer here
            if tmax[0] < tmax[1]:
                if tmax[0] < tmax[2]:
                    cur_ind[0] += step[0]
                    dv1 = self.get_val(v_pos, v_dir, tmax[0], tf)
                    dt = fmin(tmax[0], 1.0) - enter_t
                    dv0 = dv1
                    enter_t = tmax[0]
                    tmax[0] += tdelta[0]
                else:
                    cur_ind[2] += step[2]
                    dv1 = self.get_val(v_pos, v_dir, tmax[2], tf)
                    dt = fmin(tmax[2], 1.0) - enter_t
                    dv0 = dv1
                    enter_t = tmax[2]
                    tmax[2] += tdelta[2]
            else:
                if tmax[1] < tmax[2]:
                    cur_ind[1] += step[1]
                    dv1 = self.get_val(v_pos, v_dir, tmax[1], tf)
                    dt = fmin(tmax[1], 1.0) - enter_t
                    dv0 = dv1
                    enter_t = tmax[1]
                    tmax[1] += tdelta[1]
                else:
                    cur_ind[2] += step[2]
                    dv1 = self.get_val(v_pos, v_dir, tmax[2], tf)
                    dt = fmin(tmax[2], 1.0) - enter_t
                    dv0 = dv1
                    enter_t = tmax[2]
                    tmax[2] += tdelta[2]
            tf.eval_transfer(dv0, dv1, dt, rgba)
            if enter_t > 1.0: break
        return hit

    cdef np.float64_t get_val(self,
                              np.float64_t v_pos[3],
                              np.float64_t v_dir[3],
                              np.float64_t t,
                              TransferFunctionProxy tf):
        cdef np.float64_t cp[3], dv, dp[3], temp
        cdef int i, ci[3]
        t = fclip(t, 0.0, 1.0)
        for i in range(3):
            cp[i] = v_pos[i] + t * v_dir[i]
            temp = v_pos[i] - self.left_edge[i] + 1e-8*self.dds[i]
            ci[i] = <int> floor( temp / self.dds[i] )
            dp[i] = fmod(temp, self.dds[i])/self.dds[i]
        tf.vp.xpos[tf.ci[0]][tf.ci[1]].append(cp[0])
        tf.vp.ypos[tf.ci[0]][tf.ci[1]].append(cp[1])
        tf.vp.zpos[tf.ci[0]][tf.ci[1]].append(cp[2])
        fast_interpolate(self.dims, ci, dp, self.data)
