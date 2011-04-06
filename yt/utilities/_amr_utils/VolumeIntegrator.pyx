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
cimport kdtree_utils
cimport healpix_interface
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
    if i < a: return a
    if i > b: return b
    return i

cdef inline np.float64_t fclip(np.float64_t f,
                      np.float64_t a, np.float64_t b):
    return fmin(fmax(f, a), b)

cdef extern from "math.h":
    double exp(double x)
    float expf(float x)
    double floor(double x)
    double ceil(double x)
    double fmod(double x, double y)
    double log2(double x)
    long int lrint(double x)
    double fabs(double x)

cdef extern from "FixedInterpolator.h":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data)
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3], np.float64_t *data)
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data)
    np.float64_t eval_gradient(int *ds, int *ci, np.float64_t *dp,
                                       np.float64_t *data, np.float64_t *grad)

def hp_pix2vec_nest(long nside, long ipix):
    cdef double v[3]
    healpix_interface.pix2vec_nest(nside, ipix, v)
    cdef np.ndarray[np.float64_t, ndim=1] tr = np.empty((3,), dtype='float64')
    tr[0] = v[0]
    tr[1] = v[1]
    tr[2] = v[2]
    return tr

def arr_pix2vec_nest(long nside,
                     np.ndarray[np.int64_t, ndim=1] aipix):
    cdef int n = aipix.shape[0]
    cdef int i
    cdef double v[3]
    cdef long ipix
    cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 3), dtype='float64')
    for i in range(n):
        ipix = aipix[i]
        healpix_interface.pix2vec_nest(nside, ipix, v)
        tr[i,0] = v[0]
        tr[i,1] = v[1]
        tr[i,2] = v[2]
    return tr

def hp_vec2pix_nest(long nside, double x, double y, double z):
    cdef double v[3]
    v[0] = x
    v[1] = y
    v[2] = z
    cdef long ipix
    healpix_interface.vec2pix_nest(nside, v, &ipix)
    return ipix

def arr_vec2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] x,
                     np.ndarray[np.float64_t, ndim=1] y,
                     np.ndarray[np.float64_t, ndim=1] z):
    cdef int n = x.shape[0]
    cdef int i
    cdef double v[3]
    cdef long ipix
    cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    for i in range(n):
        v[0] = x[i]
        v[1] = y[i]
        v[2] = z[i]
        healpix_interface.vec2pix_nest(nside, v, &ipix)
        tr[i] = ipix
    return tr

def hp_pix2ang_nest(long nside, long ipnest):
    cdef double theta, phi
    healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
    return (theta, phi)

def arr_pix2ang_nest(long nside, np.ndarray[np.int64_t, ndim=1] aipnest):
    cdef int n = aipnest.shape[0]
    cdef int i
    cdef long ipnest
    cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 2), dtype='float64')
    cdef double theta, phi
    for i in range(n):
        ipnest = aipnest[i]
        healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
        tr[i,0] = theta
        tr[i,1] = phi
    return tr

def hp_ang2pix_nest(long nside, double theta, double phi):
    cdef long ipix
    healpix_interface.ang2pix_nest(nside, theta, phi, &ipix)
    return ipix

def arr_ang2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] atheta,
                     np.ndarray[np.float64_t, ndim=1] aphi):
    cdef int n = atheta.shape[0]
    cdef int i
    cdef long ipnest
    cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    cdef double theta, phi
    for i in range(n):
        theta = atheta[i]
        phi = aphi[i]
        healpix_interface.ang2pix_nest(nside, theta, phi, &ipnest)
        tr[i] = ipnest
    return tr

cdef class star_kdtree_container:
    cdef kdtree_utils.kdtree *tree
    cdef public np.float64_t sigma
    cdef public np.float64_t coeff

    def __init__(self):
        self.tree = kdtree_utils.kd_create(3)

    def add_points(self,
                   np.ndarray[np.float64_t, ndim=1] pos_x,
                   np.ndarray[np.float64_t, ndim=1] pos_y,
                   np.ndarray[np.float64_t, ndim=1] pos_z,
                   np.ndarray[np.float64_t, ndim=2] star_colors):
        cdef int i, n
        cdef np.float64_t *pointer = <np.float64_t *> star_colors.data
        for i in range(pos_x.shape[0]):
            kdtree_utils.kd_insert3(self.tree,
                pos_x[i], pos_y[i], pos_z[i], pointer + i*3)

    def __dealloc__(self):
        kdtree_utils.kd_free(self.tree)

cdef class VectorPlane

cdef struct FieldInterpolationTable:
    # Note that we make an assumption about retaining a reference to values
    # externally.
    np.float64_t *values 
    np.float64_t bounds[2]
    np.float64_t dbin
    np.float64_t idbin
    int field_id
    int weight_field_id
    int weight_table_id
    int nbins
    int pass_through

cdef void FIT_initialize_table(FieldInterpolationTable *fit, int nbins,
              np.float64_t *values, np.float64_t bounds1, np.float64_t bounds2,
              int field_id, int weight_field_id = -1, int weight_table_id = -1,
              int pass_through = 0):
    fit.bounds[0] = bounds1; fit.bounds[1] = bounds2
    fit.nbins = nbins
    fit.dbin = (fit.bounds[1] - fit.bounds[0])/fit.nbins
    fit.idbin = 1.0/fit.dbin
    # Better not pull this out from under us, yo
    fit.values = values
    fit.field_id = field_id
    fit.weight_field_id = weight_field_id
    fit.weight_table_id = weight_table_id
    fit.pass_through = pass_through

cdef np.float64_t FIT_get_value(FieldInterpolationTable *fit,
                            np.float64_t *dvs):
    cdef np.float64_t bv, dy, dd, tf
    cdef int bin_id
    if fit.pass_through == 1: return dvs[fit.field_id]
    if dvs[fit.field_id] > fit.bounds[1] or dvs[fit.field_id] < fit.bounds[0]: return 0.0
    bin_id = <int> ((dvs[fit.field_id] - fit.bounds[0]) * fit.idbin)
    dd = dvs[fit.field_id] - (fit.bounds[0] + bin_id * fit.dbin) # x - x0
    bv = fit.values[bin_id]
    dy = fit.values[bin_id + 1] - bv
    if fit.weight_field_id != -1:
        return dvs[fit.weight_field_id] * (bv + dd*dy*fit.idbin)
    return (bv + dd*dy*fit.idbin)

cdef class TransferFunctionProxy:
    cdef int n_fields
    cdef int n_field_tables
    cdef public int ns

    # These are the field tables and their affiliated storage.
    # We have one field_id for every table.  Note that a single field can
    # correspond to multiple tables, and each field table will only have
    # interpolate called once.
    cdef FieldInterpolationTable field_tables[6]
    cdef np.float64_t istorage[6]

    # Here are the field tables that correspond to each of the six channels.
    # We have three emission channels, three absorption channels.
    cdef int field_table_ids[6]

    # We store a reference to the transfer function object and to the field
    # interpolation tables
    cdef public object tf_obj
    cdef public object my_field_tables

    def __cinit__(self, tf_obj):
        # We have N fields.  We have 6 channels.  We have M field tables.
        # The idea is that we can have multiple channels corresponding to the
        # same field table.  So, we create storage for the outputs from all the
        # field tables.  We need to know which field value to pass in to the
        # field table, and we need to know which table to use for each of the
        # six channels.
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] temp
        cdef FieldInterpolationTable fit

        self.tf_obj = tf_obj

        self.n_field_tables = tf_obj.n_field_tables
        for i in range(6): self.istorage[i] = 0.0

        self.my_field_tables = []
        for i in range(self.n_field_tables):
            temp = tf_obj.tables[i].y
            FIT_initialize_table(&self.field_tables[i],
                      temp.shape[0],
                      <np.float64_t *> temp.data,
                      tf_obj.tables[i].x_bounds[0],
                      tf_obj.tables[i].x_bounds[1],
                      tf_obj.field_ids[i], tf_obj.weight_field_ids[i],
                      tf_obj.weight_table_ids[i],
                      tf_obj.tables[i].pass_through)
            self.my_field_tables.append((tf_obj.tables[i],
                                         tf_obj.tables[i].y))
            self.field_tables[i].field_id = tf_obj.field_ids[i]
            self.field_tables[i].weight_field_id = tf_obj.weight_field_ids[i]
            print "Field table", i, "corresponds to",
            print self.field_tables[i].field_id,
            print "(Weighted with ", self.field_tables[i].weight_field_id,
            print ")"

        for i in range(6):
            self.field_table_ids[i] = tf_obj.field_table_ids[i]
            print "Channel", i, "corresponds to", self.field_table_ids[i]
            
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void eval_transfer(self, np.float64_t dt, np.float64_t *dvs,
                                  np.float64_t *rgba, np.float64_t *grad):
        cdef int i, fid, use
        cdef np.float64_t ta, tf, trgba[6], dot_prod
        # NOTE: We now disable this.  I have left it to ease the process of
        # potentially, one day, re-including it.
        #use = 0
        #for i in range(self.n_field_tables):
        #    fid = self.field_tables[i].field_id
        #    if (dvs[fid] >= self.field_tables[i].bounds[0]) and \
        #       (dvs[fid] <= self.field_tables[i].bounds[1]):
        #        use = 1
        #        break
        for i in range(self.n_field_tables):
            self.istorage[i] = FIT_get_value(&self.field_tables[i], dvs)
        # We have to do this after the interpolation
        for i in range(self.n_field_tables):
            fid = self.field_tables[i].weight_table_id
            if fid != -1: self.istorage[i] *= self.istorage[fid]
        for i in range(6):
            trgba[i] = self.istorage[self.field_table_ids[i]]
            #print i, trgba[i],
        #print
        # A few words on opacity.  We're going to be integrating equation 1.23
        # from Rybicki & Lightman.  dI_\nu / ds = -\alpha_\nu I_\nu + j_\nu
        # \alpha_nu = \kappa \rho , but we leave that up to the input
        # transfer function.
        # SOoooooOOOooo, the upshot is that we are doing a rectangular
        # integration here:
        #   I_{i+1} = ds * C_i + (1.0 - ds*alpha_i) * I_i
        for i in range(3):
            # This is the new way: alpha corresponds to opacity of a given
            # slice.  Previously it was ill-defined, but represented some
            # measure of emissivity.
            ta = fmax((1.0 - dt*trgba[i+3]), 0.0)
            rgba[i  ] = dt*trgba[i  ] + ta * rgba[i  ]
            #rgba[i+3] = dt*trgba[i+3] + ta * rgba[i+3]
            # This is the old way:
            #rgba[i  ] += trgba[i] * (1.0 - rgba[i+3])*dt*trgba[i+3]
            #rgba[i+3] += trgba[i] * (1.0 - rgba[i+3])*dt*trgba[i+3]

cdef class VectorPlane:
    cdef public object avp_pos, avp_dir, acenter, aimage
    cdef np.float64_t *vp_pos, *vp_dir, *center, *image,
    cdef np.float64_t pdx, pdy, bounds[4]
    cdef int nv[2]
    cdef int vp_strides[3]
    cdef int im_strides[3]
    cdef int vd_strides[3]
    cdef public object ax_vec, ay_vec
    cdef np.float64_t *x_vec, *y_vec

    def __cinit__(self, 
                  np.ndarray[np.float64_t, ndim=3] vp_pos,
                  np.ndarray vp_dir,
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
        self.nv[0] = vp_pos.shape[0]
        self.nv[1] = vp_pos.shape[1]
        for i in range(4): self.bounds[i] = bounds[i]
        self.pdx = (self.bounds[1] - self.bounds[0])/self.nv[0]
        self.pdy = (self.bounds[3] - self.bounds[2])/self.nv[1]
        for i in range(3):
            self.vp_strides[i] = vp_pos.strides[i] / 8
            self.im_strides[i] = image.strides[i] / 8
        if vp_dir.ndim > 1:
            for i in range(3):
                self.vd_strides[i] = vp_dir.strides[i] / 8
        else:
            self.vd_strides[0] = self.vd_strides[1] = self.vd_strides[2] = -1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void get_start_stop(self, np.float64_t *ex, int *rv):
        # Extrema need to be re-centered
        cdef np.float64_t cx, cy
        cdef int i
        cx = cy = 0.0
        for i in range(3):
            cx += self.center[i] * self.x_vec[i]
            cy += self.center[i] * self.y_vec[i]
        rv[0] = lrint((ex[0] - cx - self.bounds[0])/self.pdx)
        rv[1] = rv[0] + lrint((ex[1] - ex[0])/self.pdx)
        rv[2] = lrint((ex[2] - cy - self.bounds[2])/self.pdy)
        rv[3] = rv[2] + lrint((ex[3] - ex[2])/self.pdy)

    cdef inline void copy_into(self, np.float64_t *fv, np.float64_t *tv,
                        int i, int j, int nk, int strides[3]):
        # We know the first two dimensions of our from-vector, and our
        # to-vector is flat and 'ni' long
        cdef int k
        cdef int offset = strides[0] * i + strides[1] * j
        for k in range(nk):
            tv[k] = fv[offset + k]

    cdef inline void copy_back(self, np.float64_t *fv, np.float64_t *tv,
                        int i, int j, int nk, int strides[3]):
        cdef int k
        cdef int offset = strides[0] * i + strides[1] * j
        for k in range(nk):
            tv[offset + k] = fv[k]

cdef struct AdaptiveRayPacket

cdef class PartitionedGrid:
    cdef public object my_data
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef np.float64_t *data[6]
    cdef np.float64_t dvs[6]
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef np.float64_t dds[3]
    cdef np.float64_t idds[3]
    cdef int dims[3]
    cdef public int parent_grid_id
    cdef public int n_fields
    cdef kdtree_utils.kdtree *star_list
    cdef np.float64_t star_er
    cdef np.float64_t star_sigma_num
    cdef np.float64_t star_coeff

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self,
                  int parent_grid_id, int n_fields, data,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  np.ndarray[np.int64_t, ndim=1] dims,
                  star_kdtree_container star_tree = None):
        # The data is likely brought in via a slice, so we copy it
        cdef int i, j, k, size
        cdef np.ndarray[np.float64_t, ndim=3] tdata
        self.parent_grid_id = parent_grid_id
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
            self.dims[i] = dims[i]
            self.dds[i] = (self.right_edge[i] - self.left_edge[i])/dims[i]
            self.idds[i] = 1.0/self.dds[i]
        self.my_data = data
        self.n_fields = n_fields
        for i in range(n_fields):
            tdata = data[i]
            self.data[i] = <np.float64_t *> tdata.data
        if star_tree is None:
            self.star_list = NULL
        else:
            self.set_star_tree(star_tree)

    def set_star_tree(self, star_kdtree_container star_tree):
        self.star_list = star_tree.tree
        self.star_sigma_num = 2.0*star_tree.sigma**2.0
        self.star_er = 2.326 * star_tree.sigma
        self.star_coeff = star_tree.coeff

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def cast_plane(self, TransferFunctionProxy tf, VectorPlane vp):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, ni, nj, nn
        cdef int iter[4]
        cdef np.float64_t v_pos[3], v_dir[3], rgba[6], extrema[4]
        hit = 0
        self.calculate_extent(vp, extrema)
        vp.get_start_stop(extrema, iter)
        iter[0] = iclip(iter[0]-1, 0, vp.nv[0])
        iter[1] = iclip(iter[1]+1, 0, vp.nv[0])
        iter[2] = iclip(iter[2]-1, 0, vp.nv[1])
        iter[3] = iclip(iter[3]+1, 0, vp.nv[1])
        if vp.vd_strides[0] == -1:
            for vi in range(iter[0], iter[1]):
                for vj in range(iter[2], iter[3]):
                    vp.copy_into(vp.vp_pos, v_pos, vi, vj, 3, vp.vp_strides)
                    vp.copy_into(vp.image, rgba, vi, vj, 3, vp.im_strides)
                    self.integrate_ray(v_pos, vp.vp_dir, rgba, tf)
                    vp.copy_back(rgba, vp.image, vi, vj, 3, vp.im_strides)
        else:
            # If we do not have an orthographic projection, we have to cast all
            # our rays (until we can get an extrema calculation...)
            for vi in range(vp.nv[0]):
                for vj in range(vp.nv[1]):
                    vp.copy_into(vp.vp_pos, v_pos, vi, vj, 3, vp.vp_strides)
                    vp.copy_into(vp.image, rgba, vi, vj, 3, vp.im_strides)
                    vp.copy_into(vp.vp_dir, v_dir, vi, vj, 3, vp.vd_strides)
                    self.integrate_ray(v_pos, v_dir, rgba, tf)
                    vp.copy_back(rgba, vp.image, vi, vj, 3, vp.im_strides)
        return hit

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
                    # This should rotate it into the vector plane
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
                                 TransferFunctionProxy tf,
                                 np.float64_t *return_t = NULL):
        cdef int cur_ind[3], step[3], x, y, i, n, flat_ind, hit, direction
        cdef np.float64_t intersect_t = 1.0
        cdef np.float64_t iv_dir[3]
        cdef np.float64_t intersect[3], tmax[3], tdelta[3]
        cdef np.float64_t enter_t, dist, alpha, dt, exit_t
        cdef np.float64_t tr, tl, temp_x, temp_y, dv
        for i in range(3):
            if (v_dir[i] < 0):
                step[i] = -1
            elif (v_dir[i] == 0):
                step[i] = 1
                tmax[i] = 1e60
                iv_dir[i] = 1e60
                tdelta[i] = 1e-60
                continue
            else:
                step[i] = 1
            x = (i+1) % 3
            y = (i+2) % 3
            iv_dir[i] = 1.0/v_dir[i]
            tl = (self.left_edge[i] - v_pos[i])*iv_dir[i]
            temp_x = (v_pos[x] + tl*v_dir[x])
            temp_y = (v_pos[y] + tl*v_dir[y])
            if self.left_edge[x] <= temp_x and temp_x <= self.right_edge[x] and \
               self.left_edge[y] <= temp_y and temp_y <= self.right_edge[y] and \
               0.0 <= tl and tl < intersect_t:
                direction = i
                intersect_t = tl
            tr = (self.right_edge[i] - v_pos[i])*iv_dir[i]
            temp_x = (v_pos[x] + tr*v_dir[x])
            temp_y = (v_pos[y] + tr*v_dir[y])
            if self.left_edge[x] <= temp_x and temp_x <= self.right_edge[x] and \
               self.left_edge[y] <= temp_y and temp_y <= self.right_edge[y] and \
               0.0 <= tr and tr < intersect_t:
                direction = i
                intersect_t = tr
        if self.left_edge[0] <= v_pos[0] and v_pos[0] <= self.right_edge[0] and \
           self.left_edge[1] <= v_pos[1] and v_pos[1] <= self.right_edge[1] and \
           self.left_edge[2] <= v_pos[2] and v_pos[2] <= self.right_edge[2]:
            intersect_t = 0.0
        if not ((0.0 <= intersect_t) and (intersect_t < 1.0)): return 0
        for i in range(3):
            intersect[i] = v_pos[i] + intersect_t * v_dir[i]
            cur_ind[i] = <int> floor((intersect[i] +
                                      step[i]*1e-8*self.dds[i] -
                                      self.left_edge[i])*self.idds[i])
            tmax[i] = (((cur_ind[i]+step[i])*self.dds[i])+
                        self.left_edge[i]-v_pos[i])*iv_dir[i]
            # This deals with the asymmetry in having our indices refer to the
            # left edge of a cell, but the right edge of the brick being one
            # extra zone out.
            if cur_ind[i] == self.dims[i] and step[i] < 0:
                cur_ind[i] = self.dims[i] - 1
            if cur_ind[i] < 0 or cur_ind[i] >= self.dims[i]: return 0
            if step[i] > 0:
                tmax[i] = (((cur_ind[i]+1)*self.dds[i])
                            +self.left_edge[i]-v_pos[i])*iv_dir[i]
            if step[i] < 0:
                tmax[i] = (((cur_ind[i]+0)*self.dds[i])
                            +self.left_edge[i]-v_pos[i])*iv_dir[i]
            tdelta[i] = (self.dds[i]*iv_dir[i])
            if tdelta[i] < 0: tdelta[i] *= -1
        # We have to jumpstart our calculation
        enter_t = intersect_t
        hit = 0
        while 1:
            # dims here is one less than the dimensions of the data,
            # but we are tracing on the grid, not on the data...
            if (not (0 <= cur_ind[0] < self.dims[0])) or \
               (not (0 <= cur_ind[1] < self.dims[1])) or \
               (not (0 <= cur_ind[2] < self.dims[2])):
                break
            hit += 1
            if tmax[0] < tmax[1]:
                if tmax[0] < tmax[2]:
                    exit_t = fmin(tmax[0], 1.0)
                    self.sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind,
                                       rgba, tf)
                    cur_ind[0] += step[0]
                    enter_t = tmax[0]
                    tmax[0] += tdelta[0]
                else:
                    exit_t = fmin(tmax[2], 1.0)
                    self.sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind,
                                       rgba, tf)
                    cur_ind[2] += step[2]
                    enter_t = tmax[2]
                    tmax[2] += tdelta[2]
            else:
                if tmax[1] < tmax[2]:
                    exit_t = fmin(tmax[1], 1.0)
                    self.sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind,
                                       rgba, tf)
                    cur_ind[1] += step[1]
                    enter_t = tmax[1]
                    tmax[1] += tdelta[1]
                else:
                    exit_t = fmin(tmax[2], 1.0)
                    self.sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind,
                                       rgba, tf)
                    cur_ind[2] += step[2]
                    enter_t = tmax[2]
                    tmax[2] += tdelta[2]
            if enter_t >= 1.0: break
        if return_t != NULL: return_t[0] = exit_t
        return hit

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void sample_values(self,
                            np.float64_t v_pos[3],
                            np.float64_t v_dir[3],
                            np.float64_t enter_t,
                            np.float64_t exit_t,
                            int ci[3],
                            np.float64_t *rgba,
                            TransferFunctionProxy tf):
        cdef np.float64_t cp[3], dp[3], pos[3], dt, t, dv
        cdef np.float64_t grad[3], ds[3]
        cdef np.float64_t local_dds[3], cell_left[3]
        grad[0] = grad[1] = grad[2] = 0.0
        cdef int dti, i
        cdef kdtree_utils.kdres *ballq = NULL
        dt = (exit_t - enter_t) / tf.ns # 4 samples should be dt=0.25
        cdef int offset = ci[0] * (self.dims[1] + 1) * (self.dims[2] + 1) \
                        + ci[1] * (self.dims[2] + 1) + ci[2]
        for i in range(3):
            cell_left[i] = ci[i] * self.dds[i] + self.left_edge[i]
            # this gets us dp as the current first sample position
            pos[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
            dp[i] = pos[i] - cell_left[i]
            dp[i] *= self.idds[i]
            ds[i] = v_dir[i] * self.idds[i] * dt
            local_dds[i] = v_dir[i] * dt
        if self.star_list != NULL:
            ballq = kdtree_utils.kd_nearest_range3(
                self.star_list, cell_left[0] + self.dds[0]*0.5,
                                cell_left[1] + self.dds[1]*0.5,
                                cell_left[2] + self.dds[2]*0.5,
                                self.star_er + 0.9*self.dds[0])
                                            # ~0.866 + a bit
        for dti in range(tf.ns): 
            for i in range(self.n_fields):
                self.dvs[i] = offset_interpolate(self.dims, dp, self.data[i] + offset)
            #if (dv < tf.x_bounds[0]) or (dv > tf.x_bounds[1]):
            #    continue
            if self.star_list != NULL: self.add_stars(ballq, dt, pos, rgba)
            tf.eval_transfer(dt, self.dvs, rgba, grad)
            for i in range(3):
                dp[i] += ds[i]
                pos[i] += local_dds[i]
        if ballq != NULL: kdtree_utils.kd_res_free(ballq)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void add_stars(self, kdtree_utils.kdres *ballq,
            np.float64_t dt, np.float64_t pos[3], np.float64_t *rgba):
        cdef int i, n, ns
        cdef double px, py, pz
        cdef np.float64_t gexp, gaussian
        cdef np.float64_t* colors = NULL
        ns = kdtree_utils.kd_res_size(ballq)
        for n in range(ns):
            # We've got Dodgson here!
            kdtree_utils.kd_res_item3(ballq, &px, &py, &pz)
            colors = <np.float64_t *> kdtree_utils.kd_res_item_data(ballq)
            kdtree_utils.kd_res_next(ballq)
            gexp = (px - pos[0])*(px - pos[0]) \
                 + (py - pos[1])*(py - pos[1]) \
                 + (pz - pos[2])*(pz - pos[2])
            gaussian = self.star_coeff * exp(-gexp/self.star_sigma_num)
            for i in range(3): rgba[i] += gaussian*dt*colors[i]
        kdtree_utils.kd_res_rewind(ballq)
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def integrate_streamline(self, pos, np.float64_t h):
        cdef np.float64_t k1[3], k2[3], k3[3], k4[3]
        cdef np.float64_t newpos[3], oldpos[3]
        for i in range(3):
            newpos[i] = oldpos[i] = pos[i]
        self.get_vector_field(newpos, k1)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k1[i]*h

        if not (self.left_edge[0] < newpos[0] and newpos[0] < self.right_edge[0] and \
                self.left_edge[1] < newpos[1] and newpos[1] < self.right_edge[1] and \
                self.left_edge[2] < newpos[2] and newpos[2] < self.right_edge[2]):
            for i in range(3):
                pos[i] = newpos[i]
            return
        
        self.get_vector_field(newpos, k2)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k2[i]*h

        if not (self.left_edge[0] <= newpos[0] and newpos[0] <= self.right_edge[0] and \
                self.left_edge[1] <= newpos[1] and newpos[1] <= self.right_edge[1] and \
                self.left_edge[2] <= newpos[2] and newpos[2] <= self.right_edge[2]):
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k3)
        for i in range(3):
            newpos[i] = oldpos[i] + k3[i]*h
            
        if not (self.left_edge[0] <= newpos[0] and newpos[0] <= self.right_edge[0] and \
                self.left_edge[1] <= newpos[1] and newpos[1] <= self.right_edge[1] and \
                self.left_edge[2] <= newpos[2] and newpos[2] <= self.right_edge[2]):
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k4)

        for i in range(3):
            pos[i] = oldpos[i] + h*(k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_vector_field(self, np.float64_t pos[3],
                               np.float64_t *vel):
        cdef np.float64_t dp[3]
        cdef int ci[3] 
        cdef np.float64_t vel_mag = 0.0
        for i in range(3):
            ci[i] = (int)((pos[i]-self.left_edge[i])/self.dds[i])
            dp[i] = (pos[i] - self.left_edge[i])%(self.dds[i])

        cdef int offset = ci[0] * (self.dims[1] + 1) * (self.dims[2] + 1) \
                          + ci[1] * (self.dims[2] + 1) + ci[2]
        
        for i in range(3):
            vel[i] = offset_interpolate(self.dims, dp, self.data[i] + offset)
            vel_mag += vel[i]*vel[i]
        vel_mag = np.sqrt(vel_mag)
        if vel_mag != 0.0:
            for i in range(3):
                vel[i] /= vel_mag

cdef class GridFace:
    cdef int direction
    cdef public np.float64_t coord
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, grid, int direction, int left):
        self.direction = direction
        if left == 1:
            self.coord = grid.LeftEdge[direction]
        else:
            self.coord = grid.RightEdge[direction]
        cdef int i
        for i in range(3):
            self.left_edge[i] = grid.LeftEdge[i]
            self.right_edge[i] = grid.RightEdge[i]
        self.left_edge[direction] = self.right_edge[direction] = self.coord

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int proj_overlap(self, np.float64_t *left_edge, np.float64_t *right_edge):
        cdef int xax, yax
        xax = (self.direction + 1) % 3
        yax = (self.direction + 2) % 3
        if left_edge[xax] >= self.right_edge[xax]: return 0
        if right_edge[xax] <= self.left_edge[xax]: return 0
        if left_edge[yax] >= self.right_edge[yax]: return 0
        if right_edge[yax] <= self.left_edge[yax]: return 0
        return 1

cdef class ProtoPrism:
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef public object LeftEdge
    cdef public object RightEdge
    cdef public object subgrid_faces
    cdef public int parent_grid_id
    def __cinit__(self, int parent_grid_id,
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  subgrid_faces):
        self.parent_grid_id = parent_grid_id
        cdef int i
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
        self.subgrid_faces = subgrid_faces

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def sweep(self, int direction = 0, int stack = 0):
        cdef int i
        cdef GridFace face
        cdef np.float64_t proto_split[3]
        for i in range(3): proto_split[i] = self.right_edge[i]
        for face in self.subgrid_faces[direction]:
            proto_split[direction] = face.coord
            if proto_split[direction] <= self.left_edge[direction]:
                continue
            if proto_split[direction] == self.right_edge[direction]:
                if stack == 2: return [self]
                return self.sweep((direction + 1) % 3, stack + 1)
            if face.proj_overlap(self.left_edge, proto_split) == 1:
                left, right = self.split(proto_split, direction)
                LC = left.sweep((direction + 1) % 3)
                RC = right.sweep(direction)
                return LC + RC
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef object split(self, np.float64_t *sp, int direction):
        cdef int i
        cdef np.ndarray split_left = self.LeftEdge.copy()
        cdef np.ndarray split_right = self.RightEdge.copy()

        for i in range(3): split_left[i] = self.right_edge[i]
        split_left[direction] = sp[direction]
        left = ProtoPrism(self.parent_grid_id, self.LeftEdge, split_left,
                          self.subgrid_faces)

        for i in range(3): split_right[i] = self.left_edge[i]
        split_right[direction] = sp[direction]
        right = ProtoPrism(self.parent_grid_id, split_right, self.RightEdge,
                           self.subgrid_faces)

        return (left, right)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_brick(self, np.ndarray[np.float64_t, ndim=1] grid_left_edge,
                        np.ndarray[np.float64_t, ndim=1] grid_dds,
                        child_mask):
        # We get passed in the left edge, the dds (which gives dimensions) and
        # the data, which is already vertex-centered.
        cdef PartitionedGrid PG
        cdef int li[3], ri[3], idims[3], i
        for i in range(3):
            li[i] = lrint((self.left_edge[i] - grid_left_edge[i])/grid_dds[i])
            ri[i] = lrint((self.right_edge[i] - grid_left_edge[i])/grid_dds[i])
            idims[i] = ri[i] - li[i]
        if child_mask[li[0], li[1], li[2]] == 0: return []
        cdef np.ndarray[np.int64_t, ndim=1] dims = np.empty(3, dtype='int64')
        for i in range(3):
            dims[i] = idims[i]
        #cdef np.ndarray[np.float64_t, ndim=3] new_data
        #new_data = data[li[0]:ri[0]+1,li[1]:ri[1]+1,li[2]:ri[2]+1].copy()
        #PG = PartitionedGrid(self.parent_grid_id, new_data,
        #                     self.LeftEdge, self.RightEdge, dims)
        return ((li[0], ri[0]), (li[1], ri[1]), (li[2], ri[2]), dims)

cdef struct AdaptiveRayPacket:
    long nside
    long ipix
    np.float64_t t
    np.float64_t v_dir[3]
    np.float64_t value[4]
    np.float64_t pos[3]
    AdaptiveRayPacket *next
    AdaptiveRayPacket *prev
    AdaptiveRayPacket *brick_next
    #int cgi

cdef class AdaptiveRaySource:
    cdef np.float64_t center[3]
    cdef public np.float64_t rays_per_cell
    cdef AdaptiveRayPacket *first
    cdef public np.float64_t normalization
    cdef public int nrays
    cdef public int max_nside
    cdef AdaptiveRayPacket **packet_pointers
    cdef AdaptiveRayPacket **lpacket_pointers

    def __cinit__(self, center, rays_per_cell, initial_nside,
                  np.float64_t normalization, brick_list, int max_nside = 8192):
        cdef int i
        self.max_nside = max_nside
        self.center[0] = center[0]
        self.center[1] = center[1]
        self.center[2] = center[2]
        self.rays_per_cell = rays_per_cell
        cdef AdaptiveRayPacket *ray
        cdef AdaptiveRayPacket *last = NULL
        cdef PartitionedGrid pg
        cdef double v_dir[3]
        cdef int nbricks = len(brick_list)
        # You see, killbots have a preset kill limit. Knowing their weakness, I
        # sent wave after wave of my own men at them until they reached their
        # limit and shut down.
        self.lpacket_pointers = <AdaptiveRayPacket **> \
            malloc(sizeof(AdaptiveRayPacket*)*nbricks)
        self.packet_pointers = <AdaptiveRayPacket **> \
            malloc(sizeof(AdaptiveRayPacket*)*nbricks)
        for i in range(nbricks):
            self.lpacket_pointers[i] = self.packet_pointers[i] = NULL
        self.normalization = normalization
        self.nrays = 12*initial_nside*initial_nside
        for i in range(self.nrays):
            # Initialize rays here
            ray = <AdaptiveRayPacket *> malloc(sizeof(AdaptiveRayPacket))
            ray.prev = last
            ray.ipix = i
            ray.nside = initial_nside
            ray.t = 0.0 # We assume we are not on a brick boundary
            healpix_interface.pix2vec_nest(initial_nside, i, v_dir)
            ray.v_dir[0] = v_dir[0] * normalization
            ray.v_dir[1] = v_dir[1] * normalization
            ray.v_dir[2] = v_dir[2] * normalization
            ray.value[0] = ray.value[1] = ray.value[2] = ray.value[3] = 0.0
            ray.next = NULL
            #ray.cgi = 0
            ray.pos[0] = self.center[0]
            ray.pos[1] = self.center[1]
            ray.pos[2] = self.center[2]
            ray.brick_next = NULL
            if last != NULL:
                last.next = ray
                last.brick_next = ray
            else:
                self.first = ray
            last = ray
        self.packet_pointers[0] = self.first
        self.lpacket_pointers[0] = last

    def __dealloc__(self):
        cdef AdaptiveRayPacket *next
        cdef AdaptiveRayPacket *ray = self.first
        while ray != NULL:
            next = ray.next
            free(ray)
            ray = next
        free(self.packet_pointers)

    def get_rays(self):
        cdef AdaptiveRayPacket *ray = self.first
        cdef int count = 0
        while ray != NULL:
            count += 1
            ray = ray.next
        cdef np.ndarray[np.int64_t, ndim=2] info = np.zeros((count, 2), dtype="int64")
        cdef np.ndarray[np.float64_t, ndim=2] values = np.zeros((count, 4), dtype="float64")
        count = 0
        ray = self.first
        while ray != NULL:
            info[count, 0] = ray.nside
            info[count, 1] = ray.ipix
            values[count, 0] = ray.value[0]
            values[count, 1] = ray.value[1]
            values[count, 2] = ray.value[2]
            values[count, 3] = ray.value[3]
            count += 1
            ray = ray.next
        return info, values

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def integrate_brick(self, PartitionedGrid pg, TransferFunctionProxy tf,
                        int pgi, np.ndarray[np.float64_t, ndim=2] ledges,
                                 np.ndarray[np.float64_t, ndim=2] redges):
        cdef np.float64_t domega
        domega = self.get_domega(pg.left_edge, pg.right_edge)
        #print "dOmega", domega, self.nrays
        cdef int count = 0
        cdef int i, j
        cdef AdaptiveRayPacket *ray = self.packet_pointers[pgi]
        cdef AdaptiveRayPacket *next
        cdef int *grid_neighbors = self.find_neighbors(pgi, pg.dds[0], ledges, redges)
        cdef np.float64_t enter_t, dt, offpos[3]
        cdef int found_a_home, hit
        #print "Grid: ", pgi, "has", grid_neighbors[0], "neighbors"
        while ray != NULL:
            # Note that we may end up splitting a ray such that it ends up
            # outside the brick!  This will likely cause them to get lost.
            #print count
            count +=1
            # We don't need to check for intersection anymore, as we are the
            # Ruler of the planet Omicron Persei 8
            #if self.intersects(ray, pg):
            ray = self.refine_ray(ray, domega, pg.dds[0],
                                  pg.left_edge, pg.right_edge)
            enter_t = ray.t
            hit = pg.integrate_ray(self.center, ray.v_dir, ray.value, tf, &ray.t)
            if hit == 0: dt = 0.0
            else: dt = (ray.t - enter_t)/hit
            for i in range(3):
                ray.pos[i] = ray.v_dir[i] * ray.t + self.center[i]
                offpos[i] = ray.pos[i] + ray.v_dir[i] * 1e-5*dt
            # We set 'next' after the refinement has occurred
            next = ray.brick_next
            found_a_home = 0
            for j in range(grid_neighbors[0]):
                i = grid_neighbors[j+1]
                if ((ledges[i, 0] <= offpos[0] <= redges[i, 0]) and
                    (ledges[i, 1] <= offpos[1] <= redges[i, 1]) and
                    (ledges[i, 2] <= offpos[2] <= redges[i, 2])):
                    if self.lpacket_pointers[i] == NULL:
                        self.packet_pointers[i] = \
                        self.lpacket_pointers[i] = ray
                        ray.brick_next = NULL
                    else:
                        self.lpacket_pointers[i].brick_next = ray
                        self.lpacket_pointers[i] = ray
                        ray.brick_next = NULL
                    #ray.cgi = i
                    found_a_home = 1
                    break
            ray = next
        free(grid_neighbors)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int *find_neighbors(self, int this_grid, np.float64_t dds,
                             np.ndarray[np.float64_t, ndim=2] ledges,
                             np.ndarray[np.float64_t, ndim=2] redges):
        # Iterate once to count the number of grids, iterate a second time to
        # fill the array.  This could be better with a linked list, but I think
        # that the difference should not be substantial.
        cdef int i = 0
        cdef int count = 0
        # We grow our grids by dx in every direction, then look for overlap
        cdef np.float64_t gle[3], gre[3]
        gle[0] = ledges[this_grid, 0] - dds
        gle[1] = ledges[this_grid, 1] - dds
        gle[2] = ledges[this_grid, 2] - dds
        gre[0] = redges[this_grid, 0] + dds
        gre[1] = redges[this_grid, 1] + dds
        gre[2] = redges[this_grid, 2] + dds
        for i in range(this_grid+1, ledges.shape[0]):
            # Check for overlap
            if ((gle[0] <= redges[i, 0] and gre[0] >= ledges[i, 0]) and
                (gle[1] <= redges[i, 1] and gre[1] >= ledges[i, 1]) and
                (gle[2] <= redges[i, 2] and gre[2] >= ledges[i, 2])):
                count += 1
        cdef int *tr = <int *> malloc(sizeof(int) * (count + 1))
        tr[0] = count
        count = 0
        for i in range(this_grid+1, ledges.shape[0]):
            # Check for overlap
            if ((gle[0] <= redges[i, 0] and gre[0] >= ledges[i, 0]) and
                (gle[1] <= redges[i, 1] and gre[1] >= ledges[i, 1]) and
                (gle[2] <= redges[i, 2] and gre[2] >= ledges[i, 2])):
                tr[count + 1] = i
                count += 1
        return tr

    cdef int intersects(self, AdaptiveRayPacket *ray, PartitionedGrid pg):
        cdef np.float64_t pos[3]
        cdef int i
        for i in range(3):
            # Is this correct, for the normalized v_dir?
            if ray.pos[i] < pg.left_edge[i]: return 0
            if ray.pos[i] > pg.right_edge[i]: return 0
        return 1

    cdef np.float64_t get_domega(self, np.float64_t left_edge[3],
                                       np.float64_t right_edge[3]):
        # We should calculate the subtending angle at the maximum radius of the
        # brick
        cdef int i, j, k
        cdef np.float64_t r2[3], max_r2, domega, *edge[2]
        # We now figure out when the ray will leave the box, at worst.
        edge[0] = left_edge
        edge[1] = right_edge
        max_r2 = -1.0
        for i in range(2):
            r2[0] = (edge[i][0] - self.center[0])**2.0
            for j in range(2):
                r2[1] = r2[0] + (edge[j][1] - self.center[1])**2.0
                for k in range(3):
                    r2[2] = r2[1] + (edge[k][2] - self.center[2])**2.0
                    max_r2 = fmax(max_r2, r2[2])
        domega = 4.0 * 3.1415926 * max_r2 # Used to be / Nrays
        return domega

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef AdaptiveRayPacket *refine_ray(self, AdaptiveRayPacket *ray,
                                       np.float64_t domega,
                                       np.float64_t dx,
                                       np.float64_t left_edge[3],
                                       np.float64_t right_edge[3]):
        # This should be recursive; we are not correctly applying split
        # criteria multiple times.
        cdef long Nrays = 12 * ray.nside * ray.nside
        cdef int i, j
        if domega/Nrays < dx*dx/self.rays_per_cell:
            #print ray.nside, domega/Nrays, dx, (domega/Nrays * self.rays_per_cell)/(dx*dx)
            return ray
        if ray.nside >= self.max_nside: return ray
        #print "Refining %s from %s to %s" % (ray.ipix, ray.nside, ray.nside*2)
        # Now we make four new ones
        cdef double v_dir[3]
        cdef AdaptiveRayPacket *prev = ray.prev
        # It is important to note here that brick_prev is a local variable for
        # the newly created rays, not the previous ray in this brick, as that
        # has already been passed on to its next brick
        cdef AdaptiveRayPacket *brick_prev = NULL
        for i in range(4):
            new_ray = <AdaptiveRayPacket *> malloc(
                            sizeof(AdaptiveRayPacket))
            new_ray.nside = ray.nside * 2
            new_ray.ipix = ray.ipix * 4 + i
            new_ray.t = ray.t
            #new_ray.cgi = ray.cgi
            new_ray.prev = prev
            if new_ray.prev != NULL:
                new_ray.prev.next = new_ray
            if brick_prev != NULL:
                brick_prev.brick_next = new_ray
            prev = brick_prev = new_ray
            healpix_interface.pix2vec_nest(
                    new_ray.nside, new_ray.ipix, v_dir)
            for j in range(3):
                new_ray.v_dir[j] = v_dir[j] * self.normalization
                new_ray.value[j] = ray.value[j]
                new_ray.pos[j] = self.center[j] + ray.t * new_ray.v_dir[j]
            new_ray.value[3] = ray.value[3]

        new_ray.next = ray.next
        new_ray.brick_next = ray.brick_next
        if new_ray.next != NULL:
            new_ray.next.prev = new_ray
        if self.first == ray:
            self.first = new_ray.prev.prev.prev
        free(ray)
        self.nrays += 3
        return new_ray.prev.prev.prev

# From Enzo:
#   dOmega = 4 pi r^2/Nrays
#   if (dOmega > RaysPerCell * dx^2) then split
