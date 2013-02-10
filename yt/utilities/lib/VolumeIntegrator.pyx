"""
Simple integrators for the radiative transfer equation

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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
from libc.stdlib cimport malloc, free, abs

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
    long double expl(long double x)
    double floor(double x)
    double ceil(double x)
    double fmod(double x, double y)
    double log2(double x)
    long int lrint(double x)
    double fabs(double x)
    double cos(double x)
    double sin(double x)
    double asin(double x)
    double acos(double x)

cdef struct Triangle:
    Triangle *next
    np.float64_t p[3][3]
    np.float64_t val

cdef struct TriangleCollection:
    int count
    Triangle *first
    Triangle *current

cdef Triangle *AddTriangle(Triangle *self,
                    np.float64_t p0[3], np.float64_t p1[3], np.float64_t p2[3]):
    cdef Triangle *nn = <Triangle *> malloc(sizeof(Triangle))
    if self != NULL:
        self.next = nn
    cdef int i
    for i in range(3):
        nn.p[0][i] = p0[i]
    for i in range(3):
        nn.p[1][i] = p1[i]
    for i in range(3):
        nn.p[2][i] = p2[i]
    nn.next = NULL
    return nn

cdef int CountTriangles(Triangle *first):
    cdef int count = 0
    cdef Triangle *this = first
    while this != NULL:
        count += 1
        this = this.next
    return count

cdef void FillTriangleValues(np.ndarray[np.float64_t, ndim=1] values,
                             Triangle *first):
    cdef Triangle *this = first
    cdef Triangle *last
    cdef int i = 0
    while this != NULL:
        values[i] = this.val
        i += 1
        last = this
        this = this.next

cdef void WipeTriangles(Triangle *first):
    cdef Triangle *this = first
    cdef Triangle *last
    while this != NULL:
        last = this
        this = this.next
        free(last)

cdef void FillAndWipeTriangles(np.ndarray[np.float64_t, ndim=2] vertices,
                               Triangle *first):
    cdef int count = 0
    cdef Triangle *this = first
    cdef Triangle *last
    cdef int i, j
    while this != NULL:
        for i in range(3):
            for j in range(3):
                vertices[count, j] = this.p[i][j]
            count += 1 # Do it at the end because it's an index
        last = this
        this = this.next
        free(last)

cdef extern from "FixedInterpolator.h":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data)
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3], np.float64_t *data)
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data)
    void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                       np.float64_t grad[3])
    void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval)
    void vertex_interp(np.float64_t v1, np.float64_t v2, np.float64_t isovalue,
                       np.float64_t vl[3], np.float64_t dds[3],
                       np.float64_t x, np.float64_t y, np.float64_t z,
                       int vind1, int vind2)

#cdef extern int *edge_table
#cdef extern int **tri_table

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

def arr_fisheye_vectors(int resolution, np.float64_t fov, int nimx=1, int
        nimy=1, int nimi=0, int nimj=0, np.float64_t off_theta=0.0, np.float64_t
        off_phi=0.0):
    # We now follow figures 4-7 of:
    # http://paulbourke.net/miscellaneous/domefisheye/fisheye/
    # ...but all in Cython.
    cdef np.ndarray[np.float64_t, ndim=3] vp
    cdef int i, j, k
    cdef np.float64_t r, phi, theta, px, py
    cdef np.float64_t pi = 3.1415926
    cdef np.float64_t fov_rad = fov * pi / 180.0
    cdef int nx = resolution/nimx
    cdef int ny = resolution/nimy
    vp = np.zeros((nx,ny, 3), dtype="float64")
    for i in range(nx):
        px = 2.0 * (nimi*nx + i) / (resolution) - 1.0
        for j in range(ny):
            py = 2.0 * (nimj*ny + j) / (resolution) - 1.0
            r = (px*px + py*py)**0.5
            if r == 0.0:
                phi = 0.0
            elif px < 0:
                phi = pi - asin(py / r)
            else:
                phi = asin(py / r)
            theta = r * fov_rad / 2.0
            theta += off_theta
            phi += off_phi
            vp[i,j,0] = sin(theta) * cos(phi)
            vp[i,j,1] = sin(theta) * sin(phi)
            vp[i,j,2] = cos(theta)
    return vp

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
    cdef np.float64_t bv, dy, dd, tf, rv
    cdef int bin_id
    if fit.pass_through == 1:
        rv = dvs[fit.field_id]
        if fit.weight_field_id != -1: rv *= dvs[fit.weight_field_id]
        return rv
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
            #print "Field table", i, "corresponds to",
            #print self.field_tables[i].field_id,
            #print "(Weighted with ", self.field_tables[i].weight_field_id,
            #print ")"

        for i in range(6):
            self.field_table_ids[i] = tf_obj.field_table_ids[i]
            #print "Channel", i, "corresponds to", self.field_table_ids[i]
            
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void eval_transfer(self, np.float64_t dt, np.float64_t *dvs,
                                  np.float64_t *rgba, np.float64_t *grad):
        cdef int i, fid, use
        cdef np.float64_t ta, tf, istorage[6], trgba[6], dot_prod
        # NOTE: We now disable this.  I have left it to ease the process of
        # potentially, one day, re-including it.
        #use = 0
        #for i in range(self.n_field_tables):
        #    fid = self.field_tables[i].field_id
        #    if (dvs[fid] >= self.field_tables[i].bounds[0]) and \
        #       (dvs[fid] <= self.field_tables[i].bounds[1]):
        #        use = 1
        #        break
        for i in range(6): istorage[i] = 0.0
        for i in range(self.n_field_tables):
            istorage[i] = FIT_get_value(&self.field_tables[i], dvs)
        # We have to do this after the interpolation
        for i in range(self.n_field_tables):
            fid = self.field_tables[i].weight_table_id
            if fid != -1: istorage[i] *= istorage[fid]
        for i in range(6):
            trgba[i] = istorage[self.field_table_ids[i]]
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
                                 np.float64_t *return_t = NULL,
                                 np.float64_t enter_t = -1.0):
        cdef int cur_ind[3], step[3], x, y, i, n, flat_ind, hit, direction
        cdef np.float64_t intersect_t = 1.0
        cdef np.float64_t iv_dir[3]
        cdef np.float64_t intersect[3], tmax[3], tdelta[3]
        cdef np.float64_t dist, alpha, dt, exit_t
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
        if enter_t >= 0.0: intersect_t = enter_t
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
        # The initial and final values can be linearly interpolated between; so
        # we just have to calculate our initial and final values.
        cdef np.float64_t slopes[6]
        for i in range(3):
            dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
            dp[i] -= ci[i] * self.dds[i] + self.left_edge[i]
            dp[i] *= self.idds[i]
            ds[i] = v_dir[i] * self.idds[i] * dt
        for i in range(self.n_fields):
            slopes[i] = offset_interpolate(self.dims, dp,
                            self.data[i] + offset)
        for i in range(3):
            dp[i] += ds[i] * tf.ns
        cdef np.float64_t temp
        for i in range(self.n_fields):
            temp = slopes[i]
            slopes[i] -= offset_interpolate(self.dims, dp,
                             self.data[i] + offset)
            slopes[i] *= -1.0/tf.ns
            self.dvs[i] = temp
        if self.star_list != NULL:
            for i in range(3):
                cell_left[i] = ci[i] * self.dds[i] + self.left_edge[i]
                # this gets us dp as the current first sample position
                pos[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
                dp[i] -= tf.ns * ds[i]
                local_dds[i] = v_dir[i] * dt
            ballq = kdtree_utils.kd_nearest_range3(
                self.star_list, cell_left[0] + self.dds[0]*0.5,
                                cell_left[1] + self.dds[1]*0.5,
                                cell_left[2] + self.dds[2]*0.5,
                                self.star_er + 0.9*self.dds[0])
                                            # ~0.866 + a bit
        for dti in range(tf.ns): 
            #if (dv < tf.x_bounds[0]) or (dv > tf.x_bounds[1]):
            #    continue
            if self.star_list != NULL:
                self.add_stars(ballq, dt, pos, rgba)
                for i in range(3):
                    dp[i] += ds[i]
                    pos[i] += local_dds[i]
            tf.eval_transfer(dt, self.dvs, rgba, grad)
            for i in range(self.n_fields):
                self.dvs[i] += slopes[i]
        if ballq != NULL: kdtree_utils.kd_res_free(ballq)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void add_stars(self, kdtree_utils.kdres *ballq,
            np.float64_t dt, np.float64_t pos[3], np.float64_t *rgba):
        cdef int i, n, ns
        cdef np.float64_t px, py, pz
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
            gaussian = self.star_coeff * expl(-gexp/self.star_sigma_num)
            for i in range(3): rgba[i] += gaussian*dt*colors[i]
        kdtree_utils.kd_res_rewind(ballq)
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def integrate_streamline(self, pos, np.float64_t h, mag):
        cdef np.float64_t cmag[1]
        cdef np.float64_t k1[3], k2[3], k3[3], k4[3]
        cdef np.float64_t newpos[3], oldpos[3]
        for i in range(3):
            newpos[i] = oldpos[i] = pos[i]
        self.get_vector_field(newpos, k1, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k1[i]*h

        if not (self.left_edge[0] < newpos[0] and newpos[0] < self.right_edge[0] and \
                self.left_edge[1] < newpos[1] and newpos[1] < self.right_edge[1] and \
                self.left_edge[2] < newpos[2] and newpos[2] < self.right_edge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return
        
        self.get_vector_field(newpos, k2, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + 0.5*k2[i]*h

        if not (self.left_edge[0] <= newpos[0] and newpos[0] <= self.right_edge[0] and \
                self.left_edge[1] <= newpos[1] and newpos[1] <= self.right_edge[1] and \
                self.left_edge[2] <= newpos[2] and newpos[2] <= self.right_edge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k3, cmag)
        for i in range(3):
            newpos[i] = oldpos[i] + k3[i]*h
            
        if not (self.left_edge[0] <= newpos[0] and newpos[0] <= self.right_edge[0] and \
                self.left_edge[1] <= newpos[1] and newpos[1] <= self.right_edge[1] and \
                self.left_edge[2] <= newpos[2] and newpos[2] <= self.right_edge[2]):
            if mag is not None:
                mag[0] = cmag[0]
            for i in range(3):
                pos[i] = newpos[i]
            return

        self.get_vector_field(newpos, k4, cmag)

        for i in range(3):
            pos[i] = oldpos[i] + h*(k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0)

        if mag is not None:
            for i in range(3):
                newpos[i] = pos[i]
            self.get_vector_field(newpos, k4, cmag)
            mag[0] = cmag[0]
            
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void get_vector_field(self, np.float64_t pos[3],
                               np.float64_t *vel, np.float64_t *vel_mag):
        cdef np.float64_t dp[3]
        cdef int ci[3] 
        for i in range(3):
            ci[i] = (int)((pos[i]-self.left_edge[i])/self.dds[i])
            dp[i] = (pos[i] - self.left_edge[i])%(self.dds[i])

        cdef int offset = ci[0] * (self.dims[1] + 1) * (self.dims[2] + 1) \
                          + ci[1] * (self.dims[2] + 1) + ci[2]
        
        vel_mag[0] = 0.0
        for i in range(3):
            vel[i] = offset_interpolate(self.dims, dp, self.data[i] + offset)
            vel_mag[0] += vel[i]*vel[i]
        vel_mag[0] = np.sqrt(vel_mag[0])
        if vel_mag[0] != 0.0:
            for i in range(3):
                vel[i] /= vel_mag[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_isocontour_triangles(self, np.float64_t isovalue, int field_id = 0):
        # Much of this was inspired by code from Paul Bourke's website:
        # http://paulbourke.net/geometry/polygonise/
        # Cython makes us toss this in here, which I think will change in a
        # future release.

        cdef int i, j, k, n
        cdef int offset
        cdef np.float64_t gv[8]
        cdef int cubeindex
        cdef np.float64_t *intdata = NULL
        cdef np.float64_t x, y, z
        cdef np.float64_t mu
        cdef TriangleCollection triangles
        triangles.first = triangles.current = NULL
        triangles.count = 0
        x = self.left_edge[0]
        for i in range(self.dims[0]):
            y = self.left_edge[1]
            for j in range(self.dims[1]):
                z = self.left_edge[2]
                for k in range(self.dims[2]):
                    offset = i * (self.dims[1] + 1) * (self.dims[2] + 1) \
                           + j * (self.dims[2] + 1) + k
                    intdata = self.data[field_id] + offset
                    offset_fill(self.dims, intdata, gv)
                    march_cubes(gv, isovalue, self.dds, x, y, z,
                                &triangles)
                    z += self.dds[2]
                y += self.dds[1]
            x += self.dds[0]
        # Hallo, we are all done.
        cdef np.ndarray[np.float64_t, ndim=2] vertices 
        vertices = np.zeros((triangles.count*3,3), dtype='float64')
        FillAndWipeTriangles(vertices, triangles.first)
        return vertices

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def march_cubes_grid(np.float64_t isovalue,
                     np.ndarray[np.float64_t, ndim=3] values,
                     np.ndarray[np.int32_t, ndim=3] mask,
                     np.ndarray[np.float64_t, ndim=1] left_edge,
                     np.ndarray[np.float64_t, ndim=1] dxs,
                     obj_sample = None):
    cdef int dims[3]
    cdef int i, j, k, n, m, nt
    cdef int offset
    cdef np.float64_t gv[8], pos[3], point[3], idds[3]
    cdef np.float64_t *intdata = NULL
    cdef np.float64_t *sdata = NULL
    cdef np.float64_t x, y, z, do_sample
    cdef np.ndarray[np.float64_t, ndim=3] sample
    cdef np.ndarray[np.float64_t, ndim=1] sampled
    cdef TriangleCollection triangles
    cdef Triangle *last, *current
    if obj_sample is not None:
        sample = obj_sample
        sdata = <np.float64_t *> sample.data
        do_sample = 1
    else:
        do_sample = 0
    for i in range(3):
        dims[i] = values.shape[i] - 1
        idds[i] = 1.0 / dxs[i]
    triangles.first = triangles.current = NULL
    last = current = NULL
    triangles.count = 0
    cdef np.float64_t *data = <np.float64_t *> values.data
    cdef np.float64_t *dds = <np.float64_t *> dxs.data
    pos[0] = left_edge[0]
    for i in range(dims[0]):
        pos[1] = left_edge[1]
        for j in range(dims[1]):
            pos[2] = left_edge[2]
            for k in range(dims[2]):
                if mask[i,j,k] == 1:
                    offset = i * (dims[1] + 1) * (dims[2] + 1) \
                           + j * (dims[2] + 1) + k
                    intdata = data + offset
                    offset_fill(dims, intdata, gv)
                    nt = march_cubes(gv, isovalue, dds, pos[0], pos[1], pos[2],
                                &triangles)
                    if do_sample == 1 and nt > 0:
                        # At each triangle's center, sample our secondary field
                        if last == NULL and triangles.first != NULL:
                            current = triangles.first
                            last = NULL
                        elif last != NULL:
                            current = last.next
                        while current != NULL:
                            for n in range(3):
                                point[n] = 0.0
                            for n in range(3):
                                for m in range(3):
                                    point[m] += (current.p[n][m]-pos[m])*idds[m]
                            for n in range(3):
                                point[n] /= 3.0
                            current.val = offset_interpolate(dims, point,
                                                             sdata + offset)
                            last = current
                            if current.next == NULL: break
                            current = current.next
                pos[2] += dds[2]
            pos[1] += dds[1]
        pos[0] += dds[0]
    # Hallo, we are all done.
    cdef np.ndarray[np.float64_t, ndim=2] vertices 
    vertices = np.zeros((triangles.count*3,3), dtype='float64')
    if do_sample == 1:
        sampled = np.zeros(triangles.count, dtype='float64')
        FillTriangleValues(sampled, triangles.first)
        FillAndWipeTriangles(vertices, triangles.first)
        return vertices, sampled
    FillAndWipeTriangles(vertices, triangles.first)
    return vertices

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def march_cubes_grid_flux(
                     np.float64_t isovalue,
                     np.ndarray[np.float64_t, ndim=3] values,
                     np.ndarray[np.float64_t, ndim=3] v1,
                     np.ndarray[np.float64_t, ndim=3] v2,
                     np.ndarray[np.float64_t, ndim=3] v3,
                     np.ndarray[np.float64_t, ndim=3] flux_field,
                     np.ndarray[np.int32_t, ndim=3] mask,
                     np.ndarray[np.float64_t, ndim=1] left_edge,
                     np.ndarray[np.float64_t, ndim=1] dxs):
    cdef int dims[3]
    cdef int i, j, k, n, m
    cdef int offset
    cdef np.float64_t gv[8]
    cdef np.float64_t *intdata = NULL
    cdef TriangleCollection triangles
    cdef Triangle *current = NULL
    cdef Triangle *last = NULL
    cdef np.float64_t *data = <np.float64_t *> values.data
    cdef np.float64_t *v1data = <np.float64_t *> v1.data
    cdef np.float64_t *v2data = <np.float64_t *> v2.data
    cdef np.float64_t *v3data = <np.float64_t *> v3.data
    cdef np.float64_t *fdata = <np.float64_t *> flux_field.data
    cdef np.float64_t *dds = <np.float64_t *> dxs.data
    cdef np.float64_t flux = 0.0
    cdef np.float64_t center[3], point[3], wval, temp, area, s
    cdef np.float64_t cell_pos[3], fv[3], idds[3], normal[3]
    for i in range(3):
        dims[i] = values.shape[i] - 1
        idds[i] = 1.0 / dds[i]
    triangles.first = triangles.current = NULL
    triangles.count = 0
    cell_pos[0] = left_edge[0]
    for i in range(dims[0]):
        cell_pos[1] = left_edge[1]
        for j in range(dims[1]):
            cell_pos[2] = left_edge[2]
            for k in range(dims[2]):
                if mask[i,j,k] == 1:
                    offset = i * (dims[1] + 1) * (dims[2] + 1) \
                           + j * (dims[2] + 1) + k
                    intdata = data + offset
                    offset_fill(dims, intdata, gv)
                    march_cubes(gv, isovalue, dds,
                                cell_pos[0], cell_pos[1], cell_pos[2],
                                &triangles)
                    # Now our triangles collection has a bunch.  We now
                    # calculate fluxes for each.
                    if last == NULL and triangles.first != NULL:
                        current = triangles.first
                        last = NULL
                    elif last != NULL:
                        current = last.next
                    while current != NULL:
                        # Calculate the center of the triangle
                        wval = 0.0
                        for n in range(3):
                            center[n] = 0.0
                        for n in range(3):
                            for m in range(3):
                                point[m] = (current.p[n][m]-cell_pos[m])*idds[m]
                            # Now we calculate the value at this point
                            temp = offset_interpolate(dims, point, intdata)
                            #print "something", temp, point[0], point[1], point[2]
                            wval += temp
                            for m in range(3):
                                center[m] += temp * point[m]
                        # Now we divide by our normalizing factor
                        for n in range(3):
                            center[n] /= wval
                        # We have our center point of the triangle, in 0..1
                        # coordinates.  So now we interpolate our three
                        # fields.
                        fv[0] = offset_interpolate(dims, center, v1data + offset)
                        fv[1] = offset_interpolate(dims, center, v2data + offset)
                        fv[2] = offset_interpolate(dims, center, v3data + offset)
                        # We interpolate again the actual value data
                        wval = offset_interpolate(dims, center, fdata + offset)
                        # Now we have our flux vector and our field value!
                        # We just need a normal vector with which we can
                        # dot it.  The normal should be equal to the gradient
                        # in the center of the triangle, or thereabouts.
                        eval_gradient(dims, center, intdata, normal)
                        temp = 0.0
                        for n in range(3):
                            temp += normal[n]*normal[n]
                        # Take the negative, to ensure it points inwardly
                        temp = -(temp**0.5)
                        # Dump this somewhere for now
                        temp = wval * (fv[0] * normal[0] +
                                       fv[1] * normal[1] +
                                       fv[2] * normal[2])/temp
                        # Now we need the area of the triangle.  This will take
                        # a lot of time to calculate compared to the rest.
                        # We use Heron's formula.
                        for n in range(3):
                            fv[n] = 0.0
                        for n in range(3):
                            fv[0] += (current.p[0][n] - current.p[2][n])**2.0
                            fv[1] += (current.p[1][n] - current.p[0][n])**2.0
                            fv[2] += (current.p[2][n] - current.p[1][n])**2.0
                        s = 0.0
                        for n in range(3):
                            fv[n] = fv[n]**0.5
                            s += 0.5 * fv[n]
                        area = (s*(s-fv[0])*(s-fv[1])*(s-fv[2]))
                        area = area**0.5
                        flux += temp*area
                        last = current
                        if current.next == NULL: break
                        current = current.next
                cell_pos[2] += dds[2]
            cell_pos[1] += dds[1]
        cell_pos[0] += dds[0]
    # Hallo, we are all done.
    WipeTriangles(triangles.first)
    return flux

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int march_cubes(
                 np.float64_t gv[8], np.float64_t isovalue,
                 np.float64_t dds[3],
                 np.float64_t x, np.float64_t y, np.float64_t z,
                 TriangleCollection *triangles):
    cdef int *edge_table=[
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   ]

    cdef int **tri_table = \
    [[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
    [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
    [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
    [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
    [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
    [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
    [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
    [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
    [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
    [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
    [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
    [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
    [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
    [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
    [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
    [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
    [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
    [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
    [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
    [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
    [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
    [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
    [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
    [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
    [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
    [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
    [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
    [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
    [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
    [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
    [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
    [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
    [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
    [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
    [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
    [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
    [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
    [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
    [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
    [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
    [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
    [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
    [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
    [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
    [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
    [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
    [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
    [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
    [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
    [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
    [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
    [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
    [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
    [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
    [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
    [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
    [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
    [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
    [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
    [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
    [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
    [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
    [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
    [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
    [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
    [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
    [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
    [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
    [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
    [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
    [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
    [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
    [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
    [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
    [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
    [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
    [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
    [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
    [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
    [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
    [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
    [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
    [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
    [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
    [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
    [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
    [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
    [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
    [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
    [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
    [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
    [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
    [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
    [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
    [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
    [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
    [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
    [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
    [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
    [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
    [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
    [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
    [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
    [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
    [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
    [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
    [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
    [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
    [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
    [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
    [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
    [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
    [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
    [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
    [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
    [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
    [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
    [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
    [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
    [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
    [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
    [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
    [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
    [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
    [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
    [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
    [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
    [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
    [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
    [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
    [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
    [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
    [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
    [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
    [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
    [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
    [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
    [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
    [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
    [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
    [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
    [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
    [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
    [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
    [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
    [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
    [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
    [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
    [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
    [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
    [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
    [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
    [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
    [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
    [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
    [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]]
    cdef np.float64_t vertlist[12][3]
    cdef int cubeindex = 0
    cdef int n
    cdef int nt = 0
    for n in range(8):
        if gv[n] < isovalue:
            cubeindex |= (1 << n)
    if edge_table[cubeindex] == 0:
        return 0
    if (edge_table[cubeindex] & 1): # 0,0,0 with 1,0,0
        vertex_interp(gv[0], gv[1], isovalue, vertlist[0],
                      dds, x, y, z, 0, 1)
    if (edge_table[cubeindex] & 2): # 1,0,0 with 1,1,0
        vertex_interp(gv[1], gv[2], isovalue, vertlist[1],
                      dds, x, y, z, 1, 2)
    if (edge_table[cubeindex] & 4): # 1,1,0 with 0,1,0
        vertex_interp(gv[2], gv[3], isovalue, vertlist[2],
                      dds, x, y, z, 2, 3)
    if (edge_table[cubeindex] & 8): # 0,1,0 with 0,0,0
        vertex_interp(gv[3], gv[0], isovalue, vertlist[3],
                      dds, x, y, z, 3, 0)
    if (edge_table[cubeindex] & 16): # 0,0,1 with 1,0,1
        vertex_interp(gv[4], gv[5], isovalue, vertlist[4],
                      dds, x, y, z, 4, 5)
    if (edge_table[cubeindex] & 32): # 1,0,1 with 1,1,1
        vertex_interp(gv[5], gv[6], isovalue, vertlist[5],
                      dds, x, y, z, 5, 6)
    if (edge_table[cubeindex] & 64): # 1,1,1 with 0,1,1
        vertex_interp(gv[6], gv[7], isovalue, vertlist[6],
                      dds, x, y, z, 6, 7)
    if (edge_table[cubeindex] & 128): # 0,1,1 with 0,0,1
        vertex_interp(gv[7], gv[4], isovalue, vertlist[7],
                      dds, x, y, z, 7, 4)
    if (edge_table[cubeindex] & 256): # 0,0,0 with 0,0,1
        vertex_interp(gv[0], gv[4], isovalue, vertlist[8],
                      dds, x, y, z, 0, 4)
    if (edge_table[cubeindex] & 512): # 1,0,0 with 1,0,1
        vertex_interp(gv[1], gv[5], isovalue, vertlist[9],
                      dds, x, y, z, 1, 5)
    if (edge_table[cubeindex] & 1024): # 1,1,0 with 1,1,1
        vertex_interp(gv[2], gv[6], isovalue, vertlist[10],
                      dds, x, y, z, 2, 6)
    if (edge_table[cubeindex] & 2048): # 0,1,0 with 0,1,1
        vertex_interp(gv[3], gv[7], isovalue, vertlist[11],
                      dds, x, y, z, 3, 7)
    n = 0
    while 1:
        triangles.current = AddTriangle(triangles.current,
                    vertlist[tri_table[cubeindex][n  ]],
                    vertlist[tri_table[cubeindex][n+1]],
                    vertlist[tri_table[cubeindex][n+2]])
        triangles.count += 1
        nt += 1
        if triangles.first == NULL:
            triangles.first = triangles.current
        n += 3
        if tri_table[cubeindex][n] == -1: break
    return nt
    

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
                        np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask):
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
    int pgi

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
                  np.float64_t normalization, brick_list, 
                  np.ndarray[np.float64_t, ndim=2] ledges,
                  np.ndarray[np.float64_t, ndim=2] redges,
                  int max_nside = 8192):
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
        cdef int *grid_neighbors = <int*> malloc(sizeof(int) * (nbricks+1))
        grid_neighbors[0] = nbricks
        for i in range(nbricks):
            grid_neighbors[i+1] = i
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
            ray.pgi = -1
            if last != NULL:
                last.next = ray
            else:
                self.first = ray
            self.send_ray_home(ray, ledges, redges, grid_neighbors, 0)
            last = ray

    def __dealloc__(self):
        cdef AdaptiveRayPacket *next
        cdef AdaptiveRayPacket *ray = self.first
        while ray != NULL:
            next = ray.next
            free(ray)
            ray = next
        free(self.packet_pointers)
        free(self.lpacket_pointers)

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
            if ray.t < 0.5:
                print "PROBLEM",
                print count, ray.ipix, ray.nside, ray.t,
                print "vd", ray.v_dir[0], ray.v_dir[1], ray.v_dir[2],
                print "pos", ray.pos[0], ray.pos[1], ray.pos[2],
                print "pgi", ray.pgi
            count += 1
            ray = ray.next
        return info, values


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef np.float64_t integrate_ray(self, AdaptiveRayPacket *ray,
                      PartitionedGrid pg, TransferFunctionProxy tf):
        cdef np.float64_t self_center[3], ray_v_dir[3], ray_value[4]
        self_center[0] = self.center[0]
        self_center[1] = self.center[1]
        self_center[2] = self.center[2]
        enter_t = ray.t
        ray_v_dir[0] = ray.v_dir[0]
        ray_v_dir[1] = ray.v_dir[1]
        ray_v_dir[2] = ray.v_dir[2]
        ray_value[0] = ray.value[0]
        ray_value[1] = ray.value[1]
        ray_value[2] = ray.value[2]
        ray_value[3] = ray.value[3]
        hit = pg.integrate_ray(self_center, ray_v_dir, ray_value, tf, &ray.t,
                               ray.t)
        ray.value[0] = ray_value[0]
        ray.value[1] = ray_value[1]
        ray.value[2] = ray_value[2]
        ray.value[3] = ray_value[3]
        if hit == 0: dt = 0.0
        else: dt = (ray.t - enter_t)/hit
        for i in range(3):
            ray.pos[i] = ray.v_dir[i] * ray.t + self.center[i]
        return dt

    cdef send_ray_home(self, AdaptiveRayPacket *ray,
                       np.ndarray[np.float64_t, ndim=2] ledges,
                       np.ndarray[np.float64_t, ndim=2] redges,
                       int *grid_neighbors, np.float64_t dt,
                       int skip_append = 0):
        cdef int found_a_home = 0
        cdef int i, j, npgi
        cdef np.float64_t offpos[3]
        for i in range(3):
            offpos[i] = ray.pos[i] + ray.v_dir[i] * 1e-6*dt
        for j in range(grid_neighbors[0]):
            i = grid_neighbors[j+1]
            if ((ledges[i, 0] <= offpos[0] <= redges[i, 0]) and
                (ledges[i, 1] <= offpos[1] <= redges[i, 1]) and
                (ledges[i, 2] <= offpos[2] <= redges[i, 2]) and
                ray.pgi != i):
                if not skip_append: self.append_to_packets(i, ray)
                ray.pgi = i
                npgi = i
                found_a_home = 1
                break
        if found_a_home == 0:
            #print "Non-neighboring area", ray.pgi, ray.ipix, ray.nside
            for i in range(ledges.shape[0]):
                if ((ledges[i, 0] <= offpos[0] <= redges[i, 0]) and
                    (ledges[i, 1] <= offpos[1] <= redges[i, 1]) and
                    (ledges[i, 2] <= offpos[2] <= redges[i, 2])):
                    #print "Found a home!", i, ray.ipix, ray.nside, ray.pgi
                    if not skip_append: self.append_to_packets(i, ray)
                    ray.pgi = i
                    npgi = i
                    found_a_home = 1
                    break
            if found_a_home == 0:
                raise RuntimeError

    cdef append_to_packets(self, int pgi, AdaptiveRayPacket *ray):
        # packet_pointers are pointers to the *first* packet in a given brick
        # lpacket_pointers point to the *final* packet in a given brick, for
        # easy appending.
        if self.lpacket_pointers[pgi] == NULL or \
           self.packet_pointers[pgi] == NULL:
            self.packet_pointers[pgi] = \
            self.lpacket_pointers[pgi] = ray
            ray.brick_next = NULL
        else:
            self.lpacket_pointers[pgi].brick_next = ray
            self.lpacket_pointers[pgi] = ray
            ray.brick_next = NULL

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def integrate_brick(self, PartitionedGrid pg, TransferFunctionProxy tf,
                        int pgi, np.ndarray[np.float64_t, ndim=2] ledges,
                                 np.ndarray[np.float64_t, ndim=2] redges,
                                 pgs, int inside = -1):
        cdef np.float64_t domega, dt
        cdef PartitionedGrid pgn
        #print "dOmega", domega, self.nrays
        cdef intersects
        cdef int i, j, npgi, refined
        cdef AdaptiveRayPacket *ray2, *ray = self.packet_pointers[pgi]
        cdef AdaptiveRayPacket *next
        cdef AdaptiveRayPacket **pray
        cdef int *grid_neighbors = self.find_neighbors(pgi, pg.dds[0], ledges, redges)
        cdef int *grid2_neighbors
        cdef np.float64_t enter_t, offpos[3]
        cdef int found_a_home, hit
        while ray != NULL:
            #print "Integrating", pgi, ray.pgi, ray.ipix, ray.nside
            if pgi != ray.pgi:
                self.send_ray_home(ray, ledges, redges, grid_neighbors, dt)
                if ray.pgi != pgi:
                    ray = ray.brick_next
                    continue
            # We start in this brick, and then we integrate to the edge
            self.packet_pointers[pgi] = next = ray.brick_next
            if ray.t >= 1.0:
                ray = next
                continue
            dt = self.integrate_ray(ray, pg, tf)
            # Now the ray has moved, so we grab .brick_next first, then we
            # move it to its new home
            self.send_ray_home(ray, ledges, redges, grid_neighbors, dt, 1)
            # We now are moving into a new PG, which we check for refinement
            pgn = pgs[ray.pgi]
            domega = self.get_domega(pgn.left_edge, pgn.right_edge)
            pray = &ray
            refined = self.refine_ray(pray, domega, pgn.dds[0],
                                      pgn.left_edge, pgn.right_edge)
            if refined == 0: self.append_to_packets(ray.pgi, ray)
            # At this point we can no longer access ray, as it is no longer
            # safe.
            ray2 = pray[0]
            for i in range(refined*4):
                # If we have been refined, send the ray to its appropriate
                # location.
                self.send_ray_home(ray2, ledges, redges, grid_neighbors, 0.0, 1)
                # If it wants to go back in time that is fine but it needs to
                # make sure it gets forward in time eventually
                while ray2.pgi <= pgi and ray2.t < 1.0:
                    #print "Recursing", ray2.pgi, pgi, ray2.t, ray2.nside, ray2.ipix, dt
                    # Now we grab a new set of neighbors and whatnot
                    pgn = pgs[ray2.pgi]
                    grid2_neighbors = self.find_neighbors(ray2.pgi, pgn.dds[0],
                                                          ledges, redges)
                    # We just integrate, we don't bother with the full brick
                    # integration.  This means no recursive refinement, and
                    # potential undersampling
                    dt = self.integrate_ray(ray2, pgn, tf)
                    # Now we send this ray home.  Hopefully it'll once again be
                    # forward in time.
                    self.send_ray_home(ray2, ledges, redges, grid2_neighbors,
                                       dt, 1)
                    free(grid2_neighbors)
                # This tosses us to the next one in line, of the four..
                self.append_to_packets(ray2.pgi, ray2)
                ray2 = ray2.next
            # We use this because it's been set previously.
            ray = next
            # We check to see if anything has been *added* to the queue, via a
            # send_ray_home call, here.  Otherwise we might end up in the
            # situation that the final ray is refined, thus next is NULL, but
            # there are more rays to work on because they have been added via
            # refinement.
            if ray == NULL and self.packet_pointers[pgi] != NULL:
                ray = self.packet_pointers[pgi]
                #print "Packet pointers!", ray.ipix
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
        for i in range(ledges.shape[0]):
            # Check for overlap
            if i == this_grid: continue
            if ((gle[0] <= redges[i, 0] and gre[0] >= ledges[i, 0]) and
                (gle[1] <= redges[i, 1] and gre[1] >= ledges[i, 1]) and
                (gle[2] <= redges[i, 2] and gre[2] >= ledges[i, 2])):
                count += 1
        cdef int *tr = <int *> malloc(sizeof(int) * (count + 1))
        tr[0] = count
        count = 0
        for i in range(ledges.shape[0]):
            # Check for overlap
            if i == this_grid: continue
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

    cdef int find_owner(self, AdaptiveRayPacket *ray,
                        int *neighbors,
                        np.ndarray[np.float64_t, ndim=2] ledges,
                        np.ndarray[np.float64_t, ndim=2] redges):
        cdef int i, pgi = -1
        for i in range(ledges.shape[0]):
            pgi = i
            if ((ray.pos[0] <= redges[i, 0] and ray.pos[0] >= ledges[i, 0]) and
                (ray.pos[1] <= redges[i, 1] and ray.pos[1] >= ledges[i, 1]) and
                (ray.pos[2] <= redges[i, 2] and ray.pos[2] >= ledges[i, 2])):
                    return pgi
        return -1

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
                for k in range(2):
                    r2[2] = r2[1] + (edge[k][2] - self.center[2])**2.0
                    max_r2 = fmax(max_r2, r2[2])
        domega = 4.0 * 3.1415926 * max_r2 # Used to be / Nrays
        return domega

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int refine_ray(self, AdaptiveRayPacket **pray,
                        np.float64_t domega, np.float64_t dx,
                        np.float64_t left_edge[3],
                        np.float64_t right_edge[3]):
        cdef AdaptiveRayPacket *ray = pray[0]
        cdef AdaptiveRayPacket *new_rays[4]
        cdef long Nrays = 12 * ray.nside * ray.nside
        cdef int i, j
        if domega/Nrays < dx*dx/self.rays_per_cell:
            return 0
        if ray.nside >= self.max_nside: return 0
        cdef double v_dir[3]
        # We need a record of the previous one because we're inserting into a
        # linked list.
        for i in range(4):
            new_rays[i] = <AdaptiveRayPacket *> malloc(
                            sizeof(AdaptiveRayPacket))
            new_rays[i].nside = ray.nside * 2
            new_rays[i].ipix = ray.ipix * 4 + i
            new_rays[i].t = ray.t
            healpix_interface.pix2vec_nest(
                    new_rays[i].nside, new_rays[i].ipix, v_dir)
            for j in range(3):
                new_rays[i].v_dir[j] = v_dir[j] * self.normalization
                new_rays[i].value[j] = ray.value[j]
                new_rays[i].pos[j] = self.center[j] + ray.t * new_rays[i].v_dir[j]
            new_rays[i].value[3] = ray.value[3]
        # Insert into the external list
        if ray.prev != NULL:
            ray.prev.next = new_rays[0]
        new_rays[0].prev = ray.prev
        new_rays[3].next = ray.next
        if ray.next != NULL:
            ray.next.prev = new_rays[3]
        for i in range(3):
            # Connect backward and forward
            new_rays[i].next = new_rays[i+1]
            new_rays[3-i].prev = new_rays[2-i]
        if self.first == ray:
            self.first = new_rays[0]
        self.nrays += 3
        free(ray)
        pray[0] = new_rays[0]
        return 1

# From Enzo:
#   dOmega = 4 pi r^2/Nrays
#   if (dOmega > RaysPerCell * dx^2) then split
