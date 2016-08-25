"""
Simple integrators for the radiative transfer equation



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
#cimport healpix_interface
from libc.stdlib cimport malloc, calloc, free, abs
from libc.math cimport exp, floor, log2, \
    fabs, atan, atan2, asin, cos, sin, sqrt, acos, M_PI
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip
from field_interpolation_tables cimport \
    FieldInterpolationTable, FIT_initialize_table, FIT_eval_transfer,\
    FIT_eval_transfer_with_light
from fixed_interpolator cimport *

from cython.parallel import prange, parallel, threadid
from cpython.exc cimport PyErr_CheckSignals

from .image_samplers cimport \
    ImageSampler, \
    ImageContainer, \
    VolumeRenderAccumulator

DEF Nch = 4

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int walk_volume(VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     sampler_function *sample,
                     void *data,
                     np.float64_t *return_t = NULL,
                     np.float64_t max_t = 1.0) nogil:
    cdef int cur_ind[3]
    cdef int step[3]
    cdef int x, y, i, hit, direction
    cdef np.float64_t intersect_t = 1.1
    cdef np.float64_t iv_dir[3]
    cdef np.float64_t tmax[3]
    cdef np.float64_t tdelta[3]
    cdef np.float64_t exit_t = -1.0, enter_t = -1.0
    cdef np.float64_t tl, temp_x, temp_y = -1
    if max_t > 1.0: max_t = 1.0
    direction = -1
    if vc.left_edge[0] <= v_pos[0] and v_pos[0] < vc.right_edge[0] and \
       vc.left_edge[1] <= v_pos[1] and v_pos[1] < vc.right_edge[1] and \
       vc.left_edge[2] <= v_pos[2] and v_pos[2] < vc.right_edge[2]:
        intersect_t = 0.0
        direction = 3
    for i in range(3):
        if (v_dir[i] < 0):
            step[i] = -1
        elif (v_dir[i] == 0.0):
            step[i] = 0
            continue
        else:
            step[i] = 1
        iv_dir[i] = 1.0/v_dir[i]
        if direction == 3: continue
        x = (i+1) % 3
        y = (i+2) % 3
        if step[i] > 0:
            tl = (vc.left_edge[i] - v_pos[i])*iv_dir[i]
        else:
            tl = (vc.right_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tl*v_dir[x])
        temp_y = (v_pos[y] + tl*v_dir[y])
        if fabs(temp_x - vc.left_edge[x]) < 1e-10*vc.dds[x]:
            temp_x = vc.left_edge[x]
        elif fabs(temp_x - vc.right_edge[x]) < 1e-10*vc.dds[x]:
            temp_x = vc.right_edge[x]
        if fabs(temp_y - vc.left_edge[y]) < 1e-10*vc.dds[y]:
            temp_y = vc.left_edge[y]
        elif fabs(temp_y - vc.right_edge[y]) < 1e-10*vc.dds[y]:
            temp_y = vc.right_edge[y]
        if vc.left_edge[x] <= temp_x and temp_x <= vc.right_edge[x] and \
           vc.left_edge[y] <= temp_y and temp_y <= vc.right_edge[y] and \
           0.0 <= tl and tl < intersect_t:
            direction = i
            intersect_t = tl
    if enter_t >= 0.0: intersect_t = enter_t 
    if not ((0.0 <= intersect_t) and (intersect_t < max_t)): return 0
    for i in range(3):
        # Two things have to be set inside this loop.
        # cur_ind[i], the current index of the grid cell the ray is in
        # tmax[i], the 't' until it crosses out of the grid cell
        tdelta[i] = step[i] * iv_dir[i] * vc.dds[i]
        if i == direction and step[i] > 0:
            # Intersection with the left face in this direction
            cur_ind[i] = 0
        elif i == direction and step[i] < 0:
            # Intersection with the right face in this direction
            cur_ind[i] = vc.dims[i] - 1
        else:
            # We are somewhere in the middle of the face
            temp_x = intersect_t * v_dir[i] + v_pos[i] # current position
            temp_y = ((temp_x - vc.left_edge[i])*vc.idds[i])
            # There are some really tough cases where we just within a couple
            # least significant places of the edge, and this helps prevent
            # killing the calculation through a segfault in those cases.
            if -1 < temp_y < 0 and step[i] > 0:
                temp_y = 0.0
            elif vc.dims[i] - 1 < temp_y < vc.dims[i] and step[i] < 0:
                temp_y = vc.dims[i] - 1
            cur_ind[i] =  <int> (floor(temp_y))
        if step[i] > 0:
            temp_y = (cur_ind[i] + 1) * vc.dds[i] + vc.left_edge[i]
        elif step[i] < 0:
            temp_y = cur_ind[i] * vc.dds[i] + vc.left_edge[i]
        tmax[i] = (temp_y - v_pos[i]) * iv_dir[i]
        if step[i] == 0:
            tmax[i] = 1e60
    # We have to jumpstart our calculation
    for i in range(3):
        if cur_ind[i] == vc.dims[i] and step[i] >= 0:
            return 0
        if cur_ind[i] == -1 and step[i] <= -1:
            return 0
    enter_t = intersect_t
    hit = 0
    while 1:
        hit += 1
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                i = 0
            else:
                i = 2
        else:
            if tmax[1] < tmax[2]:
                i = 1
            else:
                i = 2
        exit_t = fmin(tmax[i], max_t)
        sample(vc, v_pos, v_dir, enter_t, exit_t, cur_ind, data)
        cur_ind[i] += step[i]
        enter_t = tmax[i]
        tmax[i] += tdelta[i]
        if cur_ind[i] < 0 or cur_ind[i] >= vc.dims[i] or enter_t >= max_t:
            break
    if return_t != NULL: return_t[0] = exit_t
    return hit

def hp_pix2vec_nest(long nside, long ipix):
    raise NotImplementedError
    # cdef double v[3]
    # healpix_interface.pix2vec_nest(nside, ipix, v)
    # cdef np.ndarray[np.float64_t, ndim=1] tr = np.empty((3,), dtype='float64')
    # tr[0] = v[0]
    # tr[1] = v[1]
    # tr[2] = v[2]
    # return tr

def arr_pix2vec_nest(long nside,
                     np.ndarray[np.int64_t, ndim=1] aipix):
    raise NotImplementedError
    # cdef int n = aipix.shape[0]
    # cdef int i
    # cdef double v[3]
    # cdef long ipix
    # cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 3), dtype='float64')
    # for i in range(n):
    #     ipix = aipix[i]
    #     healpix_interface.pix2vec_nest(nside, ipix, v)
    #     tr[i,0] = v[0]
    #     tr[i,1] = v[1]
    #     tr[i,2] = v[2]
    # return tr

def hp_vec2pix_nest(long nside, double x, double y, double z):
    raise NotImplementedError
    # cdef double v[3]
    # v[0] = x
    # v[1] = y
    # v[2] = z
    # cdef long ipix
    # healpix_interface.vec2pix_nest(nside, v, &ipix)
    # return ipix

def arr_vec2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] x,
                     np.ndarray[np.float64_t, ndim=1] y,
                     np.ndarray[np.float64_t, ndim=1] z):
    raise NotImplementedError
    # cdef int n = x.shape[0]
    # cdef int i
    # cdef double v[3]
    # cdef long ipix
    # cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    # for i in range(n):
    #     v[0] = x[i]
    #     v[1] = y[i]
    #     v[2] = z[i]
    #     healpix_interface.vec2pix_nest(nside, v, &ipix)
    #     tr[i] = ipix
    # return tr

def hp_pix2ang_nest(long nside, long ipnest):
    raise NotImplementedError
    # cdef double theta, phi
    # healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
    # return (theta, phi)

def arr_pix2ang_nest(long nside, np.ndarray[np.int64_t, ndim=1] aipnest):
    raise NotImplementedError
    # cdef int n = aipnest.shape[0]
    # cdef int i
    # cdef long ipnest
    # cdef np.ndarray[np.float64_t, ndim=2] tr = np.zeros((n, 2), dtype='float64')
    # cdef double theta, phi
    # for i in range(n):
    #     ipnest = aipnest[i]
    #     healpix_interface.pix2ang_nest(nside, ipnest, &theta, &phi)
    #     tr[i,0] = theta
    #     tr[i,1] = phi
    # return tr

def hp_ang2pix_nest(long nside, double theta, double phi):
    raise NotImplementedError
    # cdef long ipix
    # healpix_interface.ang2pix_nest(nside, theta, phi, &ipix)
    # return ipix

def arr_ang2pix_nest(long nside,
                     np.ndarray[np.float64_t, ndim=1] atheta,
                     np.ndarray[np.float64_t, ndim=1] aphi):
    raise NotImplementedError
    # cdef int n = atheta.shape[0]
    # cdef int i
    # cdef long ipnest
    # cdef np.ndarray[np.int64_t, ndim=1] tr = np.zeros(n, dtype='int64')
    # cdef double theta, phi
    # for i in range(n):
    #     theta = atheta[i]
    #     phi = aphi[i]
    #     healpix_interface.ang2pix_nest(nside, theta, phi, &ipnest)
    #     tr[i] = ipnest
    # return tr

@cython.boundscheck(False)
@cython.cdivision(False)
@cython.wraparound(False)
def pixelize_healpix(long nside,
                     np.ndarray[np.float64_t, ndim=1] values,
                     long ntheta, long nphi,
                     np.ndarray[np.float64_t, ndim=2] irotation):
    raise NotImplementedError
    # # We will first to pix2vec, rotate, then calculate the angle
    # cdef int i, j, thetai, phii
    # cdef long ipix
    # cdef double v0[3], v1[3]
    # cdef double pi = 3.1415926
    # cdef np.float64_t pi2 = pi/2.0
    # cdef np.float64_t phi, theta
    # cdef np.ndarray[np.float64_t, ndim=2] results
    # cdef np.ndarray[np.int32_t, ndim=2] count
    # results = np.zeros((ntheta, nphi), dtype="float64")
    # count = np.zeros((ntheta, nphi), dtype="int32")

    # cdef np.float64_t phi0 = 0
    # cdef np.float64_t dphi = 2.0 * pi/(nphi-1)

    # cdef np.float64_t theta0 = 0
    # cdef np.float64_t dtheta = pi/(ntheta-1)
    # # We assume these are the rotated theta and phi
    # for thetai in range(ntheta):
    #     theta = theta0 + dtheta * thetai
    #     for phii in range(nphi):
    #         phi = phi0 + dphi * phii
    #         # We have our rotated vector
    #         v1[0] = cos(phi) * sin(theta)
    #         v1[1] = sin(phi) * sin(theta)
    #         v1[2] = cos(theta)
    #         # Now we rotate back
    #         for i in range(3):
    #             v0[i] = 0
    #             for j in range(3):
    #                 v0[i] += v1[j] * irotation[j,i]
    #         # Get the pixel this vector is inside
    #         healpix_interface.vec2pix_nest(nside, v0, &ipix)
    #         results[thetai, phii] = values[ipix]
    #         count[i, j] += 1
    # return results, count

def healpix_aitoff_proj(np.ndarray[np.float64_t, ndim=1] pix_image,
                        long nside,
                        np.ndarray[np.float64_t, ndim=2] image,
                        np.ndarray[np.float64_t, ndim=2] irotation):
    raise NotImplementedError
    # cdef double pi = np.pi
    # cdef int i, j, k, l
    # cdef np.float64_t x, y, z, zb
    # cdef np.float64_t dx, dy, inside
    # cdef double v0[3], v1[3]
    # dx = 2.0 / (image.shape[1] - 1)
    # dy = 2.0 / (image.shape[0] - 1)
    # cdef np.float64_t s2 = sqrt(2.0)
    # cdef long ipix
    # for i in range(image.shape[1]):
    #     x = (-1.0 + i*dx)*s2*2.0
    #     for j in range(image.shape[0]):
    #         y = (-1.0 + j * dy)*s2
    #         zb = (x*x/8.0 + y*y/2.0 - 1.0)
    #         if zb > 0: continue
    #         z = (1.0 - (x/4.0)**2.0 - (y/2.0)**2.0)
    #         z = z**0.5
    #         # Longitude
    #         phi = (2.0*atan(z*x/(2.0 * (2.0*z*z-1.0))) + pi)
    #         # Latitude
    #         # We shift it into co-latitude
    #         theta = (asin(z*y) + pi/2.0)
    #         # Now to account for rotation we translate into vectors
    #         v1[0] = cos(phi) * sin(theta)
    #         v1[1] = sin(phi) * sin(theta)
    #         v1[2] = cos(theta)
    #         for k in range(3):
    #             v0[k] = 0
    #             for l in range(3):
    #                 v0[k] += v1[l] * irotation[l,k]
    #         healpix_interface.vec2pix_nest(nside, v0, &ipix)
    #         #print "Rotated", v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], ipix, pix_image[ipix]
    #         image[j, i] = pix_image[ipix]

def arr_fisheye_vectors(int resolution, np.float64_t fov, int nimx=1, int
        nimy=1, int nimi=0, int nimj=0, np.float64_t off_theta=0.0, np.float64_t
        off_phi=0.0):
    # We now follow figures 4-7 of:
    # http://paulbourke.net/miscellaneous/domefisheye/fisheye/
    # ...but all in Cython.
    cdef np.ndarray[np.float64_t, ndim=3] vp
    cdef int i, j
    cdef np.float64_t r, phi, theta, px, py
    cdef np.float64_t fov_rad = fov * np.pi / 180.0
    cdef int nx = resolution/nimx
    cdef int ny = resolution/nimy
    vp = np.zeros((nx,ny, 3), dtype="float64")
    for i in range(nx):
        px = (2.0 * (nimi*nx + i)) / resolution - 1.0
        for j in range(ny):
            py = (2.0 * (nimj*ny + j)) / resolution - 1.0
            r = (px*px + py*py)**0.5
            if r > 1.01:
                vp[i,j,0] = vp[i,j,1] = vp[i,j,2] = 0.0
                continue
            phi = atan2(py, px)
            theta = r * fov_rad / 2.0
            theta += off_theta
            phi += off_phi
            vp[i,j,0] = sin(theta) * cos(phi)
            vp[i,j,1] = sin(theta) * sin(phi)
            vp[i,j,2] = cos(theta)
    return vp


