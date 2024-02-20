cimport cython
import numpy as np
cimport numpy as np

from libc.math cimport cos, sin, atan2, acos, sqrt

from numpy.math cimport PI as NPY_PI
from numpy.math cimport INFINITY as NPY_INF

from yt.utilities.lib.fp_utils cimport fmax, fmin


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef (np.float64_t, np.float64_t, np.float64_t) spherical_to_cartesian(np.float64_t r,
                           np.float64_t theta,
                           np.float64_t phi) noexcept nogil:
        # transform a single point in spherical coords to cartesian
        # r : radius
        # theta: colatitude
        # phi: azimuthal (longitudinal) angle
        cdef np.float64_t x, y, xy, z

        if r == 0.0:
            return 0.0, 0.0, 0.0

        xy = r * sin(theta)
        x = xy * cos(phi)
        y = xy * sin(phi)
        z = r * cos(theta)
        return x, y, z


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef (np.float64_t, np.float64_t, np.float64_t) cartesian_to_spherical(np.float64_t x,
                           np.float64_t y,
                           np.float64_t z) noexcept nogil:
        # transform a single point in cartesian coords to spherical, returns
        # r : radius
        # theta: colatitude
        # phi: azimuthal angle in range (0, 2pi)
        cdef np.float64_t r, theta, phi
        r = sqrt(x*x + y*y + z*z)
        theta = acos(z / r)
        phi = atan2(y, x)
        # atan2 returns -pi to pi, adjust to (0, 2pi)
        if phi < 0:
            phi = phi + 2 * NPY_PI
        return r, theta, phi


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def cartesian_points_to_spherical(np.float64_t[:] x,
                                  np.float64_t[:] y,
                                  np.float64_t[:] z):
        # transform an array of points in cartesian coords to spherical, returns
        # r : radius
        # theta: colatitude
        # phi: azimuthal angle in range (0, 2pi)
        cdef np.ndarray[np.float64_t, ndim=1] r, theta, phi
        cdef int i, n_x

        r = np.zeros(x.shape[0], dtype='float64')
        theta = np.zeros(x.shape[0], dtype='float64')
        phi = np.zeros(x.shape[0], dtype='float64')
        n_x = x.size
        with nogil:
            for i in range(n_x):
                r[i], theta[i], phi[i] = cartesian_to_spherical(x[i], y[i], z[i])

        return r, theta, phi


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def spherical_points_to_cartesian(np.float64_t[:] r,
                                  np.float64_t[:] theta,
                                  np.float64_t[:] phi):
        # transform an array of points in spherical coords to cartesian
        cdef np.ndarray[np.float64_t, ndim=1] x, y, z
        cdef int i, n_r

        x = np.zeros(r.shape[0], dtype='float64')
        y = np.zeros(r.shape[0], dtype='float64')
        z = np.zeros(r.shape[0], dtype='float64')
        n_r = r.size

        with nogil:
            for i in range(n_r):
                x[i], y[i], z[i] = spherical_to_cartesian(r[i], theta[i], phi[i])

        return x, y, z


cdef class MixedCoordBBox:
    # abstract class for calculating cartesian bounding boxes
    # from non-cartesian grid elements.
    cdef void get_cartesian_bbox(self,
                                np.float64_t pos0,
                                np.float64_t pos1,
                                np.float64_t pos2,
                                np.float64_t dpos0,
                                np.float64_t dpos1,
                                np.float64_t dpos2,
                                np.float64_t xyz_i[3],
                                np.float64_t dxyz_i[3]
                                ) noexcept nogil:
        pass


cdef class SphericalMixedCoordBBox(MixedCoordBBox):
    # Cartesian bounding boxes of spherical grid elements
    cdef void get_cartesian_bbox(self,
                        np.float64_t pos0,
                        np.float64_t pos1,
                        np.float64_t pos2,
                        np.float64_t dpos0,
                        np.float64_t dpos1,
                        np.float64_t dpos2,
                        np.float64_t xyz_i[3],
                        np.float64_t dxyz_i[3]
                        ) noexcept nogil:

        cdef np.float64_t r_i, theta_i, phi_i, dr_i, dtheta_i, dphi_i
        cdef np.float64_t h_dphi, h_dtheta, h_dr, r_r
        cdef np.float64_t xi, yi, zi, r_lr, theta_lr, phi_lr, phi_lr2, theta_lr2
        cdef np.float64_t xli, yli, zli, xri, yri, zri, r_xy, r_xy2
        cdef int isign_r, isign_ph, isign_th
        cdef np.float64_t sign_r, sign_th, sign_ph

        cdef np.float64_t NPY_PI_2 = NPY_PI / 2.0
        cdef np.float64_t NPY_PI_3_2 = 3. * NPY_PI / 2.0
        cdef np.float64_t NPY_2xPI = 2. * NPY_PI

        r_i = pos0
        theta_i = pos1
        phi_i = pos2
        dr_i = dpos0
        dtheta_i = dpos1
        dphi_i = dpos2

        # initialize the left/right values
        xli = NPY_INF
        yli = NPY_INF
        zli = NPY_INF
        xri = -1.0 * NPY_INF
        yri = -1.0 * NPY_INF
        zri = -1.0 * NPY_INF

        # find the min/max bounds over the 8 corners of the
        # spherical volume element.
        h_dphi =  dphi_i / 2.0
        h_dtheta =  dtheta_i / 2.0
        h_dr =  dr_i / 2.0
        for isign_r in range(2):
            for isign_ph in range(2):
                for isign_th in range(2):
                    sign_r = 1.0 - 2.0 * <float> isign_r
                    sign_th = 1.0 - 2.0 * <float> isign_th
                    sign_ph = 1.0 - 2.0 * <float> isign_ph
                    r_lr = r_i + sign_r * h_dr
                    theta_lr = theta_i + sign_th * h_dtheta
                    phi_lr = phi_i + sign_ph * h_dphi

                    xi, yi, zi = spherical_to_cartesian(r_lr, theta_lr, phi_lr)

                    xli = fmin(xli, xi)
                    yli = fmin(yli, yi)
                    zli = fmin(zli, zi)
                    xri = fmax(xri, xi)
                    yri = fmax(yri, yi)
                    zri = fmax(zri, zi)

        # need to correct for special cases:
        # if polar angle, phi, spans pi/2, pi or 3pi/2 then just
        # taking the min/max of the corners will miss the cusp of the
        # element. When this condition is met, the x/y min/max will
        # equal +/- the projection of the max r in the xy plane -- in this case,
        # the theta angle that gives the max projection of r in
        # the x-y plane will change depending on the whether theta < or > pi/2,
        # so the following calculates for the min/max theta value of the element
        # and takes the max.
        # ALSO note, that the following does check for when an edge aligns with the
        # phi=0/2pi, it does not handle an element spanning the periodic boundary.
        # Oh and this may break down for large elements that span whole
        # quadrants...
        phi_lr =  phi_i - h_dphi
        phi_lr2 = phi_i + h_dphi
        theta_lr = theta_i - h_dtheta
        theta_lr2 = theta_i + h_dtheta
        r_r = r_i + h_dr
        if theta_lr < NPY_PI_2 and theta_lr2 > NPY_PI_2:
            r_xy = r_r
        else:
            r_xy = r_r * sin(theta_lr)
            r_xy2 = r_r * sin(theta_lr2)
            r_xy = fmax(r_xy, r_xy2)

        if phi_lr == 0.0 or phi_lr2 == NPY_2xPI:
            # need to re-check this, for when theta spans equator
            xri = r_xy
        elif phi_lr < NPY_PI_2 and phi_lr2  > NPY_PI_2:
            yri = r_xy
        elif phi_lr < NPY_PI and phi_lr2  > NPY_PI:
            xli = -r_xy
        elif phi_lr < NPY_PI_3_2 and phi_lr2  > NPY_PI_3_2:
            yli = -r_xy

        xyz_i[0] = (xri+xli)/2.
        xyz_i[1] = (yri+yli)/2.
        xyz_i[2] = (zri+zli)/2.
        dxyz_i[0] = xri-xli
        dxyz_i[1] = yri-yli
        dxyz_i[2] = zri-zli


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def cartesian_bboxes(MixedCoordBBox bbox_handler,
                      np.float64_t[:] pos0,
                      np.float64_t[:] pos1,
                      np.float64_t[:] pos2,
                      np.float64_t[:] dpos0,
                      np.float64_t[:] dpos1,
                      np.float64_t[:] dpos2,
                      np.float64_t[:] x,
                      np.float64_t[:] y,
                      np.float64_t[:] z,
                      np.float64_t[:] dx,
                      np.float64_t[:] dy,
                      np.float64_t[:] dz,
                                      ):
    # calculates the cartesian bounding boxes around non-cartesian
    # volume elements
    #
    # bbox_handler : a MixedCoordBBox child instance
    # pos0, pos1, pos2: native coordinates of element centers
    # dpos0, dpos1, dpos2: element widths in native coordinates
    # x, y, z: cartesian centers of bounding boxes, modified in place
    # dx, dy, dz : full-widths of the cartesian bounding boxes, modified in place

    cdef int i, n_pos
    cdef np.float64_t xyz_i[3]
    cdef np.float64_t dxyz_i[3]

    n_pos = pos0.size
    with nogil:
        for i in range(n_pos):

            bbox_handler.get_cartesian_bbox(pos0[i],
                                            pos1[i],
                                            pos2[i],
                                            dpos0[i],
                                            dpos1[i],
                                            dpos2[i],
                                            xyz_i,
                                            dxyz_i)

            x[i] = xyz_i[0]
            y[i] = xyz_i[1]
            z[i] = xyz_i[2]
            dx[i] = dxyz_i[0]
            dy[i] = dxyz_i[1]
            dz[i] = dxyz_i[2]
