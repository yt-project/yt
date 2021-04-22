
# distutils: libraries = STD_LIBS
"""
Simple interpolators



"""


import numpy as np

cimport cython
cimport numpy as np

from yt.utilities.lib.fp_utils cimport fclip, fmax, fmin, iclip, imax, imin


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def UnilinearlyInterpolate(np.ndarray[np.float64_t, ndim=1] table,
                           np.ndarray[np.float64_t, ndim=1] x_vals,
                           np.ndarray[np.float64_t, ndim=1] x_bins,
                           np.ndarray[np.int32_t, ndim=1] x_is,
                           np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef int i, x_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp = (x - x_bins[x_i]) * dx_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        output[i]  = table[x_i  ] * (xm) \
                   + table[x_i+1] * (xp)

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def BilinearlyInterpolate(np.ndarray[np.float64_t, ndim=2] table,
                          np.ndarray[np.float64_t, ndim=1] x_vals,
                          np.ndarray[np.float64_t, ndim=1] y_vals,
                          np.ndarray[np.float64_t, ndim=1] x_bins,
                          np.ndarray[np.float64_t, ndim=1] y_bins,
                          np.ndarray[np.int32_t, ndim=1] x_is,
                          np.ndarray[np.int32_t, ndim=1] y_is,
                          np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef double y, yp, ym
    cdef double dx_inv, dy_inv
    cdef int i, x_i, y_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        y_i = y_is[i]
        x = x_vals[i]
        y = y_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        dy_inv = 1.0 / (y_bins[y_i+1] - y_bins[y_i])
        xp = (x - x_bins[x_i]) * dx_inv
        yp = (y - y_bins[y_i]) * dy_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        ym = (y_bins[y_i+1] - y) * dy_inv
        output[i]  = table[x_i  , y_i  ] * (xm*ym) \
                   + table[x_i+1, y_i  ] * (xp*ym) \
                   + table[x_i  , y_i+1] * (xm*yp) \
                   + table[x_i+1, y_i+1] * (xp*yp)

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def TrilinearlyInterpolate(np.ndarray[np.float64_t, ndim=3] table,
                           np.ndarray[np.float64_t, ndim=1] x_vals,
                           np.ndarray[np.float64_t, ndim=1] y_vals,
                           np.ndarray[np.float64_t, ndim=1] z_vals,
                           np.ndarray[np.float64_t, ndim=1] x_bins,
                           np.ndarray[np.float64_t, ndim=1] y_bins,
                           np.ndarray[np.float64_t, ndim=1] z_bins,
                           np.ndarray[np.int_t, ndim=1] x_is,
                           np.ndarray[np.int_t, ndim=1] y_is,
                           np.ndarray[np.int_t, ndim=1] z_is,
                           np.ndarray[np.float64_t, ndim=1] output):
    cdef double x, xp, xm
    cdef double y, yp, ym
    cdef double z, zp, zm
    cdef double dx_inv, dy_inv, dz_inv
    cdef int i, x_i, y_i, z_i
    for i in range(x_vals.shape[0]):
        x_i = x_is[i]
        y_i = y_is[i]
        z_i = z_is[i]
        x = x_vals[i]
        y = y_vals[i]
        z = z_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        dy_inv = 1.0 / (y_bins[y_i+1] - y_bins[y_i])
        dz_inv = 1.0 / (z_bins[z_i+1] - z_bins[z_i])
        xp = (x - x_bins[x_i]) * dx_inv
        yp = (y - y_bins[y_i]) * dy_inv
        zp = (z - z_bins[z_i]) * dz_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        ym = (y_bins[y_i+1] - y) * dy_inv
        zm = (z_bins[z_i+1] - z) * dz_inv
        output[i]  = table[x_i  ,y_i  ,z_i  ] * (xm*ym*zm) \
                   + table[x_i+1,y_i  ,z_i  ] * (xp*ym*zm) \
                   + table[x_i  ,y_i+1,z_i  ] * (xm*yp*zm) \
                   + table[x_i  ,y_i  ,z_i+1] * (xm*ym*zp) \
                   + table[x_i+1,y_i  ,z_i+1] * (xp*ym*zp) \
                   + table[x_i  ,y_i+1,z_i+1] * (xm*yp*zp) \
                   + table[x_i+1,y_i+1,z_i  ] * (xp*yp*zm) \
                   + table[x_i+1,y_i+1,z_i+1] * (xp*yp*zp)

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(True)
def ghost_zone_interpolate(
    int rf,
    int order,
    np.ndarray[np.float64_t, ndim=3] input_field,
    np.ndarray[np.float64_t, ndim=1] input_left,
    np.ndarray[np.float64_t, ndim=3] output_field,
    np.ndarray[np.float64_t, ndim=1] output_left
):
    # Counting indices
    cdef int io, jo, ko
    cdef int ii, ji, ki
    cdef int i, j, k
    cdef int i0, j0, k0
    # dx of the input and output fields
    cdef double dxi = 1.0
    cdef double dxo = dxi / rf
    # Number of cells in the input and output grid
    cdef int[3] nxi = [input_field.shape[0],
                       input_field.shape[1],
                       input_field.shape[2]]
    cdef int[3] nxo = [output_field.shape[0],
                       output_field.shape[1],
                       output_field.shape[2]]
    # Number of cells in the interpolation stencil
    cdef int[3] nxs = [order, order, order]
    if order == 0:  nxs = [1,1,1]
    if order == -1: nxs = [6,6,6]
    if order == -2: nxs = nxi
    # Determine the index of the stencil that is left of the center
    cdef int[3] ilc = [nxs[0]//2 - 1,
                       nxs[1]//2 - 1,
                       nxs[2]//2 - 1]
    # Determine maximum ii, ji, ki values that are allowed
    cdef int ii_max = nxi[0]-nxs[0]
    cdef int ji_max = nxi[1]-nxs[1]
    cdef int ki_max = nxi[2]-nxs[2]
    # Node locations for interpolation
    cdef double[::1] xn = np.empty((nxs[0]))
    cdef double[::1] yn = np.empty((nxs[1]))
    cdef double[::1] zn = np.empty((nxs[2]))
    # Weights at each of the nodes
    cdef double[::1] xw = np.empty((nxs[0]))
    cdef double[::1] yw = np.empty((nxs[1]))
    cdef double[::1] zw = np.empty((nxs[2]))
    # Current and initial (0th) index position with respect to input grid
    cdef double iposi,  jposi,  kposi
    cdef double iposi0, jposi0, kposi0
    # Intermediate interpolated fields
    cdef double[:,::1] xyfieldi = np.empty((nxi[0],nxi[1]))
    cdef double[::1] xfieldi = np.empty((nxi[0]))
    cdef double[::1] xfieldo = np.empty((nxo[0]))

    # Compute the leftmost cell-center index position in the grid
    iposi0 = (output_left[0]+0.5/rf) - (input_left[0]+0.5) - ilc[0]
    jposi0 = (output_left[1]+0.5/rf) - (input_left[1]+0.5) - ilc[1]
    kposi0 = (output_left[2]+0.5/rf) - (input_left[2]+0.5) - ilc[2]

    # Linear interpolation ----------------------------------------------------
    if order == 2:
        kposi = kposi0
        for ko in range(nxo[2]):
            ki = <int>kposi
            zw[1] = kposi - ki
            zw[0] = 1.0 - zw[1]

            # Interpolate to the x-y plane
            for ii in range(nxi[0]):
                for ji in range(nxi[1]):
                    xyfieldi[ii,ji] = zw[0] * input_field[ii, ji, ki] \
                                    + zw[1] * input_field[ii, ji, ki+1]

            jposi = jposi0
            for jo in range(nxo[1]):
                ji = <int>jposi
                yw[1] = jposi - ji
                yw[0] = 1.0 - yw[1]

                # Interpolate to the x line
                for ii in range(nxi[0]):
                    xfieldi[ii] = yw[0] * xyfieldi[ii, ji] \
                                + yw[1] * xyfieldi[ii, ji+1]

                iposi = iposi0
                for io in range(nxo[0]):
                    ii = <int>iposi
                    xw[1] = iposi - ii
                    xw[0] = 1.0 - xw[1]

                    # Interpolate to points in output field
                    xfieldo[io] = xw[0] * xfieldi[ii]  \
                                + xw[1] * xfieldi[ii+1]

                    iposi += dxo
                output_field[:, jo, ko] = xfieldo
                jposi += dxo
            kposi += dxo

    # Cubic and quintic interpolation -----------------------------------------
    elif order > 2:
        kposi = kposi0
        for ko in range(nxo[2]):
            ki = iclip(<int>kposi, 0, ki_max)
            for k in range(nxs[2]):
                zn[k] = ki + k - ilc[2]
            zw = lagrange_weights(zn, kposi)

            # Interpolate to the x-y plane
            for ii in range(nxi[0]):
                for ji in range(nxi[1]):
                    xyfieldi[ii,ji] = 0.0
                    for k0, k in enumerate(range(ki, ki + nxs[2])):
                        xyfieldi[ii,ji] += zw[k0] * input_field[ii, ji, k]

            jposi = jposi0
            for jo in range(nxo[1]):
                ji = iclip(<int>jposi, 0, ji_max)
                for j in range(nxs[1]):
                    yn[j] = ji + j - ilc[1]
                yw = lagrange_weights(yn, jposi)

                # Interpolate to the x line
                for ii in range(nxi[0]):
                    xfieldi[ii] = 0.0
                    for j0, j in enumerate(range(ji, ji + nxs[1])):
                        xfieldi[ii] += yw[j0] * xyfieldi[ii, j]

                iposi = iposi0
                for io in range(nxo[0]):
                    ii = iclip(<int>iposi, 0, ii_max)
                    for i in range(nxs[0]):
                        xn[i] = ii + i - ilc[0]
                    xw = lagrange_weights(xn, iposi)

                    # Interpolate to points in output field
                    xfieldo[io] = 0.0
                    for i0, i in enumerate(range(ii, ii + nxs[0])):
                        xfieldo[io] += xw[i0] * xfieldi[i]

                    iposi += dxo
                output_field[:, jo, ko] = xfieldo
                jposi += dxo
            kposi += dxo

    # Akima spline interpolation ----------------------------------------------
    elif order == -1:
        kposi = kposi0
        for ko in range(nxo[2]):
            ki = iclip(<int>kposi, 0, ki_max)
            for k in range(nxs[2]):
                zn[k] = ki + k - ilc[2]

            # Interpolate to the x-y plane
            for ii in range(nxi[0]):
                for ji in range(nxi[1]):
                    xyfieldi[ii,ji] = akima_interp(
                        zn,
                        input_field[ii, ji, ki:ki+nxs[2]],
                        kposi
                    )

            # Interpolate to the x-y plane
            # zw = lagrange_weights(zn, kposi)
            # for ii in range(nxi[0]):
            #     for ji in range(nxi[1]):
            #         xyfieldi[ii,ji] = 0.0
            #         for k0, k in enumerate(range(ki, ki + nxs[2])):
            #             xyfieldi[ii,ji] += zw[k0] * input_field[ii, ji, k]

            jposi = jposi0
            for jo in range(nxo[1]):
                ji = iclip(<int>jposi, 0, ji_max)
                for j in range(nxs[1]):
                    yn[j] = ji + j - ilc[1]

                # Interpolate to the x line
                for ii in range(nxi[0]):
                    xfieldi[ii] = akima_interp(
                        yn,
                        xyfieldi[ii, ji:ji+nxs[1]],
                        jposi
                    )

                # Interpolate to the x line
                # yw = lagrange_weights(yn, jposi)
                # for ii in range(nxi[0]):
                #     xfieldi[ii] = 0.0
                #     for j0, j in enumerate(range(ji, ji + nxs[1])):
                #         xfieldi[ii] += yw[j0] * xyfieldi[ii, j]

                iposi = iposi0
                for io in range(nxo[0]):
                    ii = iclip(<int>iposi, 0, ii_max)
                    for i in range(nxs[0]):
                        xn[i] = ii + i - ilc[0]

                    # Interpolate to points in output field
                    xfieldo[io] = akima_interp(
                        xn,
                        xfieldi[ii:ii+nxs[0]],
                        iposi
                    )

                    # Interpolate to points in output field
                    # xw = lagrange_weights(xn, iposi)
                    # xfieldo[io] = 0.0
                    # for i0, i in enumerate(range(ii, ii + nxs[0])):
                    #     xfieldo[io] += xw[i0] * xfieldi[i]


                    iposi += dxo
                output_field[:, jo, ko] = xfieldo
                jposi += dxo
            kposi += dxo


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(True)
cdef lagrange_weights(double[::1] x,
                      double xp):
    cdef int nx = x.size
    cdef double[::1] wgts = np.empty((nx))
    cdef double c1, c2, c3, c4, c5
    cdef int i, j
    wgts[0] = 1.0
    c1 = 1.0
    c4 = x[0] - xp
    for i in range(1,nx):
        c2 = 1.0
        c5 = c4
        c4 = x[i] - xp
        for j in range(i):
            c3 = x[i] - x[j]
            c2 = c2 * c3
            if j == i - 1:
                wgts[i] = -c1 * c5 * wgts[i-1] / c2
            wgts[j] = c4 * wgts[j] / c3
        c1 = c2
    return wgts

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(True)
cdef akima_interp(double[::1] x,
                  double[::1] f,
                  double xp):
    cdef int nx = x.size
    cdef double[::1] m = np.empty((nx+3))
    cdef double[::1] dm = np.empty((nx+2))
    cdef double[::1] t = np.empty((nx))
    cdef double[::1] p1 = np.empty((nx-1))
    cdef double[::1] p2 = np.empty((nx-1))
    cdef double[::1] p3 = np.empty((nx-1))
    cdef double a, b, ab
    cdef int i, idx

    # Compute derivative with extra conditions (note, dx==1 so not included)
    for i in range(2,nx+1):
        m[i] = f[i-1] - f[i-2]
    m[1] = 2.0*m[2] - m[3]
    m[0] = 2.0*m[1] - m[2]
    m[nx+1] = 2.0*m[nx]   - m[nx-1]
    m[nx+2] = 2.0*m[nx+1] - m[nx]
    # Compute differences
    for i in range(nx+2):
        dm[i] = m[i+1] - m[i]
    # Compute t at nodal values
    for i in range(nx):
        a = dm[i+2]
        b = dm[i]
        ab = a + b
        if ab > 0.0:
            t[i] = (a*m[i+1] + b*m[i+2]) / ab
        else:
            t[i] = 0.5*(m[i+1] + m[i+2])
    # Compute polynomial coefficients between nodes
    for i in range(nx-1):
        p1[i] = t[i]
        p2[i] = 3.0*m[i+2] - 2.0*t[i] - t[i+1]
        p3[i] = t[i] + t[i+1] - 2.0*m[i+2]
    # Determine value at xp
    idx = <int>(xp-x[0])
    w = xp - x[idx]
    fp = f[idx] + p1[idx]*w + p2[idx]*w**2 + p3[idx]*w**3

    return fp
