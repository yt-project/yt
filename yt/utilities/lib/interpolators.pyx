
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
@cython.boundscheck(False)
def ghost_zone_interpolate(int rf,
                           np.ndarray[np.float64_t, ndim=3] input_field,
                           np.ndarray[np.float64_t, ndim=1] input_left,
                           np.ndarray[np.float64_t, ndim=3] output_field,
                           np.ndarray[np.float64_t, ndim=1] output_left):
    cdef int oi, oj, ok
    cdef int ii, ij, ik
    cdef np.float64_t xp, xm, yp, ym, zp, zm, temp
    cdef np.float64_t ods[3]
    cdef np.float64_t ids[3]
    cdef np.float64_t iids[3]
    cdef np.float64_t opos[3]
    cdef np.float64_t ropos[3]
    cdef int i
    for i in range(3):
        temp = input_left[i] + (rf * (input_field.shape[i] - 1))
        ids[i] = (temp - input_left[i])/(input_field.shape[i]-1)
        temp = output_left[i] + output_field.shape[i] - 1
        ods[i] = (temp - output_left[i])/(output_field.shape[i]-1)
        iids[i] = 1.0/ids[i]
    opos[0] = output_left[0]
    for oi in range(output_field.shape[0]):
        ropos[0] = ((opos[0] - input_left[0]) * iids[0])
        ii = iclip(<int> ropos[0], 0, input_field.shape[0] - 2)
        xp = ropos[0] - ii
        xm = 1.0 - xp
        opos[1] = output_left[1]
        for oj in range(output_field.shape[1]):
            ropos[1] = ((opos[1] - input_left[1]) * iids[1])
            ij = iclip(<int> ropos[1], 0, input_field.shape[1] - 2)
            yp = ropos[1] - ij
            ym = 1.0 - yp
            opos[2] = output_left[2]
            for ok in range(output_field.shape[2]):
                ropos[2] = ((opos[2] - input_left[2]) * iids[2])
                ik = iclip(<int> ropos[2], 0, input_field.shape[2] - 2)
                zp = ropos[2] - ik
                zm = 1.0 - zp
                output_field[oi,oj,ok] = \
                     input_field[ii  ,ij  ,ik  ] * (xm*ym*zm) \
                   + input_field[ii+1,ij  ,ik  ] * (xp*ym*zm) \
                   + input_field[ii  ,ij+1,ik  ] * (xm*yp*zm) \
                   + input_field[ii  ,ij  ,ik+1] * (xm*ym*zp) \
                   + input_field[ii+1,ij  ,ik+1] * (xp*ym*zp) \
                   + input_field[ii  ,ij+1,ik+1] * (xm*yp*zp) \
                   + input_field[ii+1,ij+1,ik  ] * (xp*yp*zm) \
                   + input_field[ii+1,ij+1,ik+1] * (xp*yp*zp)
                opos[2] += ods[2]
            opos[1] += ods[1]
        opos[0] += ods[0]

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def replace_nonperiodic_with_extrap(np.ndarray[np.float64_t, ndim=3] field,
                                    np.ndarray[np.int_t,     ndim=1] ncell_lo,
                                    np.ndarray[np.int_t,     ndim=1] ncell_hi):
    """
    The purpose of this function is to use linear extrapolation of
    simulation data to approximate data just outside of the non-periodic
    domain. Prior to evaluating this function, the array `field` contains data
    that was obtained by assuming periodicity of the computational domain. This
    function replaces the falsely-assumed periodic cells, indicated by
    `ncell_lo` (for lower end of the array) and `ncell_hi` (for upper end of
    the array), with values obtained via extrapolation of the data that was
    not falsely-assumed periodic.

    Parameters
    ----------
    field : 3D array
        Three-dimensional data of the flowfield
    ncell_lo : 1D array of length 3
        Number of cells to replace on the low side along each axis
    ncell_hi : 1D array of length 3
        Number of cells to replace on the high side along each axis

    This is best illustrated through a one-dimensional example. Consider a
    linear function: f(x) = x , x \in [0,4], and the `field` variable contains
    values that were obtained by falsely-assuming periodicity of f(x). This
    would result in the following data:

    -----------------------------------------------
    index:        0  |  1  |  2  |  3  |  4  |  5
    -----------------------------------------------
    x or f(x):  -0.5 | 0.5 | 1.5 | 2.5 | 3.5 | 4.5
    -----------------------------------------------
    field:       3.5 | 0.5 | 1.5 | 2.5 | 3.5 | 0.5
    -----------------------------------------------

    where field.shape ---> (6,1,1), ncell_lo ---> (1,0,0), and
    ncell_hi ---> (1,0,0) . For this example, as indicated by `ncell_lo`, we
    replace the value of field[0,0,0] using linear extrapolation from
    field[1,0,0] and field[2,0,0]. We would then apply a similar technique for
    the upper boundary. Evaluating this function then yields the following
    output of `field`:

    -----------------------------------------------
    field:      -0.5 | 0.5 | 1.5 | 2.5 | 3.5 | 4.5
    -----------------------------------------------

    This is particularly important for estimating gradients near non-periodic
    boundaries in the domain.

    """
    cdef int i, j, k
    cdef int ii, ij, ik
    cdef int nx, ny, nz
    cdef int lx, ly, lz
    cdef int hx, hy, hz
    cdef np.float64_t xp, xm, yp, ym, zp, zm
    cdef np.float64_t qp, qm

    nx, ny, nz = (field.shape[0], field.shape[1], field.shape[2])
    lx, ly, lz = ncell_lo
    hx, hy, hz = ncell_hi

    # Extrapolate in the x direction
    for j in range(ly, ny-hy):
        for k in range(lz, nz-hz):
            # Fix xlo
            qm = field[lx,   j, k]
            qp = field[lx+1, j, k]
            for i in range(lx):
                xp = float(i - lx)
                xm = 1.0 - xp
                field[i, j, k] = xm*qm + xp*qp
            # Fix xhi
            qm = field[nx-hx-2, j, k]
            qp = field[nx-hx-1, j, k]
            for i in range(nx-hx, nx):
                xp = float(i - (nx-hx-2))
                xm = 1.0 - xp
                field[i, j, k] = xm*qm + xp*qp

    # Extrapolate in the y direction
    for i in range(nx):
        for k in range(lz, nz-hz):
            # Fix ylo
            qm = field[i, ly,   k]
            qp = field[i, ly+1, k]
            for j in range(ly):
                yp = float(j - ly)
                ym = 1.0 - yp
                field[i, j, k] = ym*qm + yp*qp
            # Fix yhi
            qm = field[i, ny-hy-2, k]
            qp = field[i, ny-hy-1, k]
            for j in range(ny-hy, ny):
                yp = float(j - (ny-hy-2))
                ym = 1.0 - yp
                field[i, j, k] = ym*qm + yp*qp

    # Extrapolate in the z direction
    for i in range(nx):
        for j in range(ny):
            # Fix zlo
            qm = field[i, j, lz]
            qp = field[i, j, lz+1]
            for k in range(lz):
                zp = float(k - lz)
                zm = 1.0 - zp
                field[i, j, k] = zm*qm + zp*qp
            # Fix zhi
            qm = field[i, j, nz-hz-2]
            qp = field[i, j, nz-hz-1]
            for k in range(nz-hz, nz):
                zp = float(k - (nz-hz-2))
                zm = 1.0 - zp
                field[i, j, k] = zm*qm + zp*qp
