# distutils: include_dirs = LIB_DIR
# distutils: extra_compile_args = CPP14_FLAG OMP_ARGS
# distutils: extra_link_args = CPP14_FLAG OMP_ARGS
# distutils: language = c++
# distutils: libraries = STD_LIBS
# distutils: sources = yt/utilities/lib/pixelization_constants.cpp
"""
Pixelization routines



"""


import numpy as np

cimport cython
cimport libc.math as math
cimport numpy as np
from cython.view cimport array as cvarray

from yt.utilities.lib.fp_utils cimport (
    fabs,
    fmax,
    fmin,
    i64max,
    i64min,
    iclip,
    imax,
    imin,
)

from yt.utilities.exceptions import YTElementTypeNotRecognized, YTPixelizeError

from cpython.exc cimport PyErr_CheckSignals
from cython.parallel cimport parallel, prange
from libc.stdlib cimport free, malloc
from vec3_ops cimport cross, dot, subtract

from yt.geometry.particle_deposit cimport get_kernel_func, kernel_func
from yt.utilities.lib.element_mappings cimport (
    ElementSampler,
    P1Sampler1D,
    P1Sampler2D,
    P1Sampler3D,
    Q1Sampler2D,
    Q1Sampler3D,
    Q2Sampler2D,
    S2Sampler3D,
    T2Sampler2D,
    Tet2Sampler3D,
    W1Sampler3D,
)

from yt.funcs import get_pbar

from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue
from yt.utilities.lib.cykdtree.kdtree cimport KDTree, Node, PyKDTree, uint32_t, uint64_t
from yt.utilities.lib.particle_kdtree_tools cimport (
    axes_range,
    find_neighbors,
    set_axes_range,
)


cdef int TABLE_NVALS=512

cdef extern from "pixelization_constants.hpp":
    enum:
        MAX_NUM_FACES

    int HEX_IND
    int HEX_NF
    np.uint8_t hex_face_defs[MAX_NUM_FACES][2][2]

    int TETRA_IND
    int TETRA_NF
    np.uint8_t tetra_face_defs[MAX_NUM_FACES][2][2]

    int WEDGE_IND
    int WEDGE_NF
    np.uint8_t wedge_face_defs[MAX_NUM_FACES][2][2]


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixelize_cartesian(np.float64_t[:,:] buff,
                       np.float64_t[:] px,
                       np.float64_t[:] py,
                       np.float64_t[:] pdx,
                       np.float64_t[:] pdy,
                       np.float64_t[:] data,
                       bounds,
                       int antialias = 1,
                       period = None,
                       int check_period = 1,
                       np.float64_t line_width = 0.0):
    cdef np.float64_t x_min, x_max, y_min, y_max
    cdef np.float64_t period_x = 0.0, period_y = 0.0
    cdef np.float64_t width, height, px_dx, px_dy, ipx_dx, ipx_dy
    cdef np.float64_t ld_x, ld_y, cx, cy
    cdef int i, j, p, xi, yi
    cdef int lc, lr, rc, rr
    cdef np.float64_t lypx, rypx, lxpx, rxpx, overlap1, overlap2
    # These are the temp vars we get from the arrays
    cdef np.float64_t oxsp, oysp, xsp, ysp, dxsp, dysp, dsp
    # Some periodicity helpers
    cdef int xiter[2]
    cdef int yiter[2]
    cdef np.float64_t xiterv[2]
    cdef np.float64_t yiterv[2]
    if period is not None:
        period_x = period[0]
        period_y = period[1]
    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]
    width = x_max - x_min
    height = y_max - y_min
    px_dx = width / (<np.float64_t> buff.shape[1])
    px_dy = height / (<np.float64_t> buff.shape[0])
    ipx_dx = 1.0 / px_dx
    ipx_dy = 1.0 / px_dy
    if px.shape[0] != py.shape[0] or \
       px.shape[0] != pdx.shape[0] or \
       px.shape[0] != pdy.shape[0] or \
       px.shape[0] != data.shape[0]:
        raise YTPixelizeError("Arrays are not of correct shape.")
    xiter[0] = yiter[0] = 0
    xiterv[0] = yiterv[0] = 0.0
    # Here's a basic outline of what we're going to do here.  The xiter and
    # yiter variables govern whether or not we should check periodicity -- are
    # we both close enough to the edge that it would be important *and* are we
    # periodic?
    #
    # The other variables are all either pixel positions or data positions.
    # Pixel positions will vary regularly from the left edge of the window to
    # the right edge of the window; px_dx and px_dy are the dx (cell width, not
    # half-width).  ipx_dx and ipx_dy are the inverse, for quick math.
    #
    # The values in xsp, dxsp, x_min and their y counterparts, are the
    # data-space coordinates, and are related to the data fed in.  We make some
    # modifications for periodicity.
    #
    # Inside the finest loop, we compute the "left column" (lc) and "lower row"
    # (lr) and then iterate up to "right column" (rc) and "uppeR row" (rr),
    # depositing into them the data value.  Overlap computes the relative
    # overlap of a data value with a pixel.
    #
    # NOTE ON ROWS AND COLUMNS:
    #
    #   The way that images are plotting in matplotlib is somewhat different
    #   from what most might expect.  The first axis of the array plotted is
    #   what varies along the x axis.  So for instance, if you supply
    #   origin='lower' and plot the results of an mgrid operation, at a fixed
    #   'y' value you will see the results of that array held constant in the
    #   first dimension.  Here is some example code:
    #
    #   import matplotlib.pyplot as plt
    #   import numpy as np
    #   x, y = np.mgrid[0:1:100j,0:1:100j]
    #   plt.imshow(x, interpolation='nearest', origin='lower')
    #   plt.imshow(y, interpolation='nearest', origin='lower')
    #
    #   The values in the image:
    #       lower left:  arr[0,0]
    #       lower right: arr[0,-1]
    #       upper left:  arr[-1,0]
    #       upper right: arr[-1,-1]
    #
    #   So what we want here is to fill an array such that we fill:
    #       first axis : y_min .. y_max
    #       second axis: x_min .. x_max
    with nogil:
        for p in range(px.shape[0]):
            xiter[1] = yiter[1] = 999
            xiterv[1] = yiterv[1] = 0.0
            oxsp = px[p]
            oysp = py[p]
            dxsp = pdx[p]
            dysp = pdy[p]
            dsp = data[p]
            if check_period == 1:
                if (oxsp - dxsp < x_min):
                    xiter[1] = +1
                    xiterv[1] = period_x
                elif (oxsp + dxsp > x_max):
                    xiter[1] = -1
                    xiterv[1] = -period_x
                if (oysp - dysp < y_min):
                    yiter[1] = +1
                    yiterv[1] = period_y
                elif (oysp + dysp > y_max):
                    yiter[1] = -1
                    yiterv[1] = -period_y
            overlap1 = overlap2 = 1.0
            for xi in range(2):
                if xiter[xi] == 999: continue
                xsp = oxsp + xiterv[xi]
                if (xsp + dxsp < x_min) or (xsp - dxsp > x_max): continue
                for yi in range(2):
                    if yiter[yi] == 999: continue
                    ysp = oysp + yiterv[yi]
                    if (ysp + dysp < y_min) or (ysp - dysp > y_max): continue
                    lc = <int> fmax(((xsp-dxsp-x_min)*ipx_dx),0)
                    lr = <int> fmax(((ysp-dysp-y_min)*ipx_dy),0)
                    # NOTE: This is a different way of doing it than in the C
                    # routines.  In C, we were implicitly casting the
                    # initialization to int, but *not* the conditional, which
                    # was allowed an extra value:
                    #     for(j=lc;j<rc;j++)
                    # here, when assigning lc (double) to j (int) it got
                    # truncated, but no similar truncation was done in the
                    # comparison of j to rc (double).  So give ourselves a
                    # bonus row and bonus column here.
                    rc = <int> fmin(((xsp+dxsp-x_min)*ipx_dx + 1), buff.shape[1])
                    rr = <int> fmin(((ysp+dysp-y_min)*ipx_dy + 1), buff.shape[0])
                    # Note that we're iterating here over *y* in the i
                    # direction.  See the note above about this.
                    for i in range(lr, rr):
                        lypx = px_dy * i + y_min
                        rypx = px_dy * (i+1) + y_min
                        if antialias == 1:
                            overlap2 = ((fmin(rypx, ysp+dysp)
                                       - fmax(lypx, (ysp-dysp)))*ipx_dy)
                        if overlap2 < 0.0: continue
                        for j in range(lc, rc):
                            lxpx = px_dx * j + x_min
                            rxpx = px_dx * (j+1) + x_min
                            if line_width > 0:
                                # Here, we figure out if we're within
                                # line_width*px_dx of the cell edge
                                # Midpoint of x:
                                cx = (rxpx+lxpx)*0.5
                                ld_x = fmin(fabs(cx - (xsp+dxsp)),
                                            fabs(cx - (xsp-dxsp)))
                                ld_x *= ipx_dx
                                # Midpoint of y:
                                cy = (rypx+lypx)*0.5
                                ld_y = fmin(fabs(cy - (ysp+dysp)),
                                            fabs(cy - (ysp-dysp)))
                                ld_y *= ipx_dy
                                if ld_x <= line_width or ld_y <= line_width:
                                    buff[i,j] = 1.0
                            elif antialias == 1:
                                overlap1 = ((fmin(rxpx, xsp+dxsp)
                                           - fmax(lxpx, (xsp-dxsp)))*ipx_dx)
                                if overlap1 < 0.0: continue
                                # This next line is not commented out because
                                # it's an oddity; we actually want to skip
                                # depositing if the overlap is zero, and that's
                                # how it used to work when we were more
                                # conservative about the iteration indices.
                                # This will reduce artifacts if we ever move to
                                # compositing instead of replacing bitmaps.
                                if overlap1 * overlap2 < 1.e-6: continue
                                buff[i,j] += (dsp * overlap1) * overlap2
                            else:
                                buff[i,j] = dsp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixelize_cartesian_nodal(np.float64_t[:,:] buff,
                             np.float64_t[:] px,
                             np.float64_t[:] py,
                             np.float64_t[:] pz,
                             np.float64_t[:] pdx,
                             np.float64_t[:] pdy,
                             np.float64_t[:] pdz,
                             np.float64_t[:, :] data,
                             np.float64_t coord,
                             bounds,
                             int antialias = 1,
                             period = None,
                             int check_period = 1):
    cdef np.float64_t x_min, x_max, y_min, y_max
    cdef np.float64_t period_x = 0.0, period_y = 0.0
    cdef np.float64_t width, height, px_dx, px_dy, ipx_dx, ipx_dy
    cdef np.float64_t ld_x, ld_y, cx, cy, cz
    cdef int i, j, p, xi, yi
    cdef int lc, lr, rc, rr
    cdef np.float64_t lypx, rypx, lxpx, rxpx, overlap1, overlap2
    # These are the temp vars we get from the arrays
    cdef np.float64_t oxsp, oysp, ozsp
    cdef np.float64_t xsp, ysp, zsp
    cdef np.float64_t dxsp, dysp, dzsp
    # Some periodicity helpers
    cdef int xiter[2]
    cdef int yiter[2]
    cdef int ii, jj, kk, ind
    cdef np.float64_t xiterv[2]
    cdef np.float64_t yiterv[2]
    if period is not None:
        period_x = period[0]
        period_y = period[1]
    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]
    width = x_max - x_min
    height = y_max - y_min
    px_dx = width / (<np.float64_t> buff.shape[1])
    px_dy = height / (<np.float64_t> buff.shape[0])
    ipx_dx = 1.0 / px_dx
    ipx_dy = 1.0 / px_dy
    if px.shape[0] != py.shape[0] or \
       px.shape[0] != pz.shape[0] or \
       px.shape[0] != pdx.shape[0] or \
       px.shape[0] != pdy.shape[0] or \
       px.shape[0] != pdz.shape[0] or \
       px.shape[0] != data.shape[0]:
        raise YTPixelizeError("Arrays are not of correct shape.")
    xiter[0] = yiter[0] = 0
    xiterv[0] = yiterv[0] = 0.0
    # Here's a basic outline of what we're going to do here.  The xiter and
    # yiter variables govern whether or not we should check periodicity -- are
    # we both close enough to the edge that it would be important *and* are we
    # periodic?
    #
    # The other variables are all either pixel positions or data positions.
    # Pixel positions will vary regularly from the left edge of the window to
    # the right edge of the window; px_dx and px_dy are the dx (cell width, not
    # half-width).  ipx_dx and ipx_dy are the inverse, for quick math.
    #
    # The values in xsp, dxsp, x_min and their y counterparts, are the
    # data-space coordinates, and are related to the data fed in.  We make some
    # modifications for periodicity.
    #
    # Inside the finest loop, we compute the "left column" (lc) and "lower row"
    # (lr) and then iterate up to "right column" (rc) and "uppeR row" (rr),
    # depositing into them the data value.  Overlap computes the relative
    # overlap of a data value with a pixel.
    #
    # NOTE ON ROWS AND COLUMNS:
    #
    #   The way that images are plotting in matplotlib is somewhat different
    #   from what most might expect.  The first axis of the array plotted is
    #   what varies along the x axis.  So for instance, if you supply
    #   origin='lower' and plot the results of an mgrid operation, at a fixed
    #   'y' value you will see the results of that array held constant in the
    #   first dimension.  Here is some example code:
    #
    #   import matplotlib.pyplot as plt
    #   import numpy as np
    #   x, y = np.mgrid[0:1:100j,0:1:100j]
    #   plt.imshow(x, interpolation='nearest', origin='lower')
    #   plt.imshow(y, interpolation='nearest', origin='lower')
    #
    #   The values in the image:
    #       lower left:  arr[0,0]
    #       lower right: arr[0,-1]
    #       upper left:  arr[-1,0]
    #       upper right: arr[-1,-1]
    #
    #   So what we want here is to fill an array such that we fill:
    #       first axis : y_min .. y_max
    #       second axis: x_min .. x_max
    with nogil:
        for p in range(px.shape[0]):
            xiter[1] = yiter[1] = 999
            xiterv[1] = yiterv[1] = 0.0
            oxsp = px[p]
            oysp = py[p]
            ozsp = pz[p]
            dxsp = pdx[p]
            dysp = pdy[p]
            dzsp = pdz[p]
            if check_period == 1:
                if (oxsp - dxsp < x_min):
                    xiter[1] = +1
                    xiterv[1] = period_x
                elif (oxsp + dxsp > x_max):
                    xiter[1] = -1
                    xiterv[1] = -period_x
                if (oysp - dysp < y_min):
                    yiter[1] = +1
                    yiterv[1] = period_y
                elif (oysp + dysp > y_max):
                    yiter[1] = -1
                    yiterv[1] = -period_y
            overlap1 = overlap2 = 1.0
            zsp = ozsp
            for xi in range(2):
                if xiter[xi] == 999: continue
                xsp = oxsp + xiterv[xi]
                if (xsp + dxsp < x_min) or (xsp - dxsp > x_max): continue
                for yi in range(2):
                    if yiter[yi] == 999: continue
                    ysp = oysp + yiterv[yi]
                    if (ysp + dysp < y_min) or (ysp - dysp > y_max): continue
                    lc = <int> fmax(((xsp-dxsp-x_min)*ipx_dx),0)
                    lr = <int> fmax(((ysp-dysp-y_min)*ipx_dy),0)
                    # NOTE: This is a different way of doing it than in the C
                    # routines.  In C, we were implicitly casting the
                    # initialization to int, but *not* the conditional, which
                    # was allowed an extra value:
                    #     for(j=lc;j<rc;j++)
                    # here, when assigning lc (double) to j (int) it got
                    # truncated, but no similar truncation was done in the
                    # comparison of j to rc (double).  So give ourselves a
                    # bonus row and bonus column here.
                    rc = <int> fmin(((xsp+dxsp-x_min)*ipx_dx + 1), buff.shape[1])
                    rr = <int> fmin(((ysp+dysp-y_min)*ipx_dy + 1), buff.shape[0])
                    # Note that we're iterating here over *y* in the i
                    # direction.  See the note above about this.
                    for i in range(lr, rr):
                        lypx = px_dy * i + y_min
                        rypx = px_dy * (i+1) + y_min
                        for j in range(lc, rc):
                            lxpx = px_dx * j + x_min
                            rxpx = px_dx * (j+1) + x_min

                            cx = (rxpx+lxpx)*0.5
                            cy = (rypx+lypx)*0.5
                            cz = coord

                            ii = <int> (cx - xsp + dxsp)
                            jj = <int> (cy - ysp + dysp)
                            kk = <int> (cz - zsp + dzsp)

                            ind = 4*ii + 2*jj + kk

                            buff[i,j] = data[p, ind]


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixelize_off_axis_cartesian(
                       np.float64_t[:,:] buff,
                       np.float64_t[:] x,
                       np.float64_t[:] y,
                       np.float64_t[:] z,
                       np.float64_t[:] px,
                       np.float64_t[:] py,
                       np.float64_t[:] pdx,
                       np.float64_t[:] pdy,
                       np.float64_t[:] pdz,
                       np.float64_t[:] center,
                       np.float64_t[:,:] inv_mat,
                       np.int_t[:] indices,
                       np.float64_t[:] data,
                       bounds):
    cdef np.float64_t x_min, x_max, y_min, y_max
    cdef np.float64_t width, height, px_dx, px_dy, ipx_dx, ipx_dy, md
    cdef int i, j, p, ip
    cdef int lc, lr, rc, rr
    # These are the temp vars we get from the arrays
    cdef np.float64_t xsp, ysp, zsp, dxsp, dysp, dzsp, dsp
    cdef np.float64_t pxsp, pysp, cxpx, cypx, cx, cy, cz
    # Some periodicity helpers
    cdef np.ndarray[np.int64_t, ndim=2] mask
    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]
    width = x_max - x_min
    height = y_max - y_min
    px_dx = width / (<np.float64_t> buff.shape[1])
    px_dy = height / (<np.float64_t> buff.shape[0])
    ipx_dx = 1.0 / px_dx
    ipx_dy = 1.0 / px_dy
    if px.shape[0] != py.shape[0] or \
       px.shape[0] != pdx.shape[0] or \
       px.shape[0] != pdy.shape[0] or \
       px.shape[0] != pdz.shape[0] or \
       px.shape[0] != indices.shape[0] or \
       px.shape[0] != data.shape[0]:
        raise YTPixelizeError("Arrays are not of correct shape.")
    mask = np.zeros((buff.shape[0], buff.shape[1]), "int64")
    with nogil:
        for ip in range(indices.shape[0]):
            p = indices[ip]
            xsp = x[p]
            ysp = y[p]
            zsp = z[p]
            pxsp = px[p]
            pysp = py[p]
            dxsp = pdx[p]
            dysp = pdy[p]
            dzsp = pdz[p]
            dsp = data[p]
            # Any point we want to plot is at most this far from the center
            md = 2.0 * math.sqrt(dxsp*dxsp + dysp*dysp + dzsp*dzsp)
            if pxsp + md < x_min or \
               pxsp - md > x_max or \
               pysp + md < y_min or \
               pysp - md > y_max:
                continue
            lc = <int> fmax(((pxsp - md - x_min)*ipx_dx),0)
            lr = <int> fmax(((pysp - md - y_min)*ipx_dy),0)
            rc = <int> fmin(((pxsp + md - x_min)*ipx_dx + 1), buff.shape[1])
            rr = <int> fmin(((pysp + md - y_min)*ipx_dy + 1), buff.shape[0])
            for i in range(lr, rr):
                cypx = px_dy * (i + 0.5) + y_min
                for j in range(lc, rc):
                    cxpx = px_dx * (j + 0.5) + x_min
                    cx = inv_mat[0,0]*cxpx + inv_mat[0,1]*cypx + center[0]
                    cy = inv_mat[1,0]*cxpx + inv_mat[1,1]*cypx + center[1]
                    cz = inv_mat[2,0]*cxpx + inv_mat[2,1]*cypx + center[2]
                    if fabs(xsp - cx) * 0.99 > dxsp or \
                       fabs(ysp - cy) * 0.99 > dysp or \
                       fabs(zsp - cz) * 0.99 > dzsp:
                        continue
                    mask[i, j] += 1
                    buff[i, j] += dsp
    for i in range(buff.shape[0]):
        for j in range(buff.shape[1]):
            if mask[i,j] == 0: continue
            buff[i,j] /= mask[i,j]

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixelize_cylinder(np.float64_t[:,:] buff,
                      np.float64_t[:] radius,
                      np.float64_t[:] dradius,
                      np.float64_t[:] theta,
                      np.float64_t[:] dtheta,
                      np.float64_t[:] field,
                      extents):

    cdef np.float64_t x, y, dx, dy, r0, theta0
    cdef np.float64_t rmax, x0, y0, x1, y1
    cdef np.float64_t r_i, theta_i, dr_i, dtheta_i, dthetamin
    cdef np.float64_t costheta, sintheta
    cdef int i, pi, pj

    cdef int imax = np.asarray(radius).argmax()
    rmax = radius[imax] + dradius[imax]

    x0, x1, y0, y1 = extents
    dx = (x1 - x0) / buff.shape[1]
    dy = (y1 - y0) / buff.shape[0]
    cdef np.float64_t rbounds[2]
    cdef np.float64_t corners[8]
    # Find our min and max r
    corners[0] = x0*x0+y0*y0
    corners[1] = x1*x1+y0*y0
    corners[2] = x0*x0+y1*y1
    corners[3] = x1*x1+y1*y1
    corners[4] = x0*x0
    corners[5] = x1*x1
    corners[6] = y0*y0
    corners[7] = y1*y1
    rbounds[0] = rbounds[1] = corners[0]
    for i in range(8):
        rbounds[0] = fmin(rbounds[0], corners[i])
        rbounds[1] = fmax(rbounds[1], corners[i])
    rbounds[0] = math.sqrt(rbounds[0])
    rbounds[1] = math.sqrt(rbounds[1])
    # If we include the origin in either direction, we need to have radius of
    # zero as our lower bound.
    if x0 < 0 and x1 > 0:
        rbounds[0] = 0.0
    if y0 < 0 and y1 > 0:
        rbounds[0] = 0.0
    dthetamin = dx / rmax
    for i in range(radius.shape[0]):

        r0 = radius[i]
        theta0 = theta[i]
        dr_i = dradius[i]
        dtheta_i = dtheta[i]
        # Skip out early if we're offsides, for zoomed in plots
        if r0 + dr_i < rbounds[0] or r0 - dr_i > rbounds[1]:
            continue
        theta_i = theta0 - dtheta_i
        # Buffer of 0.5 here
        dthetamin = 0.5*dx/(r0 + dr_i)
        while theta_i < theta0 + dtheta_i:
            r_i = r0 - dr_i
            costheta = math.cos(theta_i)
            sintheta = math.sin(theta_i)
            while r_i < r0 + dr_i:
                if rmax <= r_i:
                    r_i += 0.5*dx
                    continue
                y = r_i * costheta
                x = r_i * sintheta
                pi = <int>((x - x0)/dx)
                pj = <int>((y - y0)/dy)
                if pi >= 0 and pi < buff.shape[0] and \
                   pj >= 0 and pj < buff.shape[1]:
                    buff[pi, pj] = field[i]
                r_i += 0.5*dx
            theta_i += dthetamin

cdef int aitoff_thetaphi_to_xy(np.float64_t theta, np.float64_t phi,
                               np.float64_t *x, np.float64_t *y) except -1:
    cdef np.float64_t z = math.sqrt(1 + math.cos(phi) * math.cos(theta / 2.0))
    x[0] = math.cos(phi) * math.sin(theta / 2.0) / z
    y[0] = math.sin(phi) / z
    return 0

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixelize_aitoff(np.float64_t[:] theta,
                    np.float64_t[:] dtheta,
                    np.float64_t[:] phi,
                    np.float64_t[:] dphi,
                    buff_size,
                    np.float64_t[:] field,
                    extents, input_img = None,
                    np.float64_t theta_offset = 0.0,
                    np.float64_t phi_offset = 0.0):
    # http://paulbourke.net/geometry/transformationprojection/
    # (theta) longitude is -pi to pi
    # (phi) latitude is -pi/2 to pi/2
    # z^2 = 1 + cos(latitude) cos(longitude/2)
    # x = cos(latitude) sin(longitude/2) / z
    # y = sin(latitude) / z
    cdef np.ndarray[np.float64_t, ndim=2] img
    cdef int i, j, nf, fi
    cdef np.float64_t x, y, z, zb
    cdef np.float64_t dx, dy
    cdef np.float64_t theta0, phi0, theta_p, dtheta_p, phi_p, dphi_p
    cdef np.float64_t PI = np.pi
    cdef np.float64_t s2 = math.sqrt(2.0)
    cdef np.float64_t xmax, ymax, xmin, ymin
    nf = field.shape[0]

    if input_img is None:
        img = np.zeros((buff_size[0], buff_size[1]))
        img[:] = np.nan
    else:
        img = input_img
    # Okay, here's our strategy.  We compute the bounds in x and y, which will
    # be a rectangle, and then for each x, y position we check to see if it's
    # within our theta.  This will cost *more* computations of the
    # (x,y)->(theta,phi) calculation, but because we no longer have to search
    # through the theta, phi arrays, it should be faster.
    dx = 2.0 / (img.shape[0] - 1)
    dy = 2.0 / (img.shape[1] - 1)
    x = y = 0
    for fi in range(nf):
        theta_p = (theta[fi] + theta_offset) - PI
        dtheta_p = dtheta[fi]
        phi_p = (phi[fi] + phi_offset) - PI/2.0
        dphi_p = dphi[fi]
        # Four transformations
        aitoff_thetaphi_to_xy(theta_p - dtheta_p, phi_p - dphi_p, &x, &y)
        xmin = x
        xmax = x
        ymin = y
        ymax = y
        aitoff_thetaphi_to_xy(theta_p - dtheta_p, phi_p + dphi_p, &x, &y)
        xmin = fmin(xmin, x)
        xmax = fmax(xmax, x)
        ymin = fmin(ymin, y)
        ymax = fmax(ymax, y)
        aitoff_thetaphi_to_xy(theta_p + dtheta_p, phi_p - dphi_p, &x, &y)
        xmin = fmin(xmin, x)
        xmax = fmax(xmax, x)
        ymin = fmin(ymin, y)
        ymax = fmax(ymax, y)
        aitoff_thetaphi_to_xy(theta_p + dtheta_p, phi_p + dphi_p, &x, &y)
        xmin = fmin(xmin, x)
        xmax = fmax(xmax, x)
        ymin = fmin(ymin, y)
        ymax = fmax(ymax, y)
        # Now we have the (projected rectangular) bounds.
        xmin = (xmin + 1) # Get this into normalized image coords
        xmax = (xmax + 1) # Get this into normalized image coords
        ymin = (ymin + 1) # Get this into normalized image coords
        ymax = (ymax + 1) # Get this into normalized image coords
        x0 = <int> (xmin / dx)
        x1 = <int> (xmax / dx) + 1
        y0 = <int> (ymin / dy)
        y1 = <int> (ymax / dy) + 1
        for i in range(x0, x1):
            x = (-1.0 + i*dx)*s2*2.0
            for j in range(y0, y1):
                y = (-1.0 + j * dy)*s2
                zb = (x*x/8.0 + y*y/2.0 - 1.0)
                if zb > 0: continue
                z = (1.0 - (x * 0.25) * (x * 0.25) - (y * 0.5) * (y * 0.5))
                z = math.sqrt(z)
                # Longitude
                theta0 = 2.0*math.atan(z*x/(2.0 * (2.0*z*z-1.0)))
                # Latitude
                # We shift it into co-latitude
                phi0 = math.asin(z*y)
                # Now we just need to figure out which pixel contributes.
                # We do not have a fast search.
                if not (theta_p - dtheta_p <= theta0 <= theta_p + dtheta_p):
                    continue
                if not (phi_p - dphi_p <= phi0 <= phi_p + dphi_p):
                    continue
                img[i, j] = field[fi]
    return img


# This function accepts a set of vertices (for a polyhedron) that are
# assumed to be in order for bottom, then top, in the same clockwise or
# counterclockwise direction (i.e., like points 1-8 in Figure 4 of the ExodusII
# manual).  It will then either *match* or *fill* the results.  If it is
# matching, it will early terminate with a 0 or final-terminate with a 1 if the
# results match.  Otherwise, it will fill the signs with -1's and 1's to show
# the sign of the dot product of the point with the cross product of the face.
cdef int check_face_dot(int nvertices,
                        np.float64_t point[3],
                        np.float64_t **vertices,
                        np.int8_t *signs,
                        int match):
    # Because of how we are doing this, we do not *care* what the signs are or
    # how the faces are ordered, we only care if they match between the point
    # and the centroid.
    # So, let's compute these vectors.  See above where these are written out
    # for ease of use.
    cdef np.float64_t vec1[3]
    cdef np.float64_t vec2[3]
    cdef np.float64_t cp_vec[3]
    cdef np.float64_t npoint[3]
    cdef np.float64_t dp
    cdef np.uint8_t faces[MAX_NUM_FACES][2][2]
    cdef np.uint8_t nf
    if nvertices == 4:
        faces = tetra_face_defs
        nf = TETRA_NF
    elif nvertices == 6:
        faces = wedge_face_defs
        nf = WEDGE_NF
    elif nvertices == 8:
        faces = hex_face_defs
        nf = HEX_NF
    else:
        return -1
    cdef int i, j, n, vi1a, vi1b, vi2a, vi2b

    for n in range(nf):
        vi1a = faces[n][0][0]
        vi1b = faces[n][0][1]
        vi2a = faces[n][1][0]
        vi2b = faces[n][1][1]
        # Shared vertex is vi1a and vi2a
        subtract(vertices[vi1b], vertices[vi1a], vec1)
        subtract(vertices[vi2b], vertices[vi2a], vec2)
        subtract(point, vertices[vi1b], npoint)
        cross(vec1, vec2, cp_vec)
        dp = dot(cp_vec, npoint)
        if match == 0:
            if dp < 0:
                signs[n] = -1
            else:
                signs[n] = 1
        else:
            if dp <= 0 and signs[n] < 0:
                continue
            elif dp >= 0 and signs[n] > 0:
                continue
            else: # mismatch!
                return 0
    return 1


def pixelize_element_mesh(np.ndarray[np.float64_t, ndim=2] coords,
                          np.ndarray[np.int64_t, ndim=2] conn,
                          buff_size,
                          np.ndarray[np.float64_t, ndim=2] field,
                          extents,
                          int index_offset = 0):
    cdef np.ndarray[np.float64_t, ndim=3] img
    img = np.zeros(buff_size, dtype="float64")
    img[:] = np.nan
    # Two steps:
    #  1. Is image point within the mesh bounding box?
    #  2. Is image point within the mesh element?
    # Second is more intensive.  It will convert the element vertices to the
    # mapped coordinate system, and check whether the result in in-bounds or not
    # Note that we have to have a pseudo-3D pixel buffer.  One dimension will
    # always be 1.
    cdef np.float64_t pLE[3]
    cdef np.float64_t pRE[3]
    cdef np.float64_t LE[3]
    cdef np.float64_t RE[3]
    cdef int use
    cdef np.int64_t n, i, pi, pj, pk, ci, cj
    cdef np.int64_t pstart[3]
    cdef np.int64_t pend[3]
    cdef np.float64_t ppoint[3]
    cdef np.float64_t idds[3]
    cdef np.float64_t dds[3]
    cdef np.float64_t *vertices
    cdef np.float64_t *field_vals
    cdef int nvertices = conn.shape[1]
    cdef int ndim = coords.shape[1]
    cdef int num_field_vals = field.shape[1]
    cdef double[4] mapped_coord
    cdef ElementSampler sampler

    # Pick the right sampler and allocate storage for the mapped coordinate
    if ndim == 3 and nvertices == 4:
        sampler = P1Sampler3D()
    elif ndim == 3 and nvertices == 6:
        sampler = W1Sampler3D()
    elif ndim == 3 and nvertices == 8:
        sampler = Q1Sampler3D()
    elif ndim == 3 and nvertices == 20:
        sampler = S2Sampler3D()
    elif ndim == 2 and nvertices == 3:
        sampler = P1Sampler2D()
    elif ndim == 1 and nvertices == 2:
        sampler = P1Sampler1D()
    elif ndim == 2 and nvertices == 4:
        sampler = Q1Sampler2D()
    elif ndim == 2 and nvertices == 9:
        sampler = Q2Sampler2D()
    elif ndim == 2 and nvertices == 6:
        sampler = T2Sampler2D()
    elif ndim == 3 and nvertices == 10:
        sampler = Tet2Sampler3D()
    else:
        raise YTElementTypeNotRecognized(ndim, nvertices)

    # if we are in 2D land, the 1 cell thick dimension had better be 'z'
    if ndim == 2:
        if buff_size[2] != 1:
            raise RuntimeError("Slices of 2D datasets must be "
                               "perpendicular to the 'z' direction.")

    # allocate temporary storage
    vertices = <np.float64_t *> malloc(ndim * sizeof(np.float64_t) * nvertices)
    field_vals = <np.float64_t *> malloc(sizeof(np.float64_t) * num_field_vals)

    # fill the image bounds and pixel size information here
    for i in range(ndim):
        pLE[i] = extents[i][0]
        pRE[i] = extents[i][1]
        dds[i] = (pRE[i] - pLE[i])/buff_size[i]
        if dds[i] == 0.0:
            idds[i] = 0.0
        else:
            idds[i] = 1.0 / dds[i]

    with cython.boundscheck(False):
        for ci in range(conn.shape[0]):

            # Fill the vertices
            LE[0] = LE[1] = LE[2] = 1e60
            RE[0] = RE[1] = RE[2] = -1e60

            for n in range(num_field_vals):
                field_vals[n] = field[ci, n]

            for n in range(nvertices):
                cj = conn[ci, n] - index_offset
                for i in range(ndim):
                    vertices[ndim*n + i] = coords[cj, i]
                    LE[i] = fmin(LE[i], vertices[ndim*n+i])
                    RE[i] = fmax(RE[i], vertices[ndim*n+i])

            use = 1
            for i in range(ndim):
                if RE[i] < pLE[i] or LE[i] >= pRE[i]:
                    use = 0
                    break
                pstart[i] = i64max(<np.int64_t> ((LE[i] - pLE[i])*idds[i]) - 1, 0)
                pend[i] = i64min(<np.int64_t> ((RE[i] - pLE[i])*idds[i]) + 1, img.shape[i]-1)

            # override for the low-dimensional case
            if ndim < 3:
                pstart[2] = 0
                pend[2] = 0
            if ndim < 2:
                pstart[1] = 0
                pend[1] = 0

            if use == 0:
                continue

            # Now our bounding box intersects, so we get the extents of our pixel
            # region which overlaps with the bounding box, and we'll check each
            # pixel in there.
            for pi in range(pstart[0], pend[0] + 1):
                ppoint[0] = (pi + 0.5) * dds[0] + pLE[0]
                for pj in range(pstart[1], pend[1] + 1):
                    ppoint[1] = (pj + 0.5) * dds[1] + pLE[1]
                    for pk in range(pstart[2], pend[2] + 1):
                        ppoint[2] = (pk + 0.5) * dds[2] + pLE[2]
                        # Now we just need to figure out if our ppoint is within
                        # our set of vertices.
                        sampler.map_real_to_unit(mapped_coord, vertices, ppoint)
                        if not sampler.check_inside(mapped_coord):
                            continue
                        if (num_field_vals == 1):
                            img[pi, pj, pk] = field_vals[0]
                        else:
                            img[pi, pj, pk] = sampler.sample_at_unit_point(mapped_coord,
                                                                           field_vals)
    free(vertices)
    free(field_vals)
    return img

# used as a cache to avoid repeatedly creating
# instances of SPHKernelInterpolationTable
kernel_tables = {}

cdef class SPHKernelInterpolationTable:
    cdef public object kernel_name
    cdef kernel_func kernel
    cdef np.float64_t[::1] table
    cdef np.float64_t[::1] q2_vals
    cdef np.float64_t q2_range, iq2_range

    def __init__(self, kernel_name):
        self.kernel_name = kernel_name
        self.kernel = get_kernel_func(kernel_name)
        self.populate_table()

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.float64_t integrate_q2(self, np.float64_t q2) nogil:
        # See equation 30 of the SPLASH paper
        cdef int i
        # Our bounds are -sqrt(R*R - q2) and sqrt(R*R-q2)
        # And our R is always 1; note that our smoothing kernel functions
        # expect it to run from 0 .. 1, so we multiply the integrand by 2
        cdef int N = 200
        cdef np.float64_t qz
        cdef np.float64_t R = 1
        cdef np.float64_t R0 = -math.sqrt(R*R-q2)
        cdef np.float64_t R1 = math.sqrt(R*R-q2)
        cdef np.float64_t dR = (R1-R0)/N
        # Set to our bounds
        cdef np.float64_t integral = 0.0
        integral += self.kernel(math.sqrt(R0*R0 + q2))
        integral += self.kernel(math.sqrt(R1*R1 + q2))
        # We're going to manually conduct a trapezoidal integration
        for i in range(1, N):
            qz = R0 + i * dR
            integral += 2.0*self.kernel(math.sqrt(qz*qz + q2))
        integral *= (R1-R0)/(2*N)
        return integral

    def populate_table(self):
        cdef int i
        self.table = cvarray(format="d", shape=(TABLE_NVALS,),
                             itemsize=sizeof(np.float64_t))
        self.q2_vals = cvarray(format="d", shape=(TABLE_NVALS,),
                             itemsize=sizeof(np.float64_t))
        # We run from 0 to 1 here over R
        for i in range(TABLE_NVALS):
            self.q2_vals[i] = i * 1.0/(TABLE_NVALS-1)
            self.table[i] = self.integrate_q2(self.q2_vals[i])

        self.q2_range = self.q2_vals[TABLE_NVALS-1] - self.q2_vals[0]
        self.iq2_range = (TABLE_NVALS-1)/self.q2_range

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline np.float64_t interpolate(self, np.float64_t q2) nogil:
        cdef int index
        cdef np.float64_t F_interpolate
        index = <int>((q2 - self.q2_vals[0])*(self.iq2_range))
        if index >= TABLE_NVALS:
            return 0.0
        F_interpolate = self.table[index] + (
                (self.table[index+1] - self.table[index])
               *(q2 - self.q2_vals[index])*self.iq2_range)
        return F_interpolate

    def interpolate_array(self, np.float64_t[:] q2_vals):
        cdef np.float64_t[:] ret = np.empty(q2_vals.shape[0])
        cdef int i
        for i in range(q2_vals.shape[0]):
            ret[i] = self.interpolate(q2_vals[i])
        return np.array(ret)

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pixelize_sph_kernel_projection(
        np.float64_t[:, :] buff,
        np.float64_t[:] posx,
        np.float64_t[:] posy,
        np.float64_t[:] hsml,
        np.float64_t[:] pmass,
        np.float64_t[:] pdens,
        np.float64_t[:] quantity_to_smooth,
        bounds,
        kernel_name="cubic",
        weight_field=None,
        int check_period=1,
        period=None):

    cdef np.intp_t xsize, ysize
    cdef np.float64_t x_min, x_max, y_min, y_max, prefactor_j
    cdef np.int64_t xi, yi, x0, x1, y0, y1, xxi, yyi
    cdef np.float64_t q_ij2, posx_diff, posy_diff, ih_j2
    cdef np.float64_t x, y, dx, dy, idx, idy, h_j2, px, py
    cdef np.float64_t period_x = 0, period_y = 0
    cdef int index, i, j, ii, jj
    cdef np.float64_t[:] _weight_field
    cdef int * xiter
    cdef int * yiter
    cdef np.float64_t * xiterv
    cdef np.float64_t * yiterv
    cdef np.float64_t * local_buf

    if weight_field is not None:
        _weight_field = weight_field

    if period is not None:
        period_x = period[0]
        period_y = period[1]

    # we find the x and y range over which we have pixels and we find how many
    # pixels we have in each dimension
    xsize, ysize = buff.shape[0], buff.shape[1]
    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]

    dx = (x_max - x_min) / xsize
    dy = (y_max - y_min) / ysize

    idx = 1.0/dx
    idy = 1.0/dy

    if kernel_name not in kernel_tables:
        kernel_tables[kernel_name] = SPHKernelInterpolationTable(kernel_name)
    cdef SPHKernelInterpolationTable itab = kernel_tables[kernel_name]

    with nogil, parallel():
        # loop through every particle
        # NOTE: this loop can be quite time consuming. However it is easily
        # parallelizable in multiple ways, such as:
        #   1) use multiple cores to process individual particles (the outer loop)
        #   2) use multiple cores to process individual pixels for a given particle
        #      (the inner loops)
        # Depending on the ratio of particles' "sphere of influence" (a.k.a. the smoothing
        # length) to the physical width of the pixels, different parallelization
        # strategies may yield different speed-ups. Strategy #1 works better in the case
        # of lots of itty bitty particles. Strategy #2 works well when we have a
        # not-very-large-number of reasonably large-compared-to-pixels particles. We
        # currently employ #1 as its workload is more even and consistent, even though it
        # comes with a price of an additional, per thread memory for storing the
        # intermediate results.

        local_buff = <np.float64_t *> malloc(sizeof(np.float64_t) * xsize * ysize)
        xiterv = <np.float64_t *> malloc(sizeof(np.float64_t) * 2)
        yiterv = <np.float64_t *> malloc(sizeof(np.float64_t) * 2)
        xiter = <int *> malloc(sizeof(int) * 2)
        yiter = <int *> malloc(sizeof(int) * 2)
        xiter[0] = yiter[0] = 0
        xiterv[0] = yiterv[0] = 0.0
        for i in range(xsize * ysize):
            local_buff[i] = 0.0

        for j in prange(0, posx.shape[0], schedule="dynamic"):
            if j % 100000 == 0:
                with gil:
                    PyErr_CheckSignals()

            xiter[1] = yiter[1] = 999

            if check_period == 1:
                if posx[j] - hsml[j] < x_min:
                    xiter[1] = +1
                    xiterv[1] = period_x
                elif posx[j] + hsml[j] > x_max:
                    xiter[1] = -1
                    xiterv[1] = -period_x
                if posy[j] - hsml[j] < y_min:
                    yiter[1] = +1
                    yiterv[1] = period_y
                elif posy[j] + hsml[j] > y_max:
                    yiter[1] = -1
                    yiterv[1] = -period_y

            # we set the smoothing length squared with lower limit of the pixel
            h_j2 = fmax(hsml[j]*hsml[j], dx*dy)
            ih_j2 = 1.0/h_j2

            prefactor_j = pmass[j] / pdens[j] / hsml[j]**2 * quantity_to_smooth[j]
            if weight_field is not None:
                prefactor_j *= _weight_field[j]

            for ii in range(2):
                if xiter[ii] == 999: continue
                px = posx[j] + xiterv[ii]
                if (px + hsml[j] < x_min) or (px - hsml[j] > x_max): continue
                for jj in range(2):
                    if yiter[jj] == 999: continue
                    py = posy[j] + yiterv[jj]
                    if (py + hsml[j] < y_min) or (py - hsml[j] > y_max): continue

                    # here we find the pixels which this particle contributes to
                    x0 = <np.int64_t> ((px - hsml[j] - x_min)*idx)
                    x1 = <np.int64_t> ((px + hsml[j] - x_min)*idx)
                    x0 = iclip(x0-1, 0, xsize)
                    x1 = iclip(x1+1, 0, xsize)

                    y0 = <np.int64_t> ((py - hsml[j] - y_min)*idy)
                    y1 = <np.int64_t> ((py + hsml[j] - y_min)*idy)
                    y0 = iclip(y0-1, 0, ysize)
                    y1 = iclip(y1+1, 0, ysize)

                    # found pixels we deposit on, loop through those pixels
                    for xi in range(x0, x1):
                        # we use the centre of the pixel to calculate contribution
                        x = (xi + 0.5) * dx + x_min

                        posx_diff = px - x
                        posx_diff = posx_diff * posx_diff

                        if posx_diff > h_j2: continue

                        for yi in range(y0, y1):
                            y = (yi + 0.5) * dy + y_min

                            posy_diff = py - y
                            posy_diff = posy_diff * posy_diff
                            if posy_diff > h_j2: continue

                            q_ij2 = (posx_diff + posy_diff) * ih_j2
                            if q_ij2 >= 1:
                                continue

                            # see equation 32 of the SPLASH paper
                            # now we just use the kernel projection
                            local_buff[xi + yi*xsize] +=  prefactor_j * itab.interpolate(q_ij2)

        with gil:
            for xxi in range(xsize):
                for yyi in range(ysize):
                    buff[xxi, yyi] += local_buff[xxi + yyi*xsize]
        free(local_buff)
        free(xiterv)
        free(yiterv)
        free(xiter)
        free(yiter)

@cython.boundscheck(False)
@cython.wraparound(False)
def interpolate_sph_positions_gather(np.float64_t[:] buff,
        np.float64_t[:, ::1] tree_positions, np.float64_t[:, ::1] field_positions,
        np.float64_t[:] hsml, np.float64_t[:] pmass, np.float64_t[:] pdens,
        np.float64_t[:] quantity_to_smooth, PyKDTree kdtree,
        int use_normalization=1, kernel_name="cubic", pbar=None,
        int num_neigh=32):

    """
    This function takes in arbitrary positions, field_positions, at which to
    perform a nearest neighbor search and perform SPH interpolation.

    The results are stored in the buffer, buff, which is in the same order as
    the field_positions are put in.
    """

    cdef np.float64_t q_ij, h_j2, ih_j2, prefactor_j, smoothed_quantity_j
    cdef np.float64_t * pos_ptr
    cdef int i, particle, index
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.float64_t[:] buff_den
    cdef KDTree * ctree = kdtree._tree

    # Which dimensions shall we use for spatial distances?
    cdef axes_range axes
    set_axes_range(&axes, -1)

    # Only allocate memory if we are using normalization
    if use_normalization:
        buff_den = np.zeros(buff.shape[0], dtype="float64")

    kernel_func = get_kernel_func(kernel_name)

    # Loop through all the positions we want to interpolate the SPH field onto
    with nogil:
        for i in range(0, buff.shape[0]):
            queue.size = 0

            # Update the current position
            pos_ptr = &field_positions[i, 0]

            # Use the KDTree to find the nearest neighbors
            find_neighbors(pos_ptr, tree_positions, queue, ctree, -1, &axes)

            # Set the smoothing length squared to the square of the distance
            # of the furthest nearest neighbor
            h_j2 = queue.heap[0]
            ih_j2 = 1.0/h_j2

            # Loop through each nearest neighbor and add contribution to the
            # buffer
            for index in range(queue.max_elements):
                particle = queue.pids[index]

                # Calculate contribution of this particle
                prefactor_j = (pmass[particle] / pdens[particle] /
                               hsml[particle]**3)
                q_ij = math.sqrt(queue.heap[index]*ih_j2)
                smoothed_quantity_j = (prefactor_j *
                                       quantity_to_smooth[particle] *
                                       kernel_func(q_ij))

                # See equations 6, 9, and 11 of the SPLASH paper
                buff[i] += smoothed_quantity_j

                if use_normalization:
                    buff_den[i] += prefactor_j * kernel_func(q_ij)

    if use_normalization:
        normalization_1d_utility(buff, buff_den)

@cython.boundscheck(False)
@cython.wraparound(False)
def interpolate_sph_grid_gather(np.float64_t[:, :, :] buff,
        np.float64_t[:, ::1] tree_positions, np.float64_t[:] bounds,
        np.float64_t[:] hsml, np.float64_t[:] pmass, np.float64_t[:] pdens,
        np.float64_t[:] quantity_to_smooth, PyKDTree kdtree,
        int use_normalization=1, kernel_name="cubic", pbar=None,
        int num_neigh=32):
    """
    This function takes in the bounds and number of cells in a grid (well,
    actually we implicity calculate this from the size of buff). Then we can
    perform nearest neighbor search and SPH interpolation at the centre of each
    cell in the grid.
    """

    cdef np.float64_t q_ij, h_j2, ih_j2, prefactor_j, smoothed_quantity_j
    cdef np.float64_t dx, dy, dz
    cdef np.float64_t[::1] pos = np.zeros(3, dtype="float64")
    cdef np.float64_t * pos_ptr = &pos[0]
    cdef int i, j, k, particle, index
    cdef BoundedPriorityQueue queue = BoundedPriorityQueue(num_neigh, True)
    cdef np.float64_t[:, :, :] buff_den
    cdef KDTree * ctree = kdtree._tree
    cdef int prog

    # Which dimensions shall we use for spatial distances?
    cdef axes_range axes
    set_axes_range(&axes, -1)

    # Only allocate memory if we are using normalization
    if use_normalization:
        buff_den = np.zeros([buff.shape[0], buff.shape[1],
                             buff.shape[2]], dtype="float64")

    kernel_func = get_kernel_func(kernel_name)
    dx = (bounds[1] - bounds[0]) / buff.shape[0]
    dy = (bounds[3] - bounds[2]) / buff.shape[1]
    dz = (bounds[5] - bounds[4]) / buff.shape[2]

    # Loop through all the positions we want to interpolate the SPH field onto
    pbar = get_pbar(title="Interpolating (gather) SPH field",
                    maxval=(buff.shape[0]*buff.shape[1]*buff.shape[2] //
                            10000)*10000)

    prog = 0
    with nogil:
        for i in range(0, buff.shape[0]):
            for j in range(0, buff.shape[1]):
                for k in range(0, buff.shape[2]):
                    prog += 1
                    if prog % 10000 == 0:
                        with gil:
                            PyErr_CheckSignals()
                            pbar.update(prog)

                    queue.size = 0

                    # Update the current position
                    pos[0] = bounds[0] + (i + 0.5) * dx
                    pos[1] = bounds[2] + (j + 0.5) * dy
                    pos[2] = bounds[4] + (k + 0.5) * dz

                    # Use the KDTree to find the nearest neighbors
                    find_neighbors(pos_ptr, tree_positions, queue, ctree, -1, &axes)

                    # Set the smoothing length squared to the square of the distance
                    # of the furthest nearest neighbor
                    h_j2 = queue.heap[0]
                    ih_j2 = 1.0/h_j2

                    # Loop through each nearest neighbor and add contribution to the
                    # buffer
                    for index in range(queue.max_elements):
                        particle = queue.pids[index]

                        # Calculate contribution of this particle
                        prefactor_j = (pmass[particle] / pdens[particle] /
                                       hsml[particle]**3)
                        q_ij = math.sqrt(queue.heap[index]*ih_j2)
                        smoothed_quantity_j = (prefactor_j *
                                               quantity_to_smooth[particle] *
                                               kernel_func(q_ij))

                        # See equations 6, 9, and 11 of the SPLASH paper
                        buff[i, j, k] += smoothed_quantity_j

                        if use_normalization:
                            buff_den[i, j, k] += prefactor_j * kernel_func(q_ij)

    if use_normalization:
        normalization_3d_utility(buff, buff_den)

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pixelize_sph_kernel_slice(
        np.float64_t[:, :] buff,
        np.float64_t[:] posx, np.float64_t[:] posy,
        np.float64_t[:] hsml, np.float64_t[:] pmass,
        np.float64_t[:] pdens,
        np.float64_t[:] quantity_to_smooth,
        bounds, kernel_name="cubic",
        int check_period=1,
        period=None):

    # similar method to pixelize_sph_kernel_projection
    cdef np.intp_t xsize, ysize
    cdef np.float64_t x_min, x_max, y_min, y_max, prefactor_j
    cdef np.int64_t xi, yi, x0, x1, y0, y1, xxi, yyi
    cdef np.float64_t q_ij, posx_diff, posy_diff, ih_j
    cdef np.float64_t x, y, dx, dy, idx, idy, h_j2, h_j, px, py
    cdef int index, i, j, ii, jj
    cdef np.float64_t period_x = 0, period_y = 0
    cdef int * xiter
    cdef int * yiter
    cdef np.float64_t * xiterv
    cdef np.float64_t * yiterv
    cdef np.float64_t * local_buf

    if period is not None:
        period_x = period[0]
        period_y = period[1]

    xsize, ysize = buff.shape[0], buff.shape[1]

    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]

    dx = (x_max - x_min) / xsize
    dy = (y_max - y_min) / ysize
    idx = 1.0/dx
    idy = 1.0/dy

    kernel_func = get_kernel_func(kernel_name)

    with nogil, parallel():
        # NOTE see note in pixelize_sph_kernel_projection
        local_buff = <np.float64_t *> malloc(sizeof(np.float64_t) * xsize * ysize)
        xiterv = <np.float64_t *> malloc(sizeof(np.float64_t) * 2)
        yiterv = <np.float64_t *> malloc(sizeof(np.float64_t) * 2)
        xiter = <int *> malloc(sizeof(int) * 2)
        yiter = <int *> malloc(sizeof(int) * 2)
        xiter[0] = yiter[0] = 0
        xiterv[0] = yiterv[0] = 0.0
        for i in range(xsize * ysize):
            local_buff[i] = 0.0

        for j in prange(0, posx.shape[0], schedule="dynamic"):
            if j % 100000 == 0:
                with gil:
                    PyErr_CheckSignals()

            xiter[1] = yiter[1] = 999

            if check_period == 1:
                if posx[j] - hsml[j] < x_min:
                    xiter[1] = +1
                    xiterv[1] = period_x
                elif posx[j] + hsml[j] > x_max:
                    xiter[1] = -1
                    xiterv[1] = -period_x
                if posy[j] - hsml[j] < y_min:
                    yiter[1] = +1
                    yiterv[1] = period_y
                elif posy[j] + hsml[j] > y_max:
                    yiter[1] = -1
                    yiterv[1] = -period_y

            h_j2 = fmax(hsml[j]*hsml[j], dx*dy)
            h_j = math.sqrt(h_j2)
            ih_j = 1.0/h_j

            prefactor_j = pmass[j] / pdens[j] / hsml[j]**3
            prefactor_j *= quantity_to_smooth[j]

            for ii in range(2):
                if xiter[ii] == 999: continue
                px = posx[j] + xiterv[ii]
                if (px + hsml[j] < x_min) or (px - hsml[j] > x_max): continue
                for jj in range(2):
                    if yiter[jj] == 999: continue
                    py = posy[j] + yiterv[jj]
                    if (py + hsml[j] < y_min) or (py - hsml[j] > y_max): continue

                    x0 = <np.int64_t> ( (px - hsml[j] - x_min) * idx)
                    x1 = <np.int64_t> ( (px + hsml[j] - x_min) * idx)
                    x0 = iclip(x0-1, 0, xsize)
                    x1 = iclip(x1+1, 0, xsize)

                    y0 = <np.int64_t> ( (py - hsml[j] - y_min) * idy)
                    y1 = <np.int64_t> ( (py + hsml[j] - y_min) * idy)
                    y0 = iclip(y0-1, 0, ysize)
                    y1 = iclip(y1+1, 0, ysize)

                    # Now we know which pixels to deposit onto for this particle,
                    # so loop over them and add this particle's contribution
                    for xi in range(x0, x1):
                        x = (xi + 0.5) * dx + x_min

                        posx_diff = px - x
                        posx_diff = posx_diff * posx_diff
                        if posx_diff > h_j2:
                            continue

                        for yi in range(y0, y1):
                            y = (yi + 0.5) * dy + y_min

                            posy_diff = py - y
                            posy_diff = posy_diff * posy_diff
                            if posy_diff > h_j2:
                                continue

                            # see equation 4 of the SPLASH paper
                            q_ij = math.sqrt(posx_diff + posy_diff) * ih_j
                            if q_ij >= 1:
                                continue

                            # see equations 6, 9, and 11 of the SPLASH paper
                            local_buff[xi + yi*xsize] += prefactor_j * kernel_func(q_ij)

        with gil:
            for xxi in range(xsize):
                for yyi in range(ysize):
                    buff[xxi, yyi] += local_buff[xxi + yyi*xsize]
        free(local_buff)
        free(xiterv)
        free(yiterv)
        free(xiter)
        free(yiter)

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pixelize_sph_kernel_arbitrary_grid(np.float64_t[:, :, :] buff,
        np.float64_t[:] posx, np.float64_t[:] posy, np.float64_t[:] posz,
        np.float64_t[:] hsml, np.float64_t[:] pmass,
        np.float64_t[:] pdens,
        np.float64_t[:] quantity_to_smooth,
        bounds, pbar=None, kernel_name="cubic",
        int check_period=1, period=None):

    cdef np.intp_t xsize, ysize, zsize
    cdef np.float64_t x_min, x_max, y_min, y_max, z_min, z_max, prefactor_j
    cdef np.int64_t xi, yi, zi, x0, x1, y0, y1, z0, z1
    cdef np.float64_t q_ij, posx_diff, posy_diff, posz_diff, px, py, pz
    cdef np.float64_t x, y, z, dx, dy, dz, idx, idy, idz, h_j3, h_j2, h_j, ih_j
    cdef int index, i, j, k, ii, jj, kk
    cdef np.float64_t period_x = 0, period_y = 0, period_z = 0

    cdef int xiter[2]
    cdef int yiter[2]
    cdef int ziter[2]
    cdef np.float64_t xiterv[2]
    cdef np.float64_t yiterv[2]
    cdef np.float64_t ziterv[2]

    xiter[0] = yiter[0] = ziter[0] = 0
    xiterv[0] = yiterv[0] = ziterv[0] = 0.0

    if period is not None:
        period_x = period[0]
        period_y = period[1]
        period_z = period[2]

    xsize, ysize, zsize = buff.shape[0], buff.shape[1], buff.shape[2]
    x_min = bounds[0]
    x_max = bounds[1]
    y_min = bounds[2]
    y_max = bounds[3]
    z_min = bounds[4]
    z_max = bounds[5]

    dx = (x_max - x_min) / xsize
    dy = (y_max - y_min) / ysize
    dz = (z_max - z_min) / zsize
    idx = 1.0/dx
    idy = 1.0/dy
    idz = 1.0/dz

    kernel_func = get_kernel_func(kernel_name)

    with nogil:
        # TODO make this parallel without using too much memory
        for j in range(0, posx.shape[0]):
            if j % 50000 == 0:
                with gil:
                    if(pbar is not None):
                        pbar.update(50000)
                    PyErr_CheckSignals()

            xiter[1] = yiter[1] = ziter[1] = 999
            xiterv[1] = yiterv[1] = ziterv[1] = 0.0

            if check_period == 1:
                if posx[j] - hsml[j] < x_min:
                    xiter[1] = +1
                    xiterv[1] = period_x
                elif posx[j] + hsml[j] > x_max:
                    xiter[1] = -1
                    xiterv[1] = -period_x
                if posy[j] - hsml[j] < y_min:
                    yiter[1] = +1
                    yiterv[1] = period_y
                elif posy[j] + hsml[j] > y_max:
                    yiter[1] = -1
                    yiterv[1] = -period_y
                if posz[j] - hsml[j] < z_min:
                    ziter[1] = +1
                    ziterv[1] = period_z
                elif posz[j] + hsml[j] > z_max:
                    ziter[1] = -1
                    ziterv[1] = -period_z

            h_j3 = fmax(hsml[j]*hsml[j]*hsml[j], dx*dy*dz)
            h_j = math.cbrt(h_j3)
            h_j2 = h_j*h_j
            ih_j = 1/h_j

            prefactor_j = pmass[j] / pdens[j] / hsml[j]**3 * quantity_to_smooth[j]

            for ii in range(2):
                if xiter[ii] == 999: continue
                px = posx[j] + xiterv[ii]
                if (px + hsml[j] < x_min) or (px - hsml[j] > x_max): continue
                for jj in range(2):
                    if yiter[jj] == 999: continue
                    py = posy[j] + yiterv[jj]
                    if (py + hsml[j] < y_min) or (py - hsml[j] > y_max): continue
                    for kk in range(2):
                        if ziter[kk] == 999: continue
                        pz = posz[j] + ziterv[kk]
                        if (pz + hsml[j] < z_min) or (pz - hsml[j] > z_max): continue

                        x0 = <np.int64_t> ( (px - hsml[j] - x_min) * idx)
                        x1 = <np.int64_t> ( (px + hsml[j] - x_min) * idx)
                        x0 = iclip(x0-1, 0, xsize)
                        x1 = iclip(x1+1, 0, xsize)

                        y0 = <np.int64_t> ( (py - hsml[j] - y_min) * idy)
                        y1 = <np.int64_t> ( (py + hsml[j] - y_min) * idy)
                        y0 = iclip(y0-1, 0, ysize)
                        y1 = iclip(y1+1, 0, ysize)

                        z0 = <np.int64_t> ( (pz - hsml[j] - z_min) * idz)
                        z1 = <np.int64_t> ( (pz + hsml[j] - z_min) * idz)
                        z0 = iclip(z0-1, 0, zsize)
                        z1 = iclip(z1+1, 0, zsize)

                        # Now we know which voxels to deposit onto for this particle,
                        # so loop over them and add this particle's contribution
                        for xi in range(x0, x1):
                            x = (xi + 0.5) * dx + x_min

                            posx_diff = px - x
                            posx_diff = posx_diff * posx_diff
                            if posx_diff > h_j2:
                                continue

                            for yi in range(y0, y1):
                                y = (yi + 0.5) * dy + y_min

                                posy_diff = py - y
                                posy_diff = posy_diff * posy_diff
                                if posy_diff > h_j2:
                                    continue

                                for zi in range(z0, z1):
                                    z = (zi + 0.5) * dz + z_min

                                    posz_diff = pz - z
                                    posz_diff = posz_diff * posz_diff
                                    if posz_diff > h_j2:
                                        continue

                                    # see equation 4 of the SPLASH paper
                                    q_ij = math.sqrt(posx_diff + posy_diff + posz_diff) * ih_j
                                    if q_ij >= 1:
                                        continue

                                    buff[xi, yi, zi] += prefactor_j * kernel_func(q_ij)


def pixelize_element_mesh_line(np.ndarray[np.float64_t, ndim=2] coords,
                               np.ndarray[np.int64_t, ndim=2] conn,
                               np.ndarray[np.float64_t, ndim=1] start_point,
                               np.ndarray[np.float64_t, ndim=1] end_point,
                               npoints,
                               np.ndarray[np.float64_t, ndim=2] field,
                               int index_offset = 0):

    # This routine chooses the correct element sampler to interpolate field
    # values at evenly spaced points along a sampling line
    cdef np.float64_t *vertices
    cdef np.float64_t *field_vals
    cdef int nvertices = conn.shape[1]
    cdef int ndim = coords.shape[1]
    cdef int num_field_vals = field.shape[1]
    cdef int num_plot_nodes = npoints
    cdef int num_intervals = npoints - 1
    cdef double[4] mapped_coord
    cdef ElementSampler sampler
    cdef np.ndarray[np.float64_t, ndim=1] lin_vec
    cdef np.ndarray[np.float64_t, ndim=1] lin_inc
    cdef np.ndarray[np.float64_t, ndim=2] lin_sample_points
    cdef np.int64_t i, n, j, k
    cdef np.ndarray[np.float64_t, ndim=1] arc_length
    cdef np.float64_t lin_length, inc_length
    cdef np.ndarray[np.float64_t, ndim=1] plot_values
    cdef np.float64_t sample_point[3]

    lin_vec = np.zeros(ndim, dtype="float64")
    lin_inc = np.zeros(ndim, dtype="float64")

    lin_sample_points = np.zeros((num_plot_nodes, ndim), dtype="float64")
    arc_length = np.zeros(num_plot_nodes, dtype="float64")
    plot_values = np.zeros(num_plot_nodes, dtype="float64")

    # Pick the right sampler and allocate storage for the mapped coordinate
    if ndim == 3 and nvertices == 4:
        sampler = P1Sampler3D()
    elif ndim == 3 and nvertices == 6:
        sampler = W1Sampler3D()
    elif ndim == 3 and nvertices == 8:
        sampler = Q1Sampler3D()
    elif ndim == 3 and nvertices == 20:
        sampler = S2Sampler3D()
    elif ndim == 2 and nvertices == 3:
        sampler = P1Sampler2D()
    elif ndim == 1 and nvertices == 2:
        sampler = P1Sampler1D()
    elif ndim == 2 and nvertices == 4:
        sampler = Q1Sampler2D()
    elif ndim == 2 and nvertices == 9:
        sampler = Q2Sampler2D()
    elif ndim == 2 and nvertices == 6:
        sampler = T2Sampler2D()
    elif ndim == 3 and nvertices == 10:
        sampler = Tet2Sampler3D()
    else:
        raise YTElementTypeNotRecognized(ndim, nvertices)

    # allocate temporary storage
    vertices = <np.float64_t *> malloc(ndim * sizeof(np.float64_t) * nvertices)
    field_vals = <np.float64_t *> malloc(sizeof(np.float64_t) * num_field_vals)

    lin_vec = end_point - start_point
    lin_length = np.linalg.norm(lin_vec)
    lin_inc = lin_vec / num_intervals
    inc_length = lin_length / num_intervals
    for j in range(ndim):
        lin_sample_points[0, j] = start_point[j]
    arc_length[0] = 0
    for i in range(1, num_intervals + 1):
        for j in range(ndim):
            lin_sample_points[i, j] = lin_sample_points[i-1, j] + lin_inc[j]
            arc_length[i] = arc_length[i-1] + inc_length

    for i in range(num_intervals + 1):
        for j in range(3):
            if j < ndim:
                sample_point[j] = lin_sample_points[i][j]
            else:
                sample_point[j] = 0
        for ci in range(conn.shape[0]):
            for n in range(num_field_vals):
                field_vals[n] = field[ci, n]

            # Fill the vertices
            for n in range(nvertices):
                cj = conn[ci, n] - index_offset
                for k in range(ndim):
                    vertices[ndim*n + k] = coords[cj, k]

            sampler.map_real_to_unit(mapped_coord, vertices, sample_point)
            if not sampler.check_inside(mapped_coord) and ci != conn.shape[0] - 1:
                continue
            elif not sampler.check_inside(mapped_coord):
                raise ValueError("Check to see that both starting and ending line points "
                                 "are within the domain of the mesh.")
            plot_values[i] = sampler.sample_at_unit_point(mapped_coord, field_vals)
            break

    free(vertices)
    free(field_vals)
    return arc_length, plot_values


@cython.boundscheck(False)
@cython.wraparound(False)
def off_axis_projection_SPH(np.float64_t[:] px,
                            np.float64_t[:] py,
                            np.float64_t[:] pz,
                            np.float64_t[:] particle_masses,
                            np.float64_t[:] particle_densities,
                            np.float64_t[:] smoothing_lengths,
                            bounds,
                            center,
                            width,
                            np.float64_t[:] quantity_to_smooth,
                            np.float64_t[:, :] projection_array,
                            normal_vector,
                            north_vector,
                            weight_field=None):
    # Do nothing in event of a 0 normal vector
    if np.allclose(normal_vector, np.array([0., 0., 0.]), rtol=1e-09):
        return

    # We want to do two rotations, one to first rotate our coordinates to have
    # the normal vector be the z-axis (i.e., the viewer's perspective), and then
    # another rotation to make the north-vector be the y-axis (i.e., north).
    # Fortunately, total_rotation_matrix = rotation_matrix_1 x rotation_matrix_2
    cdef int num_particles = np.size(px)
    cdef np.float64_t[:] z_axis = np.array([0., 0., 1.], dtype='float_')
    cdef np.float64_t[:] y_axis = np.array([0., 1., 0.], dtype='float_')
    cdef np.float64_t[:, :] normal_rotation_matrix
    cdef np.float64_t[:] transformed_north_vector
    cdef np.float64_t[:, :] north_rotation_matrix
    cdef np.float64_t[:, :] rotation_matrix

    normal_rotation_matrix = get_rotation_matrix(normal_vector, z_axis)
    transformed_north_vector = np.matmul(normal_rotation_matrix, north_vector)
    north_rotation_matrix = get_rotation_matrix(transformed_north_vector, y_axis)
    rotation_matrix = np.matmul(north_rotation_matrix, normal_rotation_matrix)

    cdef np.float64_t[:] px_rotated = np.empty(num_particles, dtype='float_')
    cdef np.float64_t[:] py_rotated = np.empty(num_particles, dtype='float_')
    cdef np.float64_t[:] coordinate_matrix = np.empty(3, dtype='float_')
    cdef np.float64_t[:] rotated_coordinates
    cdef np.float64_t[:] rotated_center
    rotated_center = rotation_matmul(
        rotation_matrix, np.array([center[0], center[1], center[2]]))

    # set up the rotated bounds
    cdef np.float64_t rot_bounds_x0 = rotated_center[0] - width[0] / 2
    cdef np.float64_t rot_bounds_x1 = rotated_center[0] + width[0] / 2
    cdef np.float64_t rot_bounds_y0 = rotated_center[1] - width[1] / 2
    cdef np.float64_t rot_bounds_y1 = rotated_center[1] + width[1] / 2

    for i in range(num_particles):
        coordinate_matrix[0] = px[i]
        coordinate_matrix[1] = py[i]
        coordinate_matrix[2] = pz[i]
        rotated_coordinates = rotation_matmul(
            rotation_matrix, coordinate_matrix)
        px_rotated[i] = rotated_coordinates[0]
        py_rotated[i] = rotated_coordinates[1]

    pixelize_sph_kernel_projection(projection_array,
                                   px_rotated,
                                   py_rotated,
                                   smoothing_lengths,
                                   particle_masses,
                                   particle_densities,
                                   quantity_to_smooth,
                                   [rot_bounds_x0, rot_bounds_x1,
                                    rot_bounds_y0, rot_bounds_y1],
                                   weight_field=weight_field,
                                   check_period=0)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.float64_t[:] rotation_matmul(np.float64_t[:, :] rotation_matrix,
                                     np.float64_t[:] coordinate_matrix):
    cdef np.float64_t[:] out = np.zeros(3)
    for i in range(3):
        for j in range(3):
            out[i] += rotation_matrix[i, j] * coordinate_matrix[j]
    return out

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.float64_t[:, :] get_rotation_matrix(np.float64_t[:] normal_vector,
                                             np.float64_t[:] final_vector):
    """ Returns a numpy rotation matrix corresponding to the
    rotation of the given normal vector to the specified final_vector.
    See https://math.stackexchange.com/a/476311 although note we return the
    inverse of what's specified there.
    """
    cdef np.float64_t[:] normal_unit_vector = normal_vector / np.linalg.norm(normal_vector)
    cdef np.float64_t[:] final_unit_vector = final_vector / np.linalg.norm(final_vector)
    cdef np.float64_t[:] v = np.cross(final_unit_vector, normal_unit_vector)
    cdef np.float64_t s = np.linalg.norm(v)
    cdef np.float64_t c = np.dot(final_unit_vector, normal_unit_vector)
    # if the normal vector is identical to the final vector, just return the
    # identity matrix
    if np.isclose(c, 1, rtol=1e-09):
        return np.identity(3, dtype='float_')
    # if the normal vector is the negative final vector, return the appropriate
    # rotation matrix for flipping your coordinate system.
    if np.isclose(s, 0, rtol=1e-09):
        return np.array([[0, -1, 0],[-1, 0, 0],[0, 0, -1]], dtype='float_')

    cdef np.float64_t[:, :] cross_product_matrix = np.array([[0, -1 * v[2], v[1]],
                                                      [v[2], 0, -1 * v[0]],
                                                      [-1 * v[1], v[0], 0]],
                                                      dtype='float_')
    return np.linalg.inv(np.identity(3, dtype='float_') + cross_product_matrix
                         + np.matmul(cross_product_matrix, cross_product_matrix)
                         * 1/(1+c))

@cython.boundscheck(False)
@cython.wraparound(False)
def normalization_3d_utility(np.float64_t[:, :, :] num,
                             np.float64_t[:, :, :] den):
    cdef int i, j, k
    for i in range(num.shape[0]):
        for j in range(num.shape[1]):
            for k in range(num.shape[2]):
                if den[i, j, k] != 0.0:
                    num[i, j, k] = num[i, j, k] / den[i, j, k]

@cython.boundscheck(False)
@cython.wraparound(False)
def normalization_2d_utility(np.float64_t[:, :] num,
                          np.float64_t[:, :] den):
    cdef int i, j
    for i in range(num.shape[0]):
        for j in range(num.shape[1]):
            if den[i, j] != 0.0:
                num[i, j] = num[i, j] / den[i, j]

@cython.boundscheck(False)
@cython.wraparound(False)
def normalization_1d_utility(np.float64_t[:] num,
                             np.float64_t[:] den):
    cdef int i
    for i in range(num.shape[0]):
        if den[i] != 0.0:
            num[i] = num[i] / den[i]
