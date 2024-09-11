
# distutils: libraries = STD_LIBS
# distutils: extra_compile_args = OMP_ARGS
# distutils: extra_link_args = OMP_ARGS
"""
Utilities for images
"""


import numpy as np

cimport numpy as np
cimport cython
from libc.math cimport atan2, ceil, floor, log2, sqrt, M_PI, M_PI_2
from libc.stdlib cimport free, malloc

from yt.utilities.lib.fp_utils cimport iclip, imin, imax, fclip, fmin, fmax
from cython.parallel import prange, parallel

cdef np.float64_t twopi = 2 * M_PI

@cython.boundscheck(False)
cdef inline void _add_point_to_greyscale_image(
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    np.float64_t px,
    np.float64_t py,
    np.float64_t pv,
    int xs,
    int ys
) noexcept nogil:
    cdef int i, j
    j = (<int> (xs * px)) % xs
    i = (<int> (ys * py)) % ys
    if (i < 0) or (i >= buffer.shape[0]) or (j < 0) or (j >= buffer.shape[1]):
        # some particles might intersect the image buffer
        # but actually be centered out of bounds. Skip those.
        # see https://github.com/yt-project/yt/issues/4603
        return
    buffer[i, j] += pv
    buffer_mask[i, j] = 1



cdef inline np.float64_t  _wrap_dist(np.float64_t dx) nogil:
    if dx > 0.5:
        return 1 - dx
    elif dx < -0.5:
        return 1 + dx
    else:
        return dx


@cython.wraparound(False)
def add_points_to_greyscale_image(
        np.ndarray[np.float64_t, ndim=2] buffer,
        np.ndarray[np.uint8_t,   ndim=2] buffer_mask,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.ndarray[np.float64_t, ndim=1] pv):
    cdef int pi
    cdef int npx = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]

    cdef np.float64_t[:, ::1] buffer_view = buffer
    cdef np.uint8_t[:, ::1] buffer_mask_view = buffer_mask

    for pi in range(npx):
        _add_point_to_greyscale_image(buffer_view, buffer_mask_view, px[pi], py[pi], pv[pi], xs, ys)
    return

cdef inline int ij2idx(const int i, const int j, const int Nx) noexcept nogil:
    return i * Nx + j

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
cdef void _add_cell_to_image_offaxis(
    const np.float64_t dx,
    const np.float64_t w,
    const np.float64_t q,
    const np.float64_t cell_max_width,
    const np.float64_t x,
    const np.float64_t y,
    const int Nx,
    const int Ny,
    const int Nsx,
    const int Nsy,
    np.float64_t* buffer,
    np.float64_t* buffer_weight,
    const int max_depth,
    const np.float64_t[:, :, ::1] stamp,
    const np.uint8_t[:, :, ::1] stamp_mask,
) noexcept nogil:
    cdef np.float64_t lx, rx, ly, ry
    cdef int j, k, depth
    cdef int jmin, jmax, kmin, kmax, jj1, jj2, kk1, kk2, itmp
    cdef np.float64_t cell_max_half_width = cell_max_width / 2
    cdef np.float64_t xx1, xx2, yy1, yy2, dvx, dvy, sw, sq, tmp, dx_loc, dy_loc
    cdef np.float64_t dx3 = dx * dx * dx * Nx * Ny

    lx = x - cell_max_half_width
    rx = x + cell_max_half_width

    ly = y - cell_max_half_width
    ry = y + cell_max_half_width

    # Compute the range of pixels that the cell may overlap
    jmin = imax(<int>floor(lx * Nx), 0)
    jmax = imax(<int>ceil(rx * Nx), Nx - 1)

    kmin = imin(<int>floor(ly * Ny), 0)
    kmax = imax(<int>ceil(ry * Ny), Ny - 1)

    # If the cell is fully within one pixel
    if (jmax == jmin + 1) and (kmax == kmin + 1):
        buffer[ij2idx(jmin, kmin, Nx)]        += q * dx3
        buffer_weight[ij2idx(jmin, kmin, Nx)] += w * dx3
        return

    # Our 'stamp' has multiple resolutions, select the one
    # that is at a higher resolution than the pixel
    # we are projecting onto with at least 4 pixels on the diagonal
    depth = iclip(
        <int> (ceil(log2(4 * sqrt(3) * dx * fmax(Nx, Ny)))),
        1,
        max_depth - 1,
    )

    jmax = imin(Nsx, 1 << depth)
    kmax = imin(Nsy, 1 << depth)

    dx_loc = cell_max_width / jmax
    dy_loc = cell_max_width / kmax

    for j in range(jmax):
        xx1 = ((j - jmax / 2.) * dx_loc + x) * Nx
        xx2 = ((j + 1 - jmax / 2.) * dx_loc + x) * Nx

        jj1 = <int> xx1
        jj2 = <int> xx2

        # The subcell is out of the projected area
        if jj2 < 0 or jj1 >= Nx: continue

        # Fraction of overlap with the pixel in x direction
        dvx = fclip((jj2 - xx1) / (xx2 - xx1), 0., 1.)

        for k in range(kmax):
            if stamp_mask[depth, j, k] == 0:
                continue

            yy1 = ((k - kmax / 2.) * dy_loc + y) * Ny
            yy2 = ((k + 1 - kmax / 2.) * dy_loc + y) * Ny

            kk1 = <int> yy1
            kk2 = <int> yy2
            # The subcell is out of the projected area
            if kk2 < 0 or kk1 >= Ny: continue

            tmp = stamp[depth, j, k] * dx3
            sw = tmp * w
            sq = tmp * q

            # Fraction of overlap with the pixel in y direction
            dvy = fclip((kk2 - yy1) / (yy2 - yy1), 0., 1.)

            if jj1 >= 0 and kk1 >= 0:
                tmp = dvx * dvy
                itmp = ij2idx(jj1, kk1, Nx)
                buffer[itmp]        += sq * tmp
                buffer_weight[itmp] += sw * tmp

            if jj1 >= 0 and kk2 < Ny:
                tmp = dvx * (1 - dvy)
                itmp = ij2idx(jj1, kk2, Nx)
                buffer[itmp]        += sq * tmp
                buffer_weight[itmp] += sw * tmp

            if jj2 < Nx and kk1 >= 0:
                tmp = (1 - dvx) * dvy
                itmp = ij2idx(jj2, kk1, Nx)
                buffer[itmp]        += sq * tmp
                buffer_weight[itmp] += sw * tmp

            if jj2 < Nx and kk2 < Ny:
                tmp = (1 - dvx) * (1 - dvy)
                itmp = ij2idx(jj2, kk2, Nx)
                buffer[itmp]        += sq * tmp
                buffer_weight[itmp] += sw * tmp

@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def add_cells_to_image_offaxis(
    *,
    const np.float64_t[:, ::1] Xp,
    const np.float64_t[::1] dXp,
    const np.float64_t[::1] qty,
    const np.float64_t[::1] weight,
    const np.float64_t[:, :] rotation,
    np.float64_t[:, ::1] buffer,
    np.float64_t[:, ::1] buffer_weight,
    const int Nx,
    const int Ny,
    const int Npix_min = 4,
):
    cdef np.ndarray[np.float64_t, ndim=1] center = np.array([0.5, 0.5, 0.5])
    cdef np.float64_t w0 = 1 / sqrt(3.)
    cdef int i, j, k

    cdef np.ndarray[np.float64_t, ndim=1] a = np.array([1., 0, 0]) * w0
    cdef np.ndarray[np.float64_t, ndim=1] b = np.array([0, 1., 0]) * w0
    cdef np.ndarray[np.float64_t, ndim=1] c = np.array([0, 0, 1.]) * w0

    a = np.dot(rotation, a)
    b = np.dot(rotation, b)
    c = np.dot(rotation, c)

    cdef np.ndarray[np.float64_t, ndim=1] o = center - (a + b + c) / 2

    cdef int Nsx, Nsy
    cdef np.float64_t dx_max = np.max(dXp)
    # The largest cell needs to be resolved by at least this number of pixels
    Nsx = max(Npix_min, int(ceil(2 * dx_max * sqrt(3) * Nx)))
    Nsy = max(Npix_min, int(ceil(2 * dx_max * sqrt(3) * Ny)))
    cdef int max_depth = int(ceil(log2(max(Nsx, Nsy))))
    cdef int depth

    cdef np.ndarray[np.float64_t, ndim=3] stamp_arr = np.zeros((max_depth, Nsx, Nsy), dtype=float)
    cdef np.ndarray[np.uint8_t, ndim=3] stamp_mask_arr = np.zeros((max_depth, Nsx, Nsy), dtype=np.uint8)
    cdef np.float64_t[:, :, ::1] stamp = stamp_arr
    cdef np.uint8_t[:, :, ::1] stamp_mask = stamp_mask_arr

    # Precompute the mip
    for depth in range(max_depth):
        if depth == 0:
            stamp[0, 0, 0] = 1
            continue
        direct_integrate_cube(
            o,
            a,
            b,
            c,
            stamp_arr[depth, :, :],
            stamp_mask_arr[depth, :, :],
            imin(1 << depth, Nsx),
            imin(1 << depth, Nsy),
        )

        stamp_arr[depth] /= np.sum(stamp[depth])

    # Iterate over all cells, applying the stamp
    cdef np.float64_t x, y, dx
    cdef np.float64_t[:, ::1] rotation_view = np.ascontiguousarray(rotation)
    cdef np.float64_t w, q, cell_max_width, sq3

    sq3 = sqrt(3.)

    # Local buffers
    cdef np.float64_t *lbuffer
    cdef np.float64_t *lbuffer_weight

    cdef int num_particles = len(Xp)

    with nogil, parallel():
        lbuffer = <np.float64_t*> malloc(sizeof(np.float64_t*) * Nx * Ny)
        lbuffer_weight = <np.float64_t*> malloc(sizeof(np.float64_t*) * Nx * Ny)
        for j in range(Nx * Ny):
            lbuffer[j] = 0
            lbuffer_weight[j] = 0

        for i in prange(num_particles, schedule="runtime"):
            dx = dXp[i]
            w = weight[i]
            q = qty[i]
            cell_max_width = dx * sq3
            x = (
                rotation_view[0, 0] * Xp[i, 0] +
                rotation_view[0, 1] * Xp[i, 1] +
                rotation_view[0, 2] * Xp[i, 2]
            ) + 0.5
            y = (
                rotation_view[1, 0] * Xp[i, 0] +
                rotation_view[1, 1] * Xp[i, 1] +
                rotation_view[1, 2] * Xp[i, 2]
            ) + 0.5

            _add_cell_to_image_offaxis(
                dx,
                w,
                q,
                cell_max_width,
                x,
                y,
                Nx,
                Ny,
                Nsx,
                Nsy,
                lbuffer,
                lbuffer_weight,
                max_depth,
                stamp,
                stamp_mask
            )

        # Copy back data in main buffer
        with gil:
            for j in range(Nx):
                for k in range(Ny):
                    buffer[j, k] += lbuffer[ij2idx(j, k, Nx)]
                    buffer_weight[j, k] += lbuffer_weight[ij2idx(j, k, Nx)]
        # Free memory
        free(lbuffer)
        free(lbuffer_weight)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.float64_t det2d(const np.float64_t[::1] a, const np.float64_t[::1] b) noexcept nogil:
    return a[0] * b[1] - a[1] * b[0]

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef bint check_in_parallelogram(
    const np.float64_t[::1] PA,
    const np.float64_t[::1] PQ,
    const np.float64_t[::1] PR,
    const int signPQ,
    const int signPR,
    np.float64_t[2] out
) noexcept nogil:
    cdef np.float64_t det_PQR = det2d(PQ, PR)
    if det_PQR == 0:
        out[0] = -1
        out[1] = -1
        return False

    out[0] = -det2d(PA, PQ) / det_PQR
    out[1] = det2d(PA, PR) / det_PQR

    if 0 <= signPQ * out[0] <= 1 and 0 <= signPR * out[1] <= 1:
        return True

    out[0] = -1
    out[1] = -1
    return False

@cython.boundscheck(False)
@cython.cdivision(True)
cdef int direct_integrate_cube(
    np.ndarray[np.float64_t, ndim=1] O,
    np.ndarray[np.float64_t, ndim=1] u,
    np.ndarray[np.float64_t, ndim=1] v,
    np.ndarray[np.float64_t, ndim=1] w,
    np.ndarray[np.float64_t, ndim=2] buffer,
    np.ndarray[np.uint8_t, ndim=2] buffer_mask,
    const int Nx,
    const int Ny,
) except -1:
    """
    Compute depth of cube from direct integration of entry/exit points of rays
    """
    cdef np.float64_t[::1] u2d = u[:2]
    cdef np.float64_t[::1] v2d = v[:2]
    cdef np.float64_t[::1] w2d = w[:2]

    cdef np.float64_t[::1] Oback = O + u + v + w

    cdef np.float64_t[::1] X = np.zeros(2)
    cdef np.float64_t[::1] OfrontA = np.zeros(2)
    cdef np.float64_t[::1] ObackA = np.zeros(2)

    cdef np.float64_t inv_dx = 1. / Nx
    cdef np.float64_t inv_dy = 1. / Ny
    cdef np.float64_t[2] nm
    cdef bint within
    cdef np.float64_t zmin, zmax, z
    cdef int Nhit, i, j
    for i in range(Nx):
        X[0] = (i + 0.5) * inv_dx

        OfrontA[0] = X[0] - O[0]
        ObackA[0] = X[0] - Oback[0]

        for j in range(Ny):
            zmin = np.inf
            zmax = -np.inf
            Nhit = 0
            X[1] = (j + 0.5) * inv_dy

            OfrontA[1] = X[1] - O[1]
            ObackA[1] = X[1] - Oback[1]

            within = check_in_parallelogram(OfrontA, v2d, u2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * u[2] + nm[1] * v[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(OfrontA, w2d, v2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * v[2] + nm[1] * w[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(OfrontA, w2d, u2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * u[2] + nm[1] * w[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(ObackA, v2d, u2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * u[2] + nm[1] * v[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(ObackA, w2d, v2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * v[2] + nm[1] * w[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(ObackA, w2d, u2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * u[2] + nm[1] * w[2]
                zmin = fmin(z, zmin)
                zmax = fmax(z, zmax)
                Nhit += 1

            if Nhit == 0:
                continue
            elif Nhit == 1:
                raise RuntimeError("This should not happen")
            else:
                buffer[i, j] += zmax - zmin
                buffer_mask[i, j] = 1

def add_points_to_image(
        np.ndarray[np.uint8_t, ndim=3] buffer,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.float64_t pv):
    cdef int i, j, k, pi
    cdef int npx = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    cdef int v
    v = iclip(<int>(pv * 255), 0, 255)
    for pi in range(npx):
        j = <int> (xs * px[pi])
        i = <int> (ys * py[pi])
        for k in range(3):
            buffer[i, j, k] = v
        buffer[i, j, 3] = 255
    return

def add_rgba_points_to_image(
        np.ndarray[np.float64_t, ndim=3] buffer,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.ndarray[np.float64_t, ndim=2] rgba,
        ):
    """
    Splat rgba points onto an image

    Given an image buffer, add colors to
    pixels defined by fractional positions px and py,
    with colors rgba.  px and py are one dimensional
    arrays, and rgba is a an array of rgba values.
    """
    cdef int i, j, k, pi
    cdef int npart = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    #iv = iclip(<int>(pv * 255), 0, 255)
    for pi in range(npart):
        j = <int> (xs * px[pi])
        i = <int> (ys * py[pi])
        if i < 0 or j < 0 or i >= xs or j >= ys:
            continue
        for k in range(4):
            buffer[i, j, k] += rgba[pi, k]
    return

@cython.boundscheck(False)
cdef np.ndarray[np.float64_t, ndim=2] sample_tetrahedron(int sample):
    cdef int i
    cdef np.float64_t a, b, c

    cdef np.ndarray[np.float64_t, ndim=2] buffer_out = np.zeros((sample, 3))
    cdef np.float64_t[:, :] buffer = buffer_out

    for i in range(sample):
        a = np.random.rand()
        b = np.random.rand()
        c = np.random.rand()

        while a + b + c > 1:
            a = np.random.rand()
            b = np.random.rand()
            c = np.random.rand()

        buffer[i, 0] = a
        buffer[i, 1] = b
        buffer[i, 2] = c

    return buffer_out


@cython.boundscheck(False)
cdef int convex_hull(np.float64_t[4][2] ABCD, int[4] hull) nogil:
    """Compute the convex hull of ABCD.
    Parameters
    ----------
    ABCD : (4, 2) array_like
        The four points
    Returns
    -------
    hull : 2D array_like
        The convex hull of ABCD
    """
    cdef int i, ilo
    ilo = 0
    for i in range(1, 4):
        if ABCD[i][0] < ABCD[ilo][0]:
            ilo = i
    hull[0] = ilo
    cdef np.uint8_t[4] mask
    mask[0] = 1
    mask[1] = 1
    mask[2] = 1
    mask[3] = 1
    mask[hull[0]] = 0

    # Iterate until we get back to the first point
    cdef np.float64_t prev_angle = M_PI_2
    cdef np.float64_t min_diff, diff, angle, best_angle
    cdef int best_j

    for i in range(1, 4):
        min_diff = twopi
        # Find point with smallest angle compared to previous segment
        for j in range(4):
            if j == hull[i-1]:  # Skip previous point
                continue

            angle = atan2(ABCD[j][1] - ABCD[hull[i-1]][1], ABCD[j][0] - ABCD[hull[i-1]][0])
            diff = (prev_angle - angle) % twopi
            if diff < min_diff:
                min_diff = diff
                best_j = j
                best_angle = angle

        if best_j == hull[0]:
            break

        hull[i] = best_j
        mask[best_j] = 0
        prev_angle = best_angle

    for i in range(4):
        if mask[i]:
            hull[3] = i
            return 3
    return 4

@cython.boundscheck(False)
cdef int integrate_quad(
    np.float64_t[2] A,
    np.float64_t[2] B,
    np.float64_t[2] C,
    np.float64_t[2] P,
    np.float64_t value,
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    int xs,
    int ys,
) nogil:
    # Find bounding box of ABC
    cdef int min_x = <int> (floor(xs * min(A[0], B[0], C[0])))
    cdef int max_x = <int> (ceil(xs * max(A[0], B[0], C[0])))
    cdef int min_y = <int> (floor(ys * min(A[1], B[1], C[1])))
    cdef int max_y = <int> (ceil(ys * max(A[1], B[1], C[1])))

    cdef np.float64_t[2] AB
    cdef np.float64_t[2] CA
    cdef np.float64_t[2] BC
    AB[0] = B[0] - A[0]
    AB[1] = B[1] - A[1]
    CA[0] = A[0] - C[0]
    CA[1] = A[1] - C[1]
    BC[0] = C[0] - B[0]
    BC[1] = C[1] - B[1]

    cdef np.float64_t[2] PA
    cdef np.float64_t[2] PB
    cdef np.float64_t[2] PC
    PA[0] = A[0] - P[0]
    PA[1] = A[1] - P[1]
    PB[0] = B[0] - P[0]
    PB[1] = B[1] - P[1]
    PC[0] = C[0] - P[0]
    PC[1] = C[1] - P[1]

    cdef np.float64_t[2] X
    X[0] = 0
    X[1] = 0

    cdef int i, j
    cdef np.float64_t[2] PX

    cdef np.float64_t aa, bb, gg, alpha, beta, gamma, v, Vtot
    Vtot = 0
    # First pass: compute integral within ABC
    for i in range(min_x, max_x):
        X[0] = (i + 0.5) / xs
        for j in range(min_y, max_y):
            X[1] = (j + 0.5) / ys
            PX[0] = X[0] - P[0]
            PX[1] = X[1] - P[1]

            # Compute alpha = (PX x BC) / (PB x BC)
            if (PB[1] * BC[0] - PB[0] * BC[1]) == 0:
                continue  # raise ValueError("PB x BC is zero")
            alpha = 1 - (PX[1] * BC[0] - PX[0] * BC[1]) / (PB[1] * BC[0] - PB[0] * BC[1])
            # Compute beta = (PX x CA) / (PC x CA)
            if (PC[1] * CA[0] - PC[0] * CA[1]) == 0:
                continue  # raise ValueError("PC x CA is zero")
            beta = 1 - (PX[1] * CA[0] - PX[0] * CA[1]) / (PC[1] * CA[0] - PC[0] * CA[1])
            # Compute gamma = (PX x AB) / (PA x AB)
            if (PA[1] * AB[0] - PA[0] * AB[1]) == 0:
                continue  # raise ValueError("PA x AB is zero")
            gamma = 1 - (PX[1] * AB[0] - PX[0] * AB[1]) / (PA[1] * AB[0] - PA[0] * AB[1])

            # Check whether X is within ABC
            aa = 0 <= alpha
            bb = 0 <= beta
            gg = 0 <= gamma

            if not (aa and bb and gg):
                continue

            v = 2
            if 0 <= alpha < 1:
                v = alpha
            if 0 <= beta < 1:
                v = min(v, beta)
            if 0 <= gamma < 1:
                v = min(v, gamma)

            Vtot += v

    # Special case: Vtot == 0
    if Vtot == 0:
        _add_point_to_greyscale_image(
            buffer, buffer_mask, P[0], P[1], value, xs, ys
        )
        return 0

    # Second pass: deposit
    for i in range(min_x, max_x):
        X[0] = (i + 0.5) / xs
        for j in range(min_y, max_y):
            X[1] = (j + 0.5) / ys
            PX[0] = X[0] - P[0]
            PX[1] = X[1] - P[1]

            # Compute alpha = (PX x BC) / (PB x BC)
            alpha = 1 - (PX[1] * BC[0] - PX[0] * BC[1]) / (PB[1] * BC[0] - PB[0] * BC[1])
            # Compute beta = (PX x CA) / (PC x CA)
            beta = 1 - (PX[1] * CA[0] - PX[0] * CA[1]) / (PC[1] * CA[0] - PC[0] * CA[1])
            # Compute gamma = (PX x AB) / (PA x AB)
            gamma = 1 - (PX[1] * AB[0] - PX[0] * AB[1]) / (PA[1] * AB[0] - PA[0] * AB[1])

            # Check whether X is within ABC
            aa = 0 <= alpha
            bb = 0 <= beta
            gg = 0 <= gamma

            if not (aa and bb and gg):
                continue

            v = 2
            if 0 <= alpha < 1:
                v = alpha
            if 0 <= beta < 1:
                v = min(v, beta)
            if 0 <= gamma < 1:
                v = min(v, gamma)

            buffer[j % ys, i % xs] += value * v / Vtot
            buffer_mask[j % ys, i % xs] = 1

@cython.boundscheck(False)
cdef int integrate_tetrahedron_proj(
    np.float64_t[3] origin,
    np.float64_t[3] a,
    np.float64_t[3] b,
    np.float64_t[3] c,
    np.float64_t value,
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    int xs,
    int ys
) nogil:
    cdef np.float64_t[4][2] ABCD_xy
    cdef np.float64_t[2] A_xy, B_xy, C_xy, D_xy, P_xy, P1_xy, P2_xy
    cdef np.float64_t det, S1, S2

    ABCD_xy[0][0] = origin[0]
    ABCD_xy[0][1] = origin[1]
    ABCD_xy[1][0] = origin[0] + a[0]
    ABCD_xy[1][1] = origin[1] + a[1]
    ABCD_xy[2][0] = origin[0] + b[0]
    ABCD_xy[2][1] = origin[1] + b[1]
    ABCD_xy[3][0] = origin[0] + c[0]
    ABCD_xy[3][1] = origin[1] + c[1]
    P_xy[0] = 0
    P_xy[1] = 0

    cdef int[4] iABC

    cdef int Nhull = convex_hull(ABCD_xy, iABC)

    A_xy[0] = ABCD_xy[iABC[0]][0]
    A_xy[1] = ABCD_xy[iABC[0]][1]
    B_xy[0] = ABCD_xy[iABC[1]][0]
    B_xy[1] = ABCD_xy[iABC[1]][1]
    C_xy[0] = ABCD_xy[iABC[2]][0]
    C_xy[1] = ABCD_xy[iABC[2]][1]
    D_xy[0] = ABCD_xy[iABC[3]][0]
    D_xy[1] = ABCD_xy[iABC[3]][1]
    if Nhull == 4:
        # Find the intersection of AC with BD
        det = ((A_xy[0] - C_xy[0]) * (B_xy[1] - D_xy[1]) - (A_xy[1] - C_xy[1]) * (B_xy[0] - D_xy[0]))
        P_xy[0] = ((A_xy[0] * C_xy[1] - A_xy[1] * C_xy[0]) * (B_xy[0] - D_xy[0]) - (A_xy[0] - C_xy[0]) * (B_xy[0] * D_xy[1] - B_xy[1] * D_xy[0])) / det
        P_xy[1] = ((A_xy[0] * C_xy[1] - A_xy[1] * C_xy[0]) * (B_xy[1] - D_xy[1]) - (A_xy[1] - C_xy[1]) * (B_xy[0] * D_xy[1] - B_xy[1] * D_xy[0])) / det

        # Compute (double the) surface of each triangle
        # Surface of triangle ABC
        S1 = (C_xy[1] - A_xy[1]) * (B_xy[0] - A_xy[0]) - (C_xy[0] - A_xy[0]) * (B_xy[1] - A_xy[1])
        # Surface of triangle ACD
        S2 = (D_xy[1] - A_xy[1]) * (C_xy[0] - A_xy[0]) - (D_xy[0] - A_xy[0]) * (C_xy[1] - A_xy[1])

        # Move slightly towards B to prevent rounding errors
        P1_xy[0] = P_xy[0] + 1e-10 * (B_xy[0] - P_xy[0])
        P1_xy[1] = P_xy[1] + 1e-10 * (B_xy[1] - P_xy[1])
        integrate_quad(A_xy, B_xy, C_xy, P1_xy, value * S1 / (S1 + S2), buffer, buffer_mask, xs, ys)

        # Move slightly towards D to prevent rounding errors
        P2_xy[0] = P_xy[0] + 1e-10 * (D_xy[0] - P_xy[0])
        P2_xy[1] = P_xy[1] + 1e-10 * (D_xy[1] - P_xy[1])
        integrate_quad(A_xy, C_xy, D_xy, P2_xy, value * S2 / (S1 + S2), buffer, buffer_mask, xs, ys)
    else:
        integrate_quad(A_xy, B_xy, C_xy, D_xy, value, buffer, buffer_mask, xs, ys)

# @cython.boundscheck(False)
def add_points_to_greyscale_image_with_lagrangian_tesselation(
    np.ndarray[np.float64_t, ndim=2] buffer,
    np.ndarray[np.uint8_t,   ndim=2] buffer_mask,
    np.ndarray[np.float64_t, ndim=4] p3d,
    np.ndarray[np.float64_t, ndim=3] pv,
):
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    cdef int Ngrid = p3d.shape[0]

    cdef np.float64_t[:, ::1] buffer_view = buffer
    cdef np.uint8_t[:, ::1] buffer_mask_view = buffer_mask
    cdef np.float64_t[:, :, :, ::1] p3d_view = p3d

    cdef np.ndarray[np.uint8_t, ndim=2] vert = np.array(
        (
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
            (0, 1, 1),
        ),
        dtype=np.uint8
    )
    cdef np.ndarray[np.uint8_t, ndim=2] conn = np.array(
        (
            (4, 0, 7, 1),
            (1, 0, 3, 7),
            (5, 1, 4, 7),
            (2, 3, 7, 1),
            (1, 5, 6, 7),
            (2, 6, 1, 7),
        ),
        dtype=np.uint8
    )

    cdef int i, m
    cdef np.uint8_t[:, ::1] off

    for m in range(len(conn)):
        # ro = sample_tetrahedron(split)
        off = vert[conn[m]]

        for i in prange(Ngrid, nogil=True):
            # print(f"{m=}/{len(conn)} {i=}/{Ngrid}")
            _lagrangian_tesselation_helper(
                Ngrid,
                i,
                p3d_view,
                off,
                buffer_view,
                buffer_mask_view,
                xs,
                ys
            )


@cython.boundscheck(False)
cdef int _lagrangian_tesselation_helper(
    const int Ngrid,
    const int i,
    const np.float64_t[:, :, :, ::1] p3d_view,
    const np.uint8_t[:, ::1] off,
    np.float64_t[:, ::1] buffer_view,
    np.uint8_t[:, ::1] buffer_mask_view,
    const int xs,
    const int ys,
) nogil:
    cdef int j, k, idim
    cdef np.float64_t[3] orig
    cdef np.float64_t[3] a
    cdef np.float64_t[3] b
    cdef np.float64_t[3] c

    for j in range(Ngrid):
        for k in range(Ngrid):
            for idim in range(3):   # stop at 2, because z is not used
                orig[idim] = p3d_view[(i + off[3][0]) % Ngrid, (j + off[3][1]) % Ngrid, (k + off[3][2]) % Ngrid, idim]
                a[idim] = _wrap_dist(p3d_view[(i + off[0][0]) % Ngrid, (j + off[0][1]) % Ngrid, (k + off[0][2]) % Ngrid, idim] - orig[idim])
                b[idim] = _wrap_dist(p3d_view[(i + off[1][0]) % Ngrid, (j + off[1][1]) % Ngrid, (k + off[1][2]) % Ngrid, idim] - orig[idim])
                c[idim] = _wrap_dist(p3d_view[(i + off[2][0]) % Ngrid, (j + off[2][1]) % Ngrid, (k + off[2][2]) % Ngrid, idim] - orig[idim])

            integrate_tetrahedron_proj(
                orig,
                a,
                b,
                c,
                1,
                buffer_view,
                buffer_mask_view,
                xs,
                ys
            )
