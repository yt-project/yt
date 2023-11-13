
# distutils: libraries = STD_LIBS
# distutils: extra_compile_args = OMP_ARGS
# distutils: extra_link_args = OMP_ARGS
"""
Utilities for images
"""


import numpy as np

cimport numpy as np
cimport cython
from libc.math cimport M_PI, M_PI_2, atan2, ceil, floor, log2, sqrt
from libc.stdlib cimport free, malloc

from yt.utilities.lib.fp_utils cimport iclip
from cython.parallel import prange, parallel

@cython.wraparound(False)
def add_points_to_greyscale_image(
        np.ndarray[np.float64_t, ndim=2] buffer,
        np.ndarray[np.uint8_t,   ndim=2] buffer_mask,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.ndarray[np.float64_t, ndim=1] pv):
    cdef int i, j, pi
    cdef int npx = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    for pi in range(npx):
        j = <int> (xs * px[pi])
        i = <int> (ys * py[pi])
        if (i < 0) or (i >= buffer.shape[0]) or (j < 0) or (j >= buffer.shape[1]):
            # some particles might intersect the image buffer
            # but actually be centered out of bounds. Skip those.
            # see https://github.com/yt-project/yt/issues/4603
            continue

        buffer[i, j] += pv[pi]
        buffer_mask[i, j] = 1
    return

# @cython.boundscheck(False)
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
    np.float64_t** buffer,
    np.float64_t** buffer_weight,
    const int max_depth,
    const np.float64_t[:, :, ::1] stamp,
    const np.uint8_t[:, :, ::1] stamp_mask,
) noexcept nogil:
    cdef np.float64_t lx, rx, ly, ry
    cdef int j, k, depth
    cdef int jmin, jmax, kmin, kmax, jj, kk
    cdef np.float64_t cell_max_half_width = cell_max_width / 2
    cdef np.float64_t xx, yy, dvx, dvy, sw, swq

    lx = x - cell_max_half_width
    rx = x + cell_max_half_width

    ly = y - cell_max_half_width
    ry = y + cell_max_half_width

    # Compute the range of pixels that the cell may overlap
    jmin = int(floor(lx * Nx))
    if jmin < 0: jmin = 0
    jmax = int(ceil(rx * Nx))
    if jmax >= Nx: jmax = Nx - 1

    kmin = int(floor(ly * Ny))
    if kmin < 0: kmin = 0
    kmax = int(ceil(ry * Ny))
    if kmax >= Ny: kmax = Ny - 1

    # If the cell is fully within one pixel
    if (jmax == jmin + 1) and (kmax == kmin + 1):
        buffer[jmin][kmin] += w * q
        buffer_weight[jmin][kmin] += w
        return

    depth = int(ceil(log2(dx * float(max(Nx, Ny)))))
    if depth < 1:
        depth = 1
    if depth >= max_depth:
        depth = max_depth - 1

    jmax = min(Nsx, 1 << depth)
    kmax = min(Nsy, 1 << depth)
    for j in range(jmax):
        xx = ((j + 0.5 - jmax / 2) / jmax * cell_max_width + x) * Nx
        jj = int(xx)
        if jj < 0 or jj >= Nx: continue

        dvx = 1 - (xx - jj)
        for k in range(kmax):
            if stamp_mask[depth, j, k] == 0:
                continue

            yy = ((k + 0.5 - kmax / 2) / kmax * cell_max_width + y) * Ny
            kk = int(yy)
            if kk < 0 or kk >= Ny: continue

            dvy = 1 - (yy - kk)

            swq = stamp[depth, j, k] * w * q
            sw =  stamp[depth, j, k] * w

            buffer[jj][kk]        += swq * (dvx * dvy)
            buffer_weight[jj][kk] +=  sw * (dvx * dvy)

            if jj < Nx - 1:
                buffer[jj + 1][kk]        += swq * (1 - dvx) * dvy
                buffer_weight[jj + 1][kk] +=  sw * (1 - dvx) * dvy

            if kk < Ny - 1:
                buffer[jj][kk + 1]        += swq * dvx * (1 - dvy)
                buffer_weight[jj][kk + 1] +=  sw * dvx * (1 - dvy)
            if jj < Nx - 1 and kk < Ny - 1:
                buffer[jj + 1][kk + 1]        += swq * (1 - dvx) * (1 - dvy)
                buffer_weight[jj + 1][kk + 1] +=  sw * (1 - dvx) * (1 - dvy)


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def add_cells_to_image_offaxis(
    *,
    const np.float64_t[:, ::1] Xp,
    const np.float64_t[::1] dXp,
    const np.float64_t[::1] qty,
    const np.float64_t[::1] weight,
    const float theta,
    const float phi,
    np.float64_t[:, ::1] buffer,
    np.float64_t[:, ::1] buffer_weight,
    const int Nx,
    const int Ny,
    const int Npix_min = 4,
):
    from scipy.spatial.transform import Rotation

    cdef np.ndarray[np.float64_t, ndim=1] center = np.array([0.5, 0.5, 0.5])
    cdef np.float64_t w0 = 1 / sqrt(3.)
    cdef int i, j, k

    cdef np.ndarray[np.float64_t, ndim=1] a = np.array([1., 0, 0]) * w0
    cdef np.ndarray[np.float64_t, ndim=1] b = np.array([0, 1., 0]) * w0
    cdef np.ndarray[np.float64_t, ndim=1] c = np.array([0, 0, 1.]) * w0
    cdef np.ndarray[np.float64_t, ndim=2] R = Rotation.from_euler('ZYX', [phi, theta, 0]).as_matrix()

    a = np.dot(R, a)
    b = np.dot(R, b)
    c = np.dot(R, c)

    cdef np.ndarray[np.float64_t, ndim=1] o = center - (a + b + c) / 2

    cdef int Nsx, Nsy
    cdef np.float64_t dx_max = np.max(dXp)
    # The largest cell needs to be resolved by at least this number of pixels
    Nsx = max(Npix_min, int(ceil(dx_max * float(Nx))))
    Nsy = max(Npix_min, int(ceil(dx_max * float(Ny))))
    cdef int max_depth = int(ceil(log2(max(Nsx, Nsy))))
    cdef int depth

    print(f"{Nsx=} {Nsy=}")

    cdef np.ndarray[np.float64_t, ndim=3] stamp_arr = np.zeros((max_depth, Nsx, Nsy), dtype=float)
    cdef np.ndarray[np.uint8_t, ndim=3] stamp_mask_arr = np.zeros((max_depth, Nsx, Nsy), dtype=np.uint8)
    cdef np.float64_t[:, :, ::1] stamp = stamp_arr
    cdef np.uint8_t[:, :, ::1] stamp_mask = stamp_mask_arr

    cdef np.float64_t[:, ::1] Xrot
    Xp_center = np.array([0.5, 0.5, 0.5])
    Xrot = np.dot(R, (Xp - Xp_center).T).T.copy() + Xp_center

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
            min(1 << depth, Nsx),
            min(1 << depth, Nsy),
        )

        stamp_arr[depth] /= np.sum(stamp[depth])

    # Iterate over all cells, applying the stamp
    cdef np.float64_t x, y, dx

    cdef np.float64_t w, q, cell_max_width, sq3

    sq3 = sqrt(3.)

    # Local buffers
    cdef np.float64_t **lbuffer
    cdef np.float64_t **lbuffer_weight

    cdef int num_particles = len(Xp)

    with nogil, parallel():
        lbuffer = <np.float64_t**> malloc(sizeof(np.float64_t*) * Nx)
        lbuffer_weight = <np.float64_t**> malloc(sizeof(np.float64_t*) * Nx)
        for j in range(Nx):
            lbuffer[j] = <np.float64_t*> malloc(sizeof(np.float64_t) * Ny)
            lbuffer_weight[j] = <np.float64_t*> malloc(sizeof(np.float64_t) * Ny)
            for k in range(Ny):
                lbuffer[j][k] = 0
                lbuffer_weight[j][k] = 0

        for i in prange(num_particles, schedule="static"):
            dx = dXp[i]
            w = weight[i]
            q = qty[i]
            cell_max_width = dx * sq3
            x = Xrot[i, 0]
            y = Xrot[i, 1]

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
                    buffer[j, k] += lbuffer[j][k]
                    buffer_weight[j, k] += lbuffer_weight[j][k]
        # Free memory
        for j in range(Nx):
            free(lbuffer[j])
            free(lbuffer_weight[j])
        free(lbuffer)
        free(lbuffer_weight)


@cython.boundscheck(False)
cdef inline np.float64_t det2d(const np.float64_t[::1] a, const np.float64_t[::1] b) noexcept nogil:
    return a[0] * b[1] - a[1] * b[0]

@cython.cdivision(True)
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
    cdef np.float64_t[::1] OfrontA = np.zeros(2), ObackA = np.zeros(2)

    cdef np.float64_t inv_dx = 1 / float(Nx)
    cdef np.float64_t inv_dy = 1 / float(Ny)
    cdef np.float64_t[2] nm
    cdef bint within
    cdef np.float64_t[::1] all_z = np.empty(6)
    cdef np.float64_t zmin, zmax, z
    cdef int Nhit, i, j
    for i in range(Nx):
        X[0] = (i + 0.5) * inv_dx
        for j in range(Ny):
            zmin = np.inf
            zmax = -np.inf
            Nhit = 0
            X[1] = (j + 0.5) * inv_dy

            OfrontA[0] = X[0] - O[0]
            OfrontA[1] = X[1] - O[1]
            all_z[:] = 0
            within = check_in_parallelogram(OfrontA, v2d, u2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * u[2] + nm[1] * v[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(OfrontA, w2d, v2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * v[2] + nm[1] * w[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(OfrontA, w2d, u2d, 1, 1, nm)
            if within:
                z = O[2] + nm[0] * u[2] + nm[1] * w[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            ObackA[0] = X[0] - Oback[0]
            ObackA[1] = X[1] - Oback[1]
            within = check_in_parallelogram(ObackA, v2d, u2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * u[2] + nm[1] * v[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(ObackA, w2d, v2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * v[2] + nm[1] * w[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            within = check_in_parallelogram(ObackA, w2d, u2d, -1, -1, nm)
            if within:
                z = Oback[2] + nm[0] * u[2] + nm[1] * w[2]
                zmin = min(z, zmin)
                zmax = max(z, zmax)
                Nhit += 1

            if Nhit == 0:
                continue
            elif Nhit == 1:
                raise RuntimeError("This should not happen")
            else:
                buffer[i, j] += zmax - zmin
                buffer_mask[i, j] |= 1

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
