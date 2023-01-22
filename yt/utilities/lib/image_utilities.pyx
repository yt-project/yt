# distutils: libraries = STD_LIBS
"""
Utilities for images
"""


import numpy as np

cimport cython
cimport numpy as np
from libc.math cimport M_PI, M_PI_2, atan2, ceil, floor

from yt.utilities.lib.fp_utils cimport fabs, iclip


@cython.boundscheck(False)
cdef inline void _add_point_to_greyscale_image(
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    np.float64_t px,
    np.float64_t py,
    np.float64_t pv,
    int xs,
    int ys
) nogil:
    cdef int i, j
    j = <int> (xs * px)
    i = <int> (ys * py)
    buffer[i, j] += pv
    buffer_mask[i, j] = 1


cdef inline np.float64_t  _wrap_dist(np.float64_t dx) nogil:
    if dx > 0.5:
        return 1 - dx
    elif dx < -0.5:
        return 1 + dx
    else:
        return dx

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

cdef np.float64_t twopi = 2 * M_PI
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
    cdef int i, imin
    imin = 0
    for i in range(1, 4):
        if ABCD[i][0] < ABCD[imin][0]:
            imin = i
    hull[0] = imin
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
cdef void integrate_quad(
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

            Vtot += v

    # Special case: Vtot == 0
    if Vtot == 0:
        i = <int> (xs * (A[0] + B[0] + C[0]) / 3) % xs
        j = <int> (ys * (A[1] + B[1] + C[1]) / 3) % ys
        buffer[j, i] += value
        buffer_mask[j, i] = 1
        return

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
cdef integrate_tetrahedron_proj(
    np.float64_t[::1] origin,
    np.float64_t[::1] a,
    np.float64_t[::1] b,
    np.float64_t[::1] c,
    np.float64_t value,
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    int xs,
    int ys
) nogil:
    # Shortcut: if the tetrahedron is smaller than a pixel
    # then just deposit the value in the center of the pixel
    cdef np.float64_t Vtetra = 1 / 6 * fabs(
        a[0] * (b[1] * c[2] - b[2] * c[1]) +
        a[1] * (b[2] * c[0] - b[0] * c[2]) +
        a[2] * (b[0] * c[1] - b[1] * c[0])
    )

    if Vtetra**2 < 1 / (xs * ys)**3:
        i = <int> (xs * (origin[0] + a[0] + b[0] + c[0]) / 4) % xs
        j = <int> (ys * (origin[1] + a[1] + b[1] + c[1]) / 4) % ys
        buffer[j, i] += value
        buffer_mask[j, i] = 1
        return
    # Find location of four vertices
    # cdef np.ndarray[np.float64_t, ndim=1] A = origin
    # cdef np.ndarray[np.float64_t, ndim=1] B = origin + a
    # cdef np.ndarray[np.float64_t, ndim=1] C = origin + b
    # cdef np.ndarray[np.float64_t, ndim=1] D = origin + c

    # Project them onto the xy plane
    cdef np.float64_t[4][2] ABCD_xy
    cdef np.float64_t[2] A_xy, B_xy, C_xy, D_xy, P_xy, P1_xy, P2_xy
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
    # print(iABC[0], iABC[1], iABC[2], iABC[3])
    cdef np.float64_t det, S1, S2
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

    # cdef np.ndarray[np.float64_t, ndim=2] sub_o = np.array((
    #     (   0,    0,    0),
    #     ( 0.5,    0,    0),
    #     (   0,  0.5,    0),
    #     (   0,    0,  0.5),
    #     ( 0.5,    0,    0), # ABCD
    #     ( 0.5,    0,    0), # ACDE
    #     (   0,  0.5,    0), # CDEF
    #     (   0,  0.5,    0), # CBDF
    # ))
    # cdef np.ndarray[np.float64_t, ndim=2] sub_a = np.array((
    #     ( 0.5,    0,    0),
    #     ( 0.5,    0,    0),
    #     ( 0.5,    0,    0),
    #     ( 0.5,    0,    0),
    #     (   0,  0.5,    0), # AB
    #     (-0.5,  0.5,    0), # AC
    #     ( 0.5, -0.5,  0.5), # CD
    #     ( 0.5,    0,    0), # CB
    # ))
    # cdef np.ndarray[np.float64_t, ndim=2] sub_b = np.array((
    #     (   0,  0.5,    0),
    #     (   0,  0.5,    0),
    #     (   0,  0.5,    0),
    #     (   0,  0.5,    0),
    #     (-0.5,  0.5,    0), # AC
    #     (   0,    0,  0.5), # AD
    #     (   0, -0.5,  0.5), # CE
    #     ( 0.5, -0.5,  0.5), # CD
    # ))
    # cdef np.ndarray[np.float64_t, ndim=2] sub_c = np.array((
    #     (   0,    0,  0.5),
    #     (   0,    0,  0.5),
    #     (   0,    0,  0.5),
    #     (   0,    0,  0.5),
    #     (   0,    0,  0.5), # AD
    #     (-0.5,    0,  0.5), # AE
    #     (   0,    0,  0.5), # CF
    #     (   0,    0,  0.5), # CF
    # ))

    cdef np.float64_t[3] orig
    cdef np.float64_t[3] a
    cdef np.float64_t[3] b
    cdef np.float64_t[3] c

    # cdef np.float64_t[3] newp

    cdef int i, j, k, m, idim

    # cdef np.float64_t f = 1.0 / 6.0
    # cdef np.float64_t Vtetra

    cdef np.uint8_t[:, ::1] off

    # cdef np.float64_t dV = 1 / (xs * ys * np.sqrt(xs * ys))

    # cdef np.float64_t* sub_o_ptr = &sub_o[0, 0]
    # cdef np.float64_t* sub_a_ptr = &sub_a[0, 0]
    # cdef np.float64_t* sub_b_ptr = &sub_b[0, 0]
    # cdef np.float64_t* sub_c_ptr = &sub_c[0, 0]

    for m in range(len(conn)):
        # ro = sample_tetrahedron(split)
        off = vert[conn[m]]

        for i in range(Ngrid):
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
                    # Vtetra = 1 / 6 * abs(
                    #     a[0] * (b[1] * c[2] - b[2] * c[1]) +
                    #     a[1] * (b[2] * c[0] - b[0] * c[2]) +
                    #     a[2] * (b[0] * c[1] - b[1] * c[0])
                    # )
                    # recursive_tetrahedron(
                    #     orig,
                    #     a,
                    #     b,
                    #     c,
                    #     dV,
                    #     f,
                    #     Vtetra,
                    #     buffer_view,
                    #     buffer_mask_view,
                    #     xs,
                    #     ys,
                    #     sub_o,
                    #     sub_a,
                    #     sub_b,
                    #     sub_c,
                    # )

    # print("Done!")

def tetra_helper(
    np.float64_t[:] orig,
    np.float64_t[:] a,
    np.float64_t[:] b,
    np.float64_t[:] c,
    np.float64_t dV,
    np.float64_t f,
    np.float64_t[:, ::1] buffer,
    np.uint8_t[:, ::1] buffer_mask,
    int xs,
    int ys,
):
    cdef np.float64_t Vtetra = 1 / 6 * abs(
        a[0] * (b[1] * c[2] - b[2] * c[1]) +
        a[1] * (b[2] * c[0] - b[0] * c[2]) +
        a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    cdef np.ndarray[np.float64_t, ndim=2] sub_o = np.array((
        (   0,    0,    0),
        ( 0.5,    0,    0),
        (   0,  0.5,    0),
        (   0,    0,  0.5),
        ( 0.5,    0,    0), # ABCD
        ( 0.5,    0,    0), # ACDE
        (   0,  0.5,    0), # CDEF
        (   0,  0.5,    0), # CBDF
    ))
    cdef np.ndarray[np.float64_t, ndim=2] sub_a = np.array((
        ( 0.5,    0,    0),
        ( 0.5,    0,    0),
        ( 0.5,    0,    0),
        ( 0.5,    0,    0),
        (   0,  0.5,    0), # AB
        (-0.5,  0.5,    0), # AC
        ( 0.5, -0.5,  0.5), # CD
        ( 0.5,    0,    0), # CB
    ))
    cdef np.ndarray[np.float64_t, ndim=2] sub_b = np.array((
        (   0,  0.5,    0),
        (   0,  0.5,    0),
        (   0,  0.5,    0),
        (   0,  0.5,    0),
        (-0.5,  0.5,    0), # AC
        (   0,    0,  0.5), # AD
        (   0, -0.5,  0.5), # CE
        ( 0.5, -0.5,  0.5), # CD
    ))
    cdef np.ndarray[np.float64_t, ndim=2] sub_c = np.array((
        (   0,    0,  0.5),
        (   0,    0,  0.5),
        (   0,    0,  0.5),
        (   0,    0,  0.5),
        (   0,    0,  0.5), # AD
        (-0.5,    0,  0.5), # AE
        (   0,    0,  0.5), # CF
        (   0,    0,  0.5), # CF
    ))
    recursive_tetrahedron(
        &orig[0],
        &a[0],
        &b[0],
        &c[0],
        dV,
        f,
        Vtetra,
        buffer,
        buffer_mask,
        xs,
        ys,
        sub_o,
        sub_a,
        sub_b,
        sub_c,
    )

@cython.boundscheck(False)
cdef void recursive_tetrahedron(
    np.float64_t* orig,
    np.float64_t* a,
    np.float64_t* b,
    np.float64_t* c,
    np.float64_t dV,
    np.float64_t f,
    np.float64_t Vtetra,
    np.float64_t[:, ::1] buffer_view,
    np.uint8_t[:, ::1] buffer_mask_view,
    int xs,
    int ys,
    np.float64_t[:, ::1] sub_o,
    np.float64_t[:, ::1] sub_a,
    np.float64_t[:, ::1] sub_b,
    np.float64_t[:, ::1] sub_c,
) nogil:
    # V = a · |b x c| / 6
    cdef int isub
    cdef np.float64_t[3] com, o2, a2, b2, c2

    if Vtetra < dV:
        # Add center of mass to buffer
        com[0] = orig[0] + 0.25 * (a[0] + b[0] + c[0])
        com[1] = orig[1] + 0.25 * (a[1] + b[1] + c[1])
        com[2] = orig[2] + 0.25 * (a[2] + b[2] + c[2])
        _add_point_to_greyscale_image(buffer_view, buffer_mask_view, com[0], com[1], f, xs, ys)

        return

    for isub in range(8):
        for idim in range(3):
            o2[idim] = sub_o[isub, 0] * a[idim] + sub_o[isub, 1] * b[idim] + sub_o[isub, 2] * c[idim] + orig[idim]
            a2[idim] = sub_a[isub, 0] * a[idim] + sub_a[isub, 1] * b[idim] + sub_a[isub, 2] * c[idim]
            b2[idim] = sub_b[isub, 0] * a[idim] + sub_b[isub, 1] * b[idim] + sub_b[isub, 2] * c[idim]
            c2[idim] = sub_c[isub, 0] * a[idim] + sub_c[isub, 1] * b[idim] + sub_c[isub, 2] * c[idim]

        recursive_tetrahedron(o2, a2, b2, c2, dV, f / 8, Vtetra / 8, buffer_view, buffer_mask_view, xs, ys, sub_o, sub_a, sub_b, sub_c)



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
