# distutils: libraries = STD_LIBS
"""
Utilities for images
"""


import numpy as np

cimport cython
cimport numpy as np

from yt.utilities.lib.fp_utils cimport iclip


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


cdef inline np.float64_t  _wrap_dist(np.float64_t dx):
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

@cython.boundscheck(False)
def add_points_to_greyscale_image_with_lagrangian_tesselation(
    np.ndarray[np.float64_t, ndim=2] buffer,
    np.ndarray[np.uint8_t,   ndim=2] buffer_mask,
    np.ndarray[np.float64_t, ndim=4] p3d,
    np.ndarray[np.float64_t, ndim=3] pv,
):
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    cdef int Ngrid = p3d.shape[0] - 1

    cdef np.float64_t[:, ::1] buffer_view = buffer
    cdef np.uint8_t[:, ::1] buffer_mask_view = buffer_mask

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

    cdef np.float64_t[3] orig
    cdef np.float64_t[3] a
    cdef np.float64_t[3] b
    cdef np.float64_t[3] c

    # cdef np.float64_t[3] newp

    cdef int i, j, k, m, idim

    cdef np.float64_t f = 1.0 / 6.0
    cdef np.float64_t Vtetra

    cdef np.ndarray[np.uint8_t, ndim=2] off

    cdef np.float64_t dV = 1 / (xs * ys * np.sqrt(xs * ys))

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
                        orig[idim] = p3d[i + off[3][0], j + off[3][1], k + off[3][2], idim]
                        a[idim] = _wrap_dist(p3d[i + off[0][0], j + off[0][1], k + off[0][2], idim] - orig[idim])
                        b[idim] = _wrap_dist(p3d[i + off[1][0], j + off[1][1], k + off[1][2], idim] - orig[idim])
                        c[idim] = _wrap_dist(p3d[i + off[2][0], j + off[2][1], k + off[2][2], idim] - orig[idim])

                    Vtetra = 1 / 6 * abs(
                        a[0] * (b[1] * c[2] - b[2] * c[1]) +
                        a[1] * (b[2] * c[0] - b[0] * c[2]) +
                        a[2] * (b[0] * c[1] - b[1] * c[0])
                    )
                    recursive_tetrahedron(
                        orig,
                        a,
                        b,
                        c,
                        dV,
                        f,
                        Vtetra,
                        buffer_view,
                        buffer_mask_view,
                        xs,
                        ys,
                        sub_o,
                        sub_a,
                        sub_b,
                        sub_c,
                    )

    print("Done!")

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
