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
    np.float64_t[:, :] buffer,
    np.uint8_t[:, :] buffer_mask,
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

    cdef np.float64_t[:, :] buffer_view = buffer
    cdef np.uint8_t[:, :] buffer_mask_view = buffer_mask

    for pi in range(npx):
        _add_point_to_greyscale_image(buffer_view, buffer_mask_view, px[pi], py[pi], pv[pi], xs, ys)
    return


@cython.boundscheck(False)
def add_points_to_greyscale_image_with_lagrangian_tesselation(
    np.ndarray[np.float64_t, ndim=2] buffer,
    np.ndarray[np.uint8_t,   ndim=2] buffer_mask,
    np.ndarray[np.float64_t, ndim=4] p3d,
    np.ndarray[np.float64_t, ndim=3] pv,
    int split=100,
):
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    cdef int Ngrid = p3d.shape[0] - 1

    cdef np.float64_t[:, :] buffer_view = buffer
    cdef np.uint8_t[:, :] buffer_mask_view = buffer_mask

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

    cdef np.float64_t[3] orig
    cdef np.float64_t[3] a
    cdef np.float64_t[3] b
    cdef np.float64_t[3] c

    cdef np.float64_t[3] newp

    cdef int i, j, k, m, s, idim

    cdef np.float64_t f = 1.0 / split / 6

    cdef np.uint8_t[:, :] off
    cdef np.float64_t[:, :] ro

    for m in range(len(conn)):
        ro = np.random.random(size=(split, 3)) # random offsets
        off = vert[conn[m]]

        for i in range(Ngrid):
            for j in range(Ngrid):
                for k in range(Ngrid):
                    for idim in range(2):   # stop at 2, because z is not used
                        orig[idim] = p3d[i + off[3][0], j + off[3][1], k + off[3][2], idim]
                        a[idim] = _wrap_dist(p3d[i + off[0][0], j + off[0][1], k + off[0][2], idim] - orig[idim])
                        b[idim] = _wrap_dist(p3d[i + off[1][0], j + off[1][1], k + off[1][2], idim] - orig[idim])
                        c[idim] = _wrap_dist(p3d[i + off[2][0], j + off[2][1], k + off[2][2], idim] - orig[idim])

                    for s in range(split):
                        for idim in range(2):  # stop at 2, because z is not used
                            newp[idim] = orig[idim] + (
                                ro[s, 0] * a[idim] +
                                ro[s, 1] * b[idim] +
                                ro[s, 2] * c[idim]
                            )
                            newp[idim] %= 1.0

                        _add_point_to_greyscale_image(
                            buffer_view, buffer_mask_view, newp[0], newp[1], f, xs, ys
                        )

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
