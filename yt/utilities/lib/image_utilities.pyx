"""
Utilities for images
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
from yt.utilities.lib.fp_utils cimport iclip

def add_points_to_greyscale_image(
        np.ndarray[np.float64_t, ndim=2] buffer,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.ndarray[np.float64_t, ndim=1] pv):
    cdef int i, j, pi
    cdef int np = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    for pi in range(np):
        j = <int> (xs * px[pi])
        i = <int> (ys * py[pi])
        buffer[i, j] += pv[pi]
    return

def add_points_to_image(
        np.ndarray[np.uint8_t, ndim=3] buffer,
        np.ndarray[np.float64_t, ndim=1] px,
        np.ndarray[np.float64_t, ndim=1] py,
        np.float64_t pv):
    cdef int i, j, k, pi
    cdef int np = px.shape[0]
    cdef int xs = buffer.shape[0]
    cdef int ys = buffer.shape[1]
    cdef int v
    v = iclip(<int>(pv * 255), 0, 255)
    for pi in range(np):
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
