"""
Functions for computing the extent of lenses and whatnot



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
from .image_samplers cimport ImageContainer

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int calculate_extent_plane_parallel(ImageContainer *image,
            VolumeContainer *vc, np.int64_t rv[4]) nogil except -1:
    # We do this for all eight corners
    cdef np.float64_t temp
    cdef np.float64_t *edges[2]
    cdef np.float64_t cx, cy
    cdef np.float64_t extrema[4]
    cdef int i, j, k
    edges[0] = vc.left_edge
    edges[1] = vc.right_edge
    extrema[0] = extrema[2] = 1e300; extrema[1] = extrema[3] = -1e300
    for i in range(2):
        for j in range(2):
            for k in range(2):
                # This should rotate it into the vector plane
                temp  = edges[i][0] * image.x_vec[0]
                temp += edges[j][1] * image.x_vec[1]
                temp += edges[k][2] * image.x_vec[2]
                if temp < extrema[0]: extrema[0] = temp
                if temp > extrema[1]: extrema[1] = temp
                temp  = edges[i][0] * image.y_vec[0]
                temp += edges[j][1] * image.y_vec[1]
                temp += edges[k][2] * image.y_vec[2]
                if temp < extrema[2]: extrema[2] = temp
                if temp > extrema[3]: extrema[3] = temp
    cx = cy = 0.0
    for i in range(3):
        cx += image.center[i] * image.x_vec[i]
        cy += image.center[i] * image.y_vec[i]
    rv[0] = lrint((extrema[0] - cx - image.bounds[0])/image.pdx)
    rv[1] = rv[0] + lrint((extrema[1] - extrema[0])/image.pdx)
    rv[2] = lrint((extrema[2] - cy - image.bounds[2])/image.pdy)
    rv[3] = rv[2] + lrint((extrema[3] - extrema[2])/image.pdy)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int calculate_extent_perspective(ImageContainer *image,
            VolumeContainer *vc, np.int64_t rv[4]) nogil except -1:

    cdef np.float64_t cam_pos[3]
    cdef np.float64_t cam_width[3]
    cdef np.float64_t north_vector[3]
    cdef np.float64_t east_vector[3]
    cdef np.float64_t normal_vector[3]
    cdef np.float64_t vertex[3]
    cdef np.float64_t pos1[3]
    cdef np.float64_t sight_vector[3]
    cdef np.float64_t sight_center[3]
    cdef np.float64_t corners[3][8]
    cdef float sight_vector_norm, sight_angle_cos, sight_length, dx, dy
    cdef int i, iv, px, py
    cdef int min_px, min_py, max_px, max_py

    min_px = SHRT_MAX
    min_py = SHRT_MAX
    max_px = -SHRT_MAX
    max_py = -SHRT_MAX

    # calculate vertices for 8 corners of vc
    corners[0][0] = vc.left_edge[0]
    corners[0][1] = vc.right_edge[0]
    corners[0][2] = vc.right_edge[0]
    corners[0][3] = vc.left_edge[0]
    corners[0][4] = vc.left_edge[0]
    corners[0][5] = vc.right_edge[0]
    corners[0][6] = vc.right_edge[0]
    corners[0][7] = vc.left_edge[0]

    corners[1][0] = vc.left_edge[1]
    corners[1][1] = vc.left_edge[1]
    corners[1][2] = vc.right_edge[1]
    corners[1][3] = vc.right_edge[1]
    corners[1][4] = vc.left_edge[1]
    corners[1][5] = vc.left_edge[1]
    corners[1][6] = vc.right_edge[1]
    corners[1][7] = vc.right_edge[1]

    corners[2][0] = vc.left_edge[2]
    corners[2][1] = vc.left_edge[2]
    corners[2][2] = vc.left_edge[2]
    corners[2][3] = vc.left_edge[2]
    corners[2][4] = vc.right_edge[2]
    corners[2][5] = vc.right_edge[2]
    corners[2][6] = vc.right_edge[2]
    corners[2][7] = vc.right_edge[2]

    # This code was ported from
    #   yt.visualization.volume_rendering.lens.PerspectiveLens.project_to_plane()
    for i in range(3):
        cam_pos[i] = image.camera_data[0, i]
        cam_width[i] = image.camera_data[1, i]
        east_vector[i] = image.camera_data[2, i]
        north_vector[i] = image.camera_data[3, i]
        normal_vector[i] = image.camera_data[4, i]

    for iv in range(8):
        vertex[0] = corners[0][iv]
        vertex[1] = corners[1][iv]
        vertex[2] = corners[2][iv]

        cam_width[1] = cam_width[0] * image.nv[1] / image.nv[0]

        subtract(vertex, cam_pos, sight_vector)
        fma(cam_width[2], normal_vector, cam_pos, sight_center)

        sight_vector_norm = L2_norm(sight_vector)
       
        if sight_vector_norm != 0:
            for i in range(3):
                sight_vector[i] /= sight_vector_norm

        sight_angle_cos = dot(sight_vector, normal_vector)
        sight_angle_cos = fclip(sight_angle_cos, -1.0, 1.0)

        if acos(sight_angle_cos) < 0.5 * M_PI and sight_angle_cos != 0.0:
            sight_length = cam_width[2] / sight_angle_cos
        else:
            sight_length = sqrt(cam_width[0]**2 + cam_width[1]**2)
            sight_length = sight_length / sqrt(1.0 - sight_angle_cos**2)

        fma(sight_length, sight_vector, cam_pos, pos1)
        subtract(pos1, sight_center, pos1)
        dx = dot(pos1, east_vector)
        dy = dot(pos1, north_vector)

        px = int(image.nv[0] * 0.5 + image.nv[0] / cam_width[0] * dx)
        py = int(image.nv[1] * 0.5 + image.nv[1] / cam_width[1] * dy)
        min_px = min(min_px, px)
        max_px = max(max_px, px)
        min_py = min(min_py, py)
        max_py = max(max_py, py)

    rv[0] = max(min_px, 0)
    rv[1] = min(max_px, image.nv[0])
    rv[2] = max(min_py, 0)
    rv[3] = min(max_py, image.nv[1])
    return 0

# We do this for a bunch of lenses.  Fallback is to grab them from the vector
# info supplied.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int calculate_extent_null(ImageContainer *image,
            VolumeContainer *vc, np.int64_t rv[4]) nogil except -1:
    rv[0] = 0
    rv[1] = image.nv[0]
    rv[2] = 0
    rv[3] = image.nv[1]
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void generate_vector_info_plane_parallel(ImageContainer *im,
            np.int64_t vi, np.int64_t vj,
            np.float64_t width[2],
            # Now outbound
            np.float64_t v_dir[3], np.float64_t v_pos[3]) nogil:
    cdef int i
    cdef np.float64_t px, py
    px = width[0] * (<np.float64_t>vi)/(<np.float64_t>im.nv[0]-1) - width[0]/2.0
    py = width[1] * (<np.float64_t>vj)/(<np.float64_t>im.nv[1]-1) - width[1]/2.0
    # atleast_3d will add to beginning and end
    v_pos[0] = im.vp_pos[0,0,0]*px + im.vp_pos[0,3,0]*py + im.vp_pos[0,9,0]
    v_pos[1] = im.vp_pos[0,1,0]*px + im.vp_pos[0,4,0]*py + im.vp_pos[0,10,0]
    v_pos[2] = im.vp_pos[0,2,0]*px + im.vp_pos[0,5,0]*py + im.vp_pos[0,11,0]
    for i in range(3): v_dir[i] = im.vp_dir[0,i,0]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void generate_vector_info_null(ImageContainer *im,
            np.int64_t vi, np.int64_t vj,
            np.float64_t width[2],
            # Now outbound
            np.float64_t v_dir[3], np.float64_t v_pos[3]) nogil:
    cdef int i
    for i in range(3):
        # Here's a funny thing: we use vi here because our *image* will be
        # flattened.  That means that im.nv will be a better one-d offset,
        # since vp_pos has funny strides.
        v_pos[i] = im.vp_pos[vi, vj, i]
        v_dir[i] = im.vp_dir[vi, vj, i]

