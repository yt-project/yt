"""
Checks for points contained in a volume



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

cdef extern from "math.h":
    double fabs(double x)

@cython.wraparound(False)
@cython.boundscheck(False)
def planar_points_in_volume(
                   np.ndarray[np.float64_t, ndim=2] points,
                   np.ndarray[np.int8_t, ndim=1] pmask,  # pixel mask
                   np.ndarray[np.float64_t, ndim=1] left_edge,
                   np.ndarray[np.float64_t, ndim=1] right_edge,
                   np.ndarray[np.int32_t, ndim=3] mask,
                   float dx):
    cdef np.ndarray[np.int8_t, ndim=1] \
         valid = np.zeros(points.shape[0], dtype='int8')
    cdef int i, dim, count
    cdef int ex
    cdef double dx_inv
    cdef unsigned int idx[3]
    count = 0
    dx_inv = 1.0 / dx
    for i in xrange(points.shape[0]):
        if pmask[i] == 0:
            continue
        ex = 1
        for dim in xrange(3):
            if points[i,dim] < left_edge[dim] or points[i,dim] > right_edge[dim]:
                valid[i] = ex = 0
                break
        if ex == 1:
            for dim in xrange(3):
                idx[dim] = <unsigned int> \
                           ((points[i,dim] - left_edge[dim]) * dx_inv)
            if mask[idx[0], idx[1], idx[2]] == 1:
                valid[i] = 1
                count += 1

    cdef np.ndarray[np.int32_t, ndim=1] result = np.empty(count, dtype='int32')
    count = 0
    for i in xrange(points.shape[0]):
        if valid[i] == 1 and pmask[i] == 1:
            result[count] = i
            count += 1

    return result

cdef inline void set_rotated_pos(
            np.float64_t cp[3], np.float64_t rdds[3][3],
            np.float64_t rorigin[3], int i, int j, int k):
    cdef int oi
    for oi in range(3):
        cp[oi] = rdds[0][oi] * (0.5 + i) \
               + rdds[1][oi] * (0.5 + j) \
               + rdds[2][oi] * (0.5 + k) \
               + rorigin[oi]

#@cython.wraparound(False)
#@cython.boundscheck(False)
def grid_points_in_volume(
                   np.ndarray[np.float64_t, ndim=1] box_lengths,
                   np.ndarray[np.float64_t, ndim=1] box_origin,
                   np.ndarray[np.float64_t, ndim=2] rot_mat,
                   np.ndarray[np.float64_t, ndim=1] grid_left_edge,
                   np.ndarray[np.float64_t, ndim=1] grid_right_edge,
                   np.ndarray[np.float64_t, ndim=1] dds,
                   np.ndarray[np.int32_t, ndim=3] mask,
                   int break_first):
    cdef int n[3]
    cdef i, j, k
    cdef np.float64_t rds[3][3]
    cdef np.float64_t cur_pos[3]
    cdef np.float64_t rorigin[3]
    for i in range(3):
        rorigin[i] = 0.0
    for i in range(3):
        n[i] = mask.shape[i]
        for j in range(3):
            # Set up our transposed dx, which has a component in every
            # direction
            rds[i][j] = dds[i] * rot_mat[j,i]
            # In our rotated coordinate system, the box origin is 0,0,0
            # so we subtract the box_origin from the grid_origin and rotate
            # that
            rorigin[j] += (grid_left_edge[i] - box_origin[i]) * rot_mat[j,i]

    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                set_rotated_pos(cur_pos, rds, rorigin, i, j, k)
                if (cur_pos[0] > box_lengths[0]): continue
                if (cur_pos[1] > box_lengths[1]): continue
                if (cur_pos[2] > box_lengths[2]): continue
                if (cur_pos[0] < 0.0): continue
                if (cur_pos[1] < 0.0): continue
                if (cur_pos[2] < 0.0): continue
                if break_first:
                    if mask[i,j,k]: return 1
                else:
                    mask[i,j,k] = 1
    return 0

cdef void normalize_vector(np.float64_t vec[3]):
    cdef int i
    cdef np.float64_t norm = 0.0
    for i in range(3):
        norm += vec[i]*vec[i]
    norm = norm**0.5
    for i in range(3):
        vec[i] /= norm

cdef void get_cross_product(np.float64_t v1[3],
                            np.float64_t v2[3],
                            np.float64_t cp[3]):
    cp[0] = v1[1]*v2[2] - v1[2]*v2[1]
    cp[1] = v1[3]*v2[0] - v1[0]*v2[3]
    cp[2] = v1[0]*v2[1] - v1[1]*v2[0]
    #print cp[0], cp[1], cp[2]

cdef int check_projected_overlap(
        np.float64_t sep_ax[3], np.float64_t sep_vec[3], int gi,
        np.float64_t b_vec[3][3], np.float64_t g_vec[3][3]):
    cdef int g_ax, b_ax
    cdef np.float64_t tba, tga, ba, ga, sep_dot
    ba = ga = sep_dot = 0.0
    for g_ax in range(3):
        # We need the grid vectors, which we'll precompute here
        tba = tga = 0.0
        for b_ax in range(3):
            tba += b_vec[g_ax][b_ax] * sep_vec[b_ax]
            tga += g_vec[g_ax][b_ax] * sep_vec[b_ax]
        ba += fabs(tba)
        ga += fabs(tga)
        sep_dot += sep_vec[g_ax] * sep_ax[g_ax]
    #print sep_vec[0], sep_vec[1], sep_vec[2],
    #print sep_ax[0], sep_ax[1], sep_ax[2]
    return (fabs(sep_dot) > ba+ga)
    # Now we do

@cython.wraparound(False)
@cython.boundscheck(False)
def find_grids_in_inclined_box(
        np.ndarray[np.float64_t, ndim=2] box_vectors,
        np.ndarray[np.float64_t, ndim=1] box_center,
        np.ndarray[np.float64_t, ndim=2] grid_left_edges,
        np.ndarray[np.float64_t, ndim=2] grid_right_edges):

    # http://www.gamasutra.com/view/feature/3383/simple_intersection_tests_for_games.php?page=5
    cdef int n = grid_right_edges.shape[0]
    cdef int g_ax, b_ax, gi
    cdef np.float64_t b_vec[3][3]
    cdef np.float64_t g_vec[3][3]
    cdef np.float64_t a_vec[3][3]
    cdef np.float64_t sep_ax[15][3]
    cdef np.float64_t sep_vec[3]
    cdef np.ndarray[np.int32_t, ndim=1] good = np.zeros(n, dtype='int32')
    cdef np.ndarray[np.float64_t, ndim=2] grid_centers
    # Fill in our axis unit vectors
    for b_ax in range(3):
        for g_ax in range(3):
            a_vec[b_ax][g_ax] = <np.float64_t> (b_ax == g_ax)
    grid_centers = (grid_right_edges + grid_left_edges)/2.0

    # Now we pre-compute our candidate separating axes, because the unit
    # vectors for all the grids are identical
    for b_ax in range(3):
        # We have 6 principal axes we already know, which are the grid (domain)
        # principal axes and the box axes
        sep_ax[b_ax][0] = sep_ax[b_ax][1] = sep_ax[b_ax][2] = 0.0
        sep_ax[b_ax][b_ax] = 1.0 # delta_ijk, for grid axes
        for g_ax in range(3):
            b_vec[b_ax][g_ax] = 0.5*box_vectors[b_ax,g_ax]
            sep_ax[b_ax + 3][g_ax] = b_vec[b_ax][g_ax] # box axes
        normalize_vector(sep_ax[b_ax + 3])
        for g_ax in range(3):
            get_cross_product(b_vec[b_ax], a_vec[g_ax], sep_ax[b_ax*3 + g_ax + 6])
            normalize_vector(sep_ax[b_ax*3 + g_ax + 6])

    for gi in range(n):
        for g_ax in range(3):
            # Calculate the separation vector
            sep_vec[g_ax] = grid_centers[gi, g_ax] - box_center[g_ax]
            # Calculate the grid axis lengths
            g_vec[g_ax][0] = g_vec[g_ax][1] = g_vec[g_ax][2] = 0.0
            g_vec[g_ax][g_ax] = 0.5 * (grid_right_edges[gi, g_ax]
                                     - grid_left_edges[gi, g_ax])
        for b_ax in range(15):
            #print b_ax,
            if check_projected_overlap(
                        sep_ax[b_ax], sep_vec, gi,
                        b_vec,  g_vec):
                good[gi] = 1
                break
    return good

def calculate_fill_grids(int fill_level, int refratio, int last_level,
                         np.ndarray[np.int64_t, ndim=1] domain_width,
                         np.ndarray[np.int64_t, ndim=1] cg_start_index,
                         np.ndarray[np.int32_t, ndim=1] cg_dims,
                         np.ndarray[np.int64_t, ndim=1] g_start_index,
                         np.ndarray[np.int32_t, ndim=1] g_dims,
                         np.ndarray[np.uint8_t, ndim=3, cast=True] g_child_mask):
    cdef np.int64_t cgstart[3]
    cdef np.int64_t gstart[3]
    cdef np.int64_t cgend[3]
    cdef np.int64_t gend[3]
    cdef np.int64_t dw[3]
    cdef np.int64_t cxi, cyi, czi, gxi, gyi, gzi, ci, cj, ck
    cdef int i, total = 0
    for i in range(3):
        dw[i] = domain_width[i]
        cgstart[i] = cg_start_index[i]
        gstart[i] = g_start_index[i]
        cgend[i] = cgstart[i] + cg_dims[i]
        gend[i] = gstart[i] + g_dims[i]
    for cxi in range(cgstart[0], cgend[0]+1):
        ci = (cxi % dw[0])
        if ci < 0: ci += dw[0]
        if ci < gstart[0]*refratio or ci >= gend[0]*refratio: continue
        gxi = (<np.int64_t> (ci / refratio)) - gstart[0]
        for cyi in range(cgstart[1], cgend[1]):
            cj = (cyi % dw[1])
            if cj < 0: cj += dw[1]
            if cj < gstart[1]*refratio or cj >= gend[1]*refratio: continue
            gyi = (<np.int64_t> (cj / refratio)) - gstart[1]
            for czi in range(cgstart[2], cgend[2]):
                ck = (czi % dw[2])
                if ck < 0: ck += dw[2]
                if ck < gstart[2]*refratio or cj >= gend[2]*refratio: continue
                gzi = (<np.int64_t> (ck / refratio)) - gstart[2]
                if last_level or g_child_mask[gxi, gyi, gzi] > 0: total += 1
    return total
