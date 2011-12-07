"""
General purpose ray casting

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


cimport numpy as np
from fp_utils cimport imax, fmax, imin, fmin, iclip, fclip

ctypedef void (*ray_sampler) (np.float64_t v_pos[3],
                              np.float64_t v_dir[3],
                              np.float64_t enter_t,
                              np.float64_t exit_t,
                              int ci[3],
                              void *rdata)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int integrate_ray(np.float64_t left_edge[3],
                       np.float64_t right_edge[3],
                       np.float64_t dds[3],
                       np.float64_t idds[3],
                       int dims[3],
                       np.float64_t v_pos[3],
                       np.float64_t v_dir[3],
                       np.float64_t *return_t,
                       np.float64_t enter_t,
                       void *rdata):
    cdef int cur_ind[3], step[3], x, y, i, n, flat_ind, hit, direction
    cdef np.float64_t intersect_t = 1.0
    cdef np.float64_t iv_dir[3]
    cdef np.float64_t intersect[3], tmax[3], tdelta[3]
    cdef np.float64_t dist, alpha, dt, exit_t
    cdef np.float64_t tr, tl, temp_x, temp_y, dv
    for i in range(3):
        if (v_dir[i] < 0):
            step[i] = -1
        elif (v_dir[i] == 0):
            step[i] = 1
            tmax[i] = 1e60
            iv_dir[i] = 1e60
            tdelta[i] = 1e-60
            continue
        else:
            step[i] = 1
        x = (i+1) % 3
        y = (i+2) % 3
        iv_dir[i] = 1.0/v_dir[i]
        tl = (left_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tl*v_dir[x])
        temp_y = (v_pos[y] + tl*v_dir[y])
        if left_edge[x] <= temp_x and temp_x <= right_edge[x] and \
           left_edge[y] <= temp_y and temp_y <= right_edge[y] and \
           0.0 <= tl and tl < intersect_t:
            direction = i
            intersect_t = tl
        tr = (right_edge[i] - v_pos[i])*iv_dir[i]
        temp_x = (v_pos[x] + tr*v_dir[x])
        temp_y = (v_pos[y] + tr*v_dir[y])
        if left_edge[x] <= temp_x and temp_x <= right_edge[x] and \
           left_edge[y] <= temp_y and temp_y <= right_edge[y] and \
           0.0 <= tr and tr < intersect_t:
            direction = i
            intersect_t = tr
    if left_edge[0] <= v_pos[0] and v_pos[0] <= right_edge[0] and \
       left_edge[1] <= v_pos[1] and v_pos[1] <= right_edge[1] and \
       left_edge[2] <= v_pos[2] and v_pos[2] <= right_edge[2]:
        intersect_t = 0.0
    if enter_t >= 0.0: intersect_t = enter_t
    if not ((0.0 <= intersect_t) and (intersect_t < 1.0)): return 0
    for i in range(3):
        intersect[i] = v_pos[i] + intersect_t * v_dir[i]
        cur_ind[i] = <int> floor((intersect[i] +
                                  step[i]*1e-8*dds[i] -
                                  left_edge[i])*idds[i])
        tmax[i] = (((cur_ind[i]+step[i])*dds[i])+
                    left_edge[i]-v_pos[i])*iv_dir[i]
        # This deals with the asymmetry in having our indices refer to the
        # left edge of a cell, but the right edge of the brick being one
        # extra zone out.
        if cur_ind[i] == dims[i] and step[i] < 0:
            cur_ind[i] = dims[i] - 1
        if cur_ind[i] < 0 or cur_ind[i] >= dims[i]: return 0
        if step[i] > 0:
            tmax[i] = (((cur_ind[i]+1)*dds[i])
                        +left_edge[i]-v_pos[i])*iv_dir[i]
        if step[i] < 0:
            tmax[i] = (((cur_ind[i]+0)*dds[i])
                        +left_edge[i]-v_pos[i])*iv_dir[i]
        tdelta[i] = (dds[i]*iv_dir[i])
        if tdelta[i] < 0: tdelta[i] *= -1
    # We have to jumpstart our calculation
    enter_t = intersect_t
    hit = 0
    while 1:
        # dims here is one less than the dimensions of the data,
        # but we are tracing on the grid, not on the data...
        if (not (0 <= cur_ind[0] < dims[0])) or \
           (not (0 <= cur_ind[1] < dims[1])) or \
           (not (0 <= cur_ind[2] < dims[2])):
            break
        hit += 1
        if tmax[0] < tmax[1]:
            if tmax[0] < tmax[2]:
                exit_t = fmin(tmax[0], 1.0)
                sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind, rdata)
                cur_ind[0] += step[0]
                enter_t = tmax[0]
                tmax[0] += tdelta[0]
            else:
                exit_t = fmin(tmax[2], 1.0)
                sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind, rdata)
                cur_ind[2] += step[2]
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
        else:
            if tmax[1] < tmax[2]:
                exit_t = fmin(tmax[1], 1.0)
                sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind, rdata)
                cur_ind[1] += step[1]
                enter_t = tmax[1]
                tmax[1] += tdelta[1]
            else:
                exit_t = fmin(tmax[2], 1.0)
                sample_values(v_pos, v_dir, enter_t, exit_t, cur_ind, rdata)
                cur_ind[2] += step[2]
                enter_t = tmax[2]
                tmax[2] += tdelta[2]
        if enter_t >= 1.0: break
    if return_t != NULL: return_t[0] = exit_t
    return hit

