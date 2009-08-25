"""
Checks for points contained in a volume

Author: John Wise <jwise77@gmail.com>
Affiliation: Princeton
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 John.  All Rights Reserved.

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


import numpy as np
cimport numpy as np
cimport cython

@cython.wraparound(False)
@cython.boundscheck(False)
def PointsInVolume(np.ndarray[np.float64_t, ndim=2] points,
                   np.ndarray[np.int8_t, ndim=1] pmask,  # pixel mask
                   np.ndarray[np.float64_t, ndim=1] left_edge,
                   np.ndarray[np.float64_t, ndim=1] right_edge,
                   np.ndarray[np.int32_t, ndim=3] mask,
                   float dx):
    cdef np.ndarray[np.int8_t, ndim=1] \
         valid = np.zeros(points.shape[0], dtype='int8')
    cdef int i, dim, count
    cdef double dx_inv
    cdef unsigned int idx[3]
    count = 0
    dx_inv = 1.0 / dx
    for i in xrange(points.shape[0]):
        if pmask[i] == 0:
            continue
        for dim in xrange(3):
            if points[i,dim] < left_edge[dim] or points[i,dim] > right_edge[dim]:
                valid[i] = 0
                break
        if valid[i] == 1:
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

                
                 
