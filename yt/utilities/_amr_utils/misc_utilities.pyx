"""
Simple utilities that don't fit anywhere else

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

def get_color_bounds(np.ndarray[np.float64_t, ndim=1] px,
                     np.ndarray[np.float64_t, ndim=1] py,
                     np.ndarray[np.float64_t, ndim=1] pdx,
                     np.ndarray[np.float64_t, ndim=1] pdy,
                     np.ndarray[np.float64_t, ndim=1] value,
                     np.float64_t leftx, np.float64_t rightx,
                     np.float64_t lefty, np.float64_t righty,
                     np.float64_t mindx = -1):
    cdef int i
    cdef np.float64_t mi = 1e100, ma = -1e100, v
    cdef int np = px.shape[0]
    for i in range(np):
        v = value[i]
        if v < mi or v > ma:
            if px[i] + pdx[i] < leftx: continue
            if px[i] - pdx[i] > rightx: continue
            if py[i] + pdy[i] < lefty: continue
            if py[i] - pdy[i] > righty: continue
            if pdx[i] < mindx or pdy[i] < mindx: continue
            if v < mi: mi = v
            if v > ma: ma = v
    return (mi, ma)
