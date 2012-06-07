"""
Shareable definitions for common fp/int Cython utilities

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
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

cimport numpy as np
cimport cython

cdef inline int imax(int i0, int i1) nogil:
    if i0 > i1: return i0
    return i1

cdef inline np.float64_t fmax(np.float64_t f0, np.float64_t f1) nogil:
    if f0 > f1: return f0
    return f1

cdef inline int imin(int i0, int i1) nogil:
    if i0 < i1: return i0
    return i1

cdef inline np.float64_t fmin(np.float64_t f0, np.float64_t f1) nogil:
    if f0 < f1: return f0
    return f1

cdef inline np.float64_t fabs(np.float64_t f0) nogil:
    if f0 < 0.0: return -f0
    return f0

cdef inline int iclip(int i, int a, int b) nogil:
    if i < a: return a
    if i > b: return b
    return i

cdef inline int i64clip(np.int64_t i, np.int64_t a, np.int64_t b) nogil:
    if i < a: return a
    if i > b: return b
    return i

cdef inline np.float64_t fclip(np.float64_t f,
                      np.float64_t a, np.float64_t b) nogil:
    return fmin(fmax(f, a), b)

