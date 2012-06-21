"""
Fixed interpolator includes

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


cdef extern from "FixedInterpolator.h":
    np.float64_t fast_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                  np.float64_t *data) nogil
    np.float64_t offset_interpolate(int ds[3], np.float64_t dp[3],
                                    np.float64_t *data) nogil
    np.float64_t trilinear_interpolate(int ds[3], int ci[3], np.float64_t dp[3],
                                       np.float64_t *data) nogil
    void eval_gradient(int ds[3], np.float64_t dp[3], np.float64_t *data,
                       np.float64_t grad[3]) nogil
    void offset_fill(int *ds, np.float64_t *data, np.float64_t *gridval) nogil
    void vertex_interp(np.float64_t v1, np.float64_t v2, np.float64_t isovalue,
                       np.float64_t vl[3], np.float64_t dds[3],
                       np.float64_t x, np.float64_t y, np.float64_t z,
                       int vind1, int vind2) nogil

