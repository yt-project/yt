"""
A light interface to kdtree, from http://code.google.com/p/kdtree/

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

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "kdtree.h":
    struct kdtree
    struct kdres

    kdtree *kd_create(int k)
    void kd_free(kdtree *tree)
    
    int kd_insert3(kdtree *tree, np.float64_t x, np.float64_t y, np.float64_t z, void *data)
    kdres *kd_nearest3(kdtree *tree, np.float64_t x, np.float64_t y,
                       np.float64_t z) nogil

    kdres *kd_nearest_range3(kdtree *tree, np.float64_t x, np.float64_t y, np.float64_t z,
                             np.float64_t range) nogil

    void kd_res_free(kdres *set) nogil
    int kd_res_size(kdres *set) nogil
    int kd_res_next(kdres *set) nogil
    void kd_res_rewind(kdres *set) nogil

    void kd_res_item3(kdres *set, np.float64_t *x, np.float64_t *y,
                      np.float64_t *z) nogil
    void *kd_res_item_data(kdres *set) nogil

    void kd_data_destructor(kdtree *tree, void (*destr)(void*))
