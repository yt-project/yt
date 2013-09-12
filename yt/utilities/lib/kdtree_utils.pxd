"""
A light interface to kdtree, from http://code.google.com/p/kdtree/



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
