"""
Cython tools for working with no local fields.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np

cimport cython

from libc.math cimport sqrt

cdef void gradient(np.float64_t field, np.float64_t * kernel_values):
    kernel_values[0] = kernel_values[0] * field
    kernel_values[1] = kernel_values[1] * field
    kernel_values[2] = kernel_values[2] * field

cdef void divergence(np.float64_t * field, np.float64_t * kernel_values):
    kernel_values[0] = kernel_values[0] * field[0]
    kernel_values[1] = kernel_values[1] * field[1]
    kernel_values[2] = kernel_values[2] * field[2]

cdef void curl(np.float64_t * field, np.float64_t * kernel_values):
    kernel_values[0] = field[1] * kernel_values[2] - field[2] * kernel_values[1]
    kernel_values[0] = field[2] * kernel_values[0] - field[0] * kernel_values[2]
    kernel_values[0] = field[0] * kernel_values[1] - field[1] * kernel_values[0]
