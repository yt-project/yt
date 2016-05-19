"""
Some simple operations for operating on ragged arrays



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython

cdef fused numpy_dt:
    np.float32_t
    np.float64_t
    np.int32_t
    np.int64_t

cdef numpy_dt r_min(numpy_dt a, numpy_dt b):
    if a < b: return a
    return b

cdef numpy_dt r_max(numpy_dt a, numpy_dt b):
    if a > b: return a
    return b

cdef numpy_dt r_add(numpy_dt a, numpy_dt b):
    return a + b

cdef numpy_dt r_subtract(numpy_dt a, numpy_dt b):
    return a - b

cdef numpy_dt r_multiply(numpy_dt a, numpy_dt b):
    return a * b

@cython.cdivision(True)
cdef numpy_dt r_divide(numpy_dt a, numpy_dt b):
    return a / b

def index_unop(np.ndarray[numpy_dt, ndim=1] values,
              np.ndarray[np.int64_t, ndim=1] indices,
              np.ndarray[np.int64_t, ndim=1] sizes,
              operation):
    cdef numpy_dt mi, ma
    if numpy_dt == np.float32_t:
        dt = "float32"
        mi = np.finfo(dt).min
        ma = np.finfo(dt).max
    elif numpy_dt == np.float64_t:
        dt = "float64"
        mi = np.finfo(dt).min
        ma = np.finfo(dt).max
    elif numpy_dt == np.int32_t:
        dt = "int32"
        mi = np.iinfo(dt).min
        ma = np.iinfo(dt).max
    elif numpy_dt == np.int64_t:
        dt = "int64"
        mi = np.iinfo(dt).min
        ma = np.iinfo(dt).max
    cdef np.ndarray[numpy_dt] out_values = np.zeros(sizes.size, dtype=dt)
    cdef numpy_dt (*func)(numpy_dt a, numpy_dt b)
    # Now we figure out our function.  At present, we only allow addition and
    # multiplication, because they are commutative and easy to bootstrap.
    cdef numpy_dt ival, val
    if operation == "sum":
        ival = 0
        func = r_add
    elif operation == "prod":
        ival = 1
        func = r_multiply
    elif operation == "max":
        ival = mi
        func = r_max
    elif operation == "min":
        ival = ma
        func = r_min
    else:
        raise NotImplementedError
    cdef np.int64_t i, ind_ind, ind_arr
    ind_ind = 0
    for i in range(sizes.size):
        # Each entry in sizes is the size of the array
        val = ival
        for _ in range(sizes[i]):
            ind_arr = indices[ind_ind]
            val = func(val, values[ind_arr])
            ind_ind += 1
        out_values[i] = val
    return out_values
