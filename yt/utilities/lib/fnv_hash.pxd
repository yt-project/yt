"""
Definitions for fnv_hash




"""



import numpy as np

cimport numpy as np


cdef np.int64_t c_fnv_hash(unsigned unsigned char[:] octets) nogil
