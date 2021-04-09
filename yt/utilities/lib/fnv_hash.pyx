# distutils: libraries = STD_LIBS
# distutils: include_dirs = LIB_DIR
"""
Fast hashing routines


"""


import numpy as np

cimport cython
cimport numpy as np


@cython.wraparound(False)
@cython.boundscheck(False)
cdef np.int64_t c_fnv_hash(unsigned char[:] octets) nogil:
    # https://bitbucket.org/yt_analysis/yt/issues/1052/field-access-tests-fail-under-python3
    # FNV hash cf. http://www.isthe.com/chongo/tech/comp/fnv/index.html
    cdef np.int64_t hash_val = 2166136261
    cdef unsigned char octet
    cdef int i
    for i in range(octets.shape[0]):
        hash_val = hash_val ^ octets[i]
        hash_val = hash_val * 16777619
    return hash_val

def fnv_hash(octets):
    """

    Create a FNV hash from a bytestring.
    Info: http://www.isthe.com/chongo/tech/comp/fnv/index.html

    Parameters
    ----------
    octets : bytestring
        The string of bytes to generate a hash from.
    """
    return c_fnv_hash(octets)
