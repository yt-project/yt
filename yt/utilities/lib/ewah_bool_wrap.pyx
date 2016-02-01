"""
Wrapper for EWAH Bool Array: https://github.com/lemire/EWAHBoolArray



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from libcpp.map cimport map
from yt.utilities.lib.ewah_bool_array cimport \
    ewah_map, ewah_bool_array

cdef np.uint64_t FLAG = ~(<np.uint64_t>0)

cdef class BoolArrayCollection:

    def __cinit__(self):
        cdef ewah_bool_array *ewah_keys = new ewah_bool_array()
        cdef ewah_bool_array *ewah_refn = new ewah_bool_array()
        cdef ewah_bool_array *ewah_coar = new ewah_bool_array()
        cdef ewah_map *ewah_coll = new ewah_map()
        self.ewah_keys = <void *> ewah_keys
        self.ewah_refn = <void *> ewah_refn
        self.ewah_coar = <void *> ewah_coar
        self.ewah_coll = <void *> ewah_coll

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_keys[0].set(i1)
        # Note the 0 here, for dereferencing
        if i2 != FLAG:
            ewah_refn[0].set(i1)
            ewah_coll[0][i1].set(i2)

    def set(self, i1, i2 = FLAG):
        self._set(i1, i2)

    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        if not ewah_keys[0].get(i1): return 0
        if not ewah_refn[0].get(i1) or (i2 == FLAG): 
            return 1
        return ewah_coll[0][i1].get(i2)

    def get(self, i1, i2 = FLAG):
        return self._get(i1, i2)

    cdef bint _contains(self, np.uint64_t i):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        return ewah_keys[0].get(i)

    def contains(self, np.uint64_t i):
        return self._contains(i)

    cdef bint _isref(self, np.uint64_t i):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        return ewah_refn[0].get(i)

    def isref(self, np.uint64_t i):
        return self._isref(i)

    cdef void _ewah_coarse(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        ewah_coar[0].reset()
        ewah_keys[0].logicalxor(ewah_refn[0],ewah_coar[0])
        return

    def ewah_coarse(self):
        return self._ewah_coarse()

    def __dealloc__(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_keys
        del ewah_refn
        del ewah_coar
        del ewah_coll

