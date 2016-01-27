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

cdef class BoolArrayCollection:

    def __cinit__(self):
        cdef ewah_bool_array *ewah_keys = new ewah_bool_array()
        cdef ewah_map *ewah_coll = new ewah_map()
        self.ewah_keys = <void *> ewah_keys
        self.ewah_coll = <void *> ewah_coll

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_keys[0].set(i1)
        # Note the 0 here, for dereferencing
        ewah_coll[0][i1].set(i2)

    def set(self, i1, i2):
        self._set(i1, i2)

    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        if not ewah_keys[0].get(i1): return 0
        return ewah_coll[0][i1].get(i2)

    def get(self, i1, i2):
        return self._get(i1, i2)

    def __dealloc__(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_keys
        del ewah_coll
