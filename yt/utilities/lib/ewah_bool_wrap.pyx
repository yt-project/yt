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

from libcpp.map cimport map as cmap
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.set cimport set as cset
from libcpp.algorithm cimport sort
from yt.utilities.lib.ewah_bool_array cimport \
    ewah_map, ewah_bool_array, sstream
from cython.operator cimport dereference, preincrement
import numpy as np

cdef extern from "<algorithm>" namespace "std" nogil:
    Iter unique[Iter](Iter first, Iter last)

cdef np.uint64_t FLAG = ~(<np.uint64_t>0)
cdef np.uint64_t MAX_VECTOR_SIZE = <np.uint64_t>1e7

# cdef class BoolArrayCollectionMap:

#     def __cinit__(self):
#         cdef map[np.uint32,BoolArrayCollection] *bool_map = 

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

    cdef int _richcmp(self, BoolArrayCollection solf, int op) except -1:

        cdef ewah_bool_array *arr1
        cdef ewah_bool_array *arr2
        cdef cmap[np.uint64_t, ewah_bool_array] *map1
        cdef cmap[np.uint64_t, ewah_bool_array] *map2
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        # == 
        if op == 2: 
            # Keys
            arr1 = <ewah_bool_array *> self.ewah_keys
            arr2 = <ewah_bool_array *> solf.ewah_keys
            if arr1[0] != arr2[0]:
                return 0
            # Refn
            arr1 = <ewah_bool_array *> self.ewah_refn
            arr2 = <ewah_bool_array *> solf.ewah_refn
            if arr1[0] != arr2[0]:
                return 0
            # Map
            map1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
            map2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
            it_map1 = map1[0].begin()
            while (it_map1 != map1[0].end()):
                it_map2 = map2[0].find(dereference(it_map1).first)
                if it_map2 == map2[0].end():
                    return 0
                if dereference(it_map1).second != dereference(it_map2).second:
                    return 0
                preincrement(it_map1)
            it_map2 =map2[0].begin()
            while (it_map2 != map2[0].end()):
                it_map1 = map1[0].find(dereference(it_map2).first)
                if it_map1 == map1[0].end():
                    return 0
                if dereference(it_map2).second != dereference(it_map1).second:
                    return 0
                preincrement(it_map2)
            # Match
            return 1
        # !=
        elif op == 3:
            if self._richcmp(solf, 2) == 1:
                return 0
            return 1
        else:
            return -1
            # options = ['<','<=','==','!=','>','>=']
            # raise NotImplementedError("Operator {} is not yet implemented.".format(options[op]))

    def __richcmp__(BoolArrayCollection self, BoolArrayCollection solf, int op):
        if self._richcmp(solf, op) == 1:
            return True
        else:
            return False

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

    cdef void _set_coarse(self, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        ewah_keys[0].set(i1)

    def set_coarse(self, i1):
        return self._set_coarse(i1)

    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        ewah_refn[0].set(i1)
        ewah_coll[0][i1].set(i2)

    cdef void _set_coarse_array(self, np.uint8_t[:] arr):
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                self._set_coarse(i1)

    cdef void _set_refined_array(self, np.uint64_t i1, np.uint8_t[:] arr):
        cdef np.uint64_t i2
        for i2 in range(arr.shape[0]):
            if arr[i2] == 1:
                self._set_refined(i1, i2)

    def set_refined(self, i1, i2):
        return self._set_refined(i1, i2)
        
    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_coll[0][i1].set(i2)

    def set_map(self, i1, i2):
        self._set_map(i1, i2)

    cdef void _set_refn(self, np.uint64_t i1):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        ewah_refn[0].set(i1)

    def set_refn(self, i1):
        self._set_refn(i1)

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

    cdef bint _get_coarse(self, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        if not ewah_keys[0].get(i1): return 0
        return 1

    def get_coarse(self, i1):
        return self._get_coarse(i1)

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

    cdef int _count_total(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef int out
        out = ewah_keys.numberOfOnes()
        return out

    def count_total(self):
        return self._count_total()

    cdef int _count_refined(self):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef int out
        out = ewah_refn.numberOfOnes()
        return out

    def count_refined(self):
        return self._count_refined()

    cdef int _count_coarse(self):
        self._ewah_coarse()
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        cdef int out
        out = ewah_coar.numberOfOnes()
        return out

    def count_coarse(self):
        return self._count_coarse()

    cdef void _append(self, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array swap, mi1_ewah1, mi1_ewah2
        cdef np.uint64_t nrefn, mi1
        # Keys
        ewah_keys1[0].logicalor(ewah_keys2[0], swap)
        ewah_keys1[0].swap(swap)
        # Refined
        ewah_refn1[0].logicalor(ewah_refn2[0], swap)
        ewah_refn1[0].swap(swap)
        # Map
        it_map2 = ewah_coll2[0].begin()
        while it_map2 != ewah_coll2[0].end():
            mi1 = dereference(it_map2).first
            mi1_ewah2 = dereference(it_map2).second
            it_map1 = ewah_coll1[0].find(mi1)
            if it_map1 == ewah_coll1[0].end():
                ewah_coll1[0][mi1] = mi1_ewah2
            else:
                mi1_ewah1 = dereference(it_map1).second
                mi1_ewah1.logicalor(mi1_ewah2, swap)
                mi1_ewah1.swap(swap)
            preincrement(it_map2)

    def append(self, solf):
        return self._append(solf)

    cdef bint _intersects(self, BoolArrayCollection solf):
        # self._ewah_coarse()
        # solf._ewah_coarse()
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        # cdef ewah_bool_array *ewah_coar1 = <ewah_bool_array *> self.ewah_coar
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        # cdef ewah_bool_array *ewah_coar2 = <ewah_bool_array *> solf.ewah_coar
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        # No intersection
        if ewah_keys1[0].intersects(ewah_keys2[0]) == 0:
            return 0
        # Intersection at coarse level
        # if ewah_coar1[0].intersects(ewah_coar2[0]) == 1:
        #     return 1
        # Intersection at refined level
        if ewah_refn1[0].intersects(ewah_refn2[0]) == 1:
            it_map1 = ewah_coll1[0].begin()
            while (it_map1 != ewah_coll1[0].end()):
                mi1 = dereference(it_map1).first
                it_map2 = ewah_coll2[0].find(mi1)
                if it_map2 != ewah_coll2[0].end():
                    mi1_ewah1 = dereference(it_map1).second
                    mi1_ewah2 = dereference(it_map2).second
                    if mi1_ewah1.intersects(mi1_ewah2):
                        return 1
                preincrement(it_map1)
        # Intersection at coarse level or refined inside coarse
        else:
            return 1
        return 0

    cdef bytes _dumps(self):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map
        cdef np.uint64_t nrefn, mi1
        cdef ewah_bool_array mi1_ewah
        # Write mi1 ewah & refinment ewah
        ewah_keys[0].write(ss,1)
        ewah_refn[0].write(ss,1)
        # Number of refined bool arrays
        nrefn = <np.uint64_t>(ewah_refn[0].numberOfOnes())
        ss.write(<const char *> &nrefn, sizeof(nrefn))
        # Loop over refined bool arrays
        it_map = ewah_coll[0].begin()
        while it_map != ewah_coll[0].end():
            mi1 = dereference(it_map).first
            mi1_ewah = dereference(it_map).second
            ss.write(<const char *> &mi1, sizeof(mi1))
            mi1_ewah.write(ss,1)
            preincrement(it_map)
        # Return type cast python bytes string
        return <bytes>ss.str()

    def dumps(self):
        return self._dumps()

    cdef void _loads(self, bytes s):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map
        cdef np.uint64_t nrefn, mi1
        cdef ewah_bool_array mi1_ewah
        cdef int i
        # Write string to string stream
        ss.write(s, len(s))
        # Read keys and refinment arrays
        ewah_keys[0].read(ss,1)
        ewah_refn[0].read(ss,1)
        # Read and check number of refined cells
        ss.read(<char *> (&nrefn), sizeof(nrefn))
        if nrefn != ewah_refn[0].numberOfOnes():
            raise Exception("Error in read. File indicates {} refinements, but bool array has {}.".format(nrefn,ewah_refn[0].numberOfOnes()))
        # Loop over refined cells
        for i in range(nrefn):
            ss.read(<char *> (&mi1), sizeof(mi1))
            ewah_coll[0][mi1].read(ss,1)
            # or...
            #mi1_ewah.read(ss,1)
            #ewah_coll[0][mi1].swap(mi1_ewah)

    def loads(self, s):
        return self._loads(s)

    def __dealloc__(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_keys
        del ewah_refn
        del ewah_coar
        del ewah_coll

    def print_info(self, prefix=''):
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(prefix,
                                                                     self._count_coarse(),
                                                                     self._count_refined(),
                                                                     self._count_total()))

# Vector version
cdef class SparseUnorderedBitmaskVector:
    def __cinit__(self):
        cdef vector[np.uint64_t] *entries = new vector[np.uint64_t]()
        self.entries = <void *> entries
        self.total = 0

    cdef void _set(self, np.uint64_t ind):
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        entries[0].push_back(ind)
        self.total += 1

    def set(self, ind):
        self._set(ind)

    cdef void _fill(self, np.uint8_t[:] mask):
        cdef np.uint64_t i, ind
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        for i in range(entries[0].size()):
            ind = entries[0][i]
            mask[ind] = 1

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        self._remove_duplicates()
        cdef np.uint64_t i, ind
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        for i in range(entries[0].size()):
            ind = entries[0][i]
            mm._set_coarse(ind)

    cdef void _reset(self):
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        entries[0].erase(entries[0].begin(), entries[0].end())
        self.total = 0

    cdef to_array(self):
        self._remove_duplicates()
        cdef np.ndarray[np.uint64_t, ndim=1] rv
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        rv = np.empty(entries[0].size(), dtype='uint64')
        for i in range(entries[0].size()):
            rv[i] = entries[0][i]
        return rv

    cdef void _remove_duplicates(self):
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        cdef vector[np.uint64_t].iterator last
        sort(entries[0].begin(), entries[0].end())
        last = unique(entries[0].begin(), entries[0].end())
        entries[0].erase(last, entries[0].end())

    cdef void _prune(self):
        if self.total > MAX_VECTOR_SIZE:
            self._remove_duplicates()
            self.total = 0

    def __dealloc__(self):
        cdef vector[np.uint64_t] *entries = <vector[np.uint64_t]*> self.entries
        del entries

# Set version
cdef class SparseUnorderedBitmaskSet:
    def __cinit__(self):
        cdef cset[np.uint64_t] *entries = new cset[np.uint64_t]()
        self.entries = <void *> entries

    cdef void _set(self, np.uint64_t ind):
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        entries[0].insert(ind)

    def set(self, ind):
        self._set(ind)

    cdef void _fill(self, np.uint8_t[:] mask):
        cdef np.uint64_t ind
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        cdef cset[np.uint64_t].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            ind = dereference(it)
            mask[ind] = 1
            preincrement(it)

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        cdef np.uint64_t ind
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        cdef cset[np.uint64_t].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            ind = dereference(it)
            mm._set_coarse(ind)
            preincrement(it)

    cdef void _reset(self):
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        entries[0].clear()

    cdef to_array(self):
        cdef np.uint64_t ind
        cdef np.ndarray[np.uint64_t, ndim=1] rv
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        cdef cset[np.uint64_t].iterator it
        rv = np.empty(entries[0].size(), dtype='uint64')
        it = entries[0].begin()
        i = 0
        while it != entries[0].end():
            ind = dereference(it)
            rv[i] = ind
            preincrement(it)
            i += 1
        return rv

    def __dealloc__(self):
        cdef cset[np.uint64_t] *entries = <cset[np.uint64_t]*> self.entries
        del entries

# vector version
cdef class SparseUnorderedRefinedBitmaskVector:
    def __cinit__(self):
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = new vector[pair[np.uint64_t,np.uint64_t]]()
        self.entries = <void *> entries
        self.total = 0

    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2):
        cdef pair[np.uint64_t,np.uint64_t] ind
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        ind.first = ind1
        ind.second = ind2
        entries[0].push_back(ind)
        self.total += 1


    def set(self, ind1, ind2):
        self._set(ind1, ind2)

    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:] mask2):
        cdef np.uint64_t i, ind
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef vector[pair[np.uint64_t,np.uint64_t]].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            ind = dereference(it).first
            mask1[ind] = 1
            ind = dereference(it).second
            mask2[ind] = 1
            preincrement(it)

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        self._remove_duplicates()
        cdef np.uint64_t mi1, mi2
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef vector[pair[np.uint64_t,np.uint64_t]].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            mi1 = dereference(it).first
            mi2 = dereference(it).second
            mm._set_refined(mi1, mi2)
            preincrement(it)

    cdef void _reset(self):
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        entries[0].erase(entries[0].begin(), entries[0].end())
        self.total = 0

    cdef to_array(self):
        cdef int i
        cdef np.ndarray[np.uint64_t, ndim=2] rv
        self._remove_duplicates()
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef vector[pair[np.uint64_t,np.uint64_t]].iterator it
        rv = np.empty((entries[0].size(),2),dtype='uint64')
        it = entries[0].begin()
        i = 0
        while it != entries[0].end():
            rv[i,0] = dereference(it).first
            rv[i,1] = dereference(it).second
            i += 1
            preincrement(it)
        return rv

    cdef void _remove_duplicates(self):
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef vector[pair[np.uint64_t,np.uint64_t]].iterator last
        sort(entries[0].begin(), entries[0].end())
        last = unique(entries[0].begin(), entries[0].end())
        entries[0].erase(last, entries[0].end())
        # http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        # cdef np.ndarray[np.uint64_t, ndim=2] rv
        # cdef np.ndarray[np.uint64_t, ndim=2] rv_uni
        # cdef np.uint64_t m
        # cdef vector[np.uint64_t].iterator last1
        # cdef vector[np.uint64_t].iterator last2
        # # cdef np.ndarray[np.uint64_t, ndim=1] _
        # cdef vector[np.uint64_t] *entries1 = <vector[np.uint64_t]*> self.entries1
        # cdef vector[np.uint64_t] *entries2 = <vector[np.uint64_t]*> self.entries2
        # rv = np.empty((entries1[0].size(),2),dtype='uint64')
        # for i in range(entries1[0].size()):
        #     rv[i,0] = entries1[0][i]
        #     rv[i,1] = entries2[0][i]
        # rv_uni = np.unique(np.ascontiguousarray(rv).view(np.dtype((np.void, rv.dtype.itemsize * rv.shape[1])))).view(rv.dtype).reshape(-1,rv.shape[1])
        # last1 = entries1[0].begin() + rv_uni.shape[0]
        # last2 = entries2[0].begin() + rv_uni.shape[0]
        # for m in range(rv_uni.shape[0]):
        #     entries1[0][m] = rv_uni[m,0]
        #     entries2[0][m] = rv_uni[m,1]
        # entries1[0].erase(last1, entries1[0].end())
        # entries2[0].erase(last2, entries2[0].end())

    cdef void _prune(self):
        if self.total > MAX_VECTOR_SIZE:
            self._remove_duplicates()
            self.total = 0

    def __dealloc__(self):
        cdef vector[pair[np.uint64_t,np.uint64_t]] *entries = <vector[pair[np.uint64_t,np.uint64_t]]*> self.entries
        del entries

# Set version
cdef class SparseUnorderedRefinedBitmaskSet:
    def __cinit__(self):
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = new cset[pair[np.uint64_t,np.uint64_t]]()
        self.entries = <void *> entries

    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2):
        cdef pair[np.uint64_t,np.uint64_t] ind
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        ind.first = ind1
        ind.second = ind2
        entries[0].insert(ind)

    def set(self, ind1, ind2):
        self._set(ind1, ind2)

    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:] mask2):
        cdef np.uint64_t ind
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef cset[pair[np.uint64_t,np.uint64_t]].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            ind = dereference(it).first
            mask1[ind] = 1
            ind = dereference(it).second
            mask2[ind] = 1
            preincrement(it)

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        cdef np.uint64_t mi1, mi2
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef cset[pair[np.uint64_t,np.uint64_t]].iterator it
        it = entries[0].begin()
        while it != entries[0].end():
            mi1 = dereference(it).first
            mi2 = dereference(it).second
            mm._set_refined(mi1, mi2)
            preincrement(it)

    cdef void _reset(self):
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        entries[0].clear()

    cdef to_array(self):
        cdef int i
        cdef np.ndarray[np.uint64_t, ndim=2] rv
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        cdef cset[pair[np.uint64_t,np.uint64_t]].iterator it
        rv = np.empty((entries[0].size(),2),dtype='uint64')
        it = entries[0].begin()
        i = 0
        while it != entries[0].end():
            rv[i,0] = dereference(it).first
            rv[i,1] = dereference(it).second
            i += 1
            preincrement(it)
        return rv

    def __dealloc__(self):
        cdef cset[pair[np.uint64_t,np.uint64_t]] *entries = <cset[pair[np.uint64_t,np.uint64_t]]*> self.entries
        del entries

