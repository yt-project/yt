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

import struct
from libcpp.map cimport map as cmap
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.set cimport set as cset
from libcpp.map cimport map
from libcpp.algorithm cimport sort
from libc.stdlib cimport malloc, free, qsort
from yt.utilities.lib.ewah_bool_array cimport \
    ewah_map, ewah_bool_array, sstream
from cython.operator cimport dereference, preincrement
import numpy as np
cimport numpy as np
cimport cython

cdef extern from "<algorithm>" namespace "std" nogil:
    Iter unique[Iter](Iter first, Iter last)

cdef np.uint64_t FLAG = ~(<np.uint64_t>0)
cdef np.uint64_t MAX_VECTOR_SIZE = <np.uint64_t>1e7

DEF UncompressedFormat = 'Pointer'

#ctypedef np.uint8_t bitarrtype
ctypedef bint bitarrtype

cdef class FileBitmasks:

    def __cinit__(self, np.uint32_t nfiles):
        cdef int i
        self.nfiles = nfiles
        cdef ewah_bool_array **ewah_keys = <ewah_bool_array **>malloc(nfiles*sizeof(ewah_bool_array*))
        cdef ewah_bool_array **ewah_refn = <ewah_bool_array **>malloc(nfiles*sizeof(ewah_bool_array*))
        cdef ewah_bool_array **ewah_owns = <ewah_bool_array **>malloc(nfiles*sizeof(ewah_bool_array*))
        cdef ewah_map **ewah_coll = <ewah_map **>malloc(nfiles*sizeof(ewah_map*))
        for i in range(nfiles):
            ewah_keys[i] = new ewah_bool_array()
            ewah_refn[i] = new ewah_bool_array()
            ewah_owns[i] = new ewah_bool_array()
            ewah_coll[i] = new ewah_map()
        self.ewah_keys = <void **>ewah_keys
        self.ewah_refn = <void **>ewah_refn
        self.ewah_owns = <void **>ewah_owns
        self.ewah_coll = <void **>ewah_coll

    def reset(self):
        self.__dealloc__()
        self.__init__(self.nfiles)

    cdef bint _iseq(self, FileBitmasks solf):
        cdef np.int32_t ifile
        cdef ewah_bool_array* arr1
        cdef ewah_bool_array* arr2
        cdef cmap[np.uint64_t, ewah_bool_array] *map1
        cdef cmap[np.uint64_t, ewah_bool_array] *map2
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        if self.nfiles != solf.nfiles:
            return 0
        for ifile in range(self.nfiles): 
            # Keys
            arr1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
            arr2 = (<ewah_bool_array **> solf.ewah_keys)[ifile]
            if arr1[0] != arr2[0]:
                return 0
            # Refn
            arr1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
            arr2 = (<ewah_bool_array **> solf.ewah_refn)[ifile]
            if arr1[0] != arr2[0]:
                return 0
            # Map
            map1 = (<cmap[np.uint64_t, ewah_bool_array] **> self.ewah_coll)[ifile]
            map2 = (<cmap[np.uint64_t, ewah_bool_array] **> solf.ewah_coll)[ifile]
            it_map1 = map1[0].begin()
            while (it_map1 != map1[0].end()):
                it_map2 = map2[0].find(dereference(it_map1).first)
                if it_map2 == map2[0].end():
                    return 0
                if dereference(it_map1).second != dereference(it_map2).second:
                    return 0
                preincrement(it_map1)
            it_map2 = map2[0].begin()
            while (it_map2 != map2[0].end()):
                it_map1 = map1[0].find(dereference(it_map2).first)
                if it_map1 == map1[0].end():
                    return 0
                if dereference(it_map2).second != dereference(it_map1).second:
                    return 0
                preincrement(it_map2)
            # Match
            return 1
        
    def iseq(self, solf):
        return self._iseq(solf)

    cdef BoolArrayCollection _get_bitmask(self, np.uint32_t ifile):
        cdef BoolArrayCollection out = BoolArrayCollection()
        cdef ewah_bool_array **ewah_keys = <ewah_bool_array **>self.ewah_keys
        cdef ewah_bool_array **ewah_refn = <ewah_bool_array **>self.ewah_refn
        cdef ewah_bool_array **ewah_owns = <ewah_bool_array **>self.ewah_owns
        cdef ewah_map **ewah_coll = <ewah_map **>self.ewah_coll
        out.ewah_keys = <void *>ewah_keys[ifile]
        out.ewah_refn = <void *>ewah_refn[ifile]
        out.ewah_owns = <void *>ewah_owns[ifile]
        out.ewah_coll = <void *>ewah_coll[ifile]
        # TODO: make sure ewah arrays are not deallocated when out is clean up
        return out

    cdef tuple _find_collisions(self, BoolArrayCollection coll, bint verbose = 0):
        cdef tuple cc, cr
        cc = self._find_collisions_coarse(coll, verbose)
        cr = self._find_collisions_refined(coll, verbose)
        return cc, cr

    cdef tuple _find_collisions_coarse(self, BoolArrayCollection coll, bint
                        verbose = 0, file_list = None):
        cdef np.int32_t ifile
        cdef ewah_bool_array arr_two, arr_swap, arr_keys, arr_refn
        cdef ewah_bool_array* iarr
        cdef ewah_bool_array* coll_keys
        cdef ewah_bool_array* coll_refn
        coll_keys = (<ewah_bool_array*> coll.ewah_keys)
        coll_refn = (<ewah_bool_array*> coll.ewah_refn)
        if file_list is None:
            file_list = range(self.nfiles)
        for ifile in file_list:
            iarr = (<ewah_bool_array **>self.ewah_keys)[ifile]
            arr_keys.logicaland(iarr[0], arr_two)
            arr_keys.logicalor(iarr[0], arr_swap)
            arr_keys.swap(arr_swap)
            arr_refn.logicalor(arr_two, arr_swap)
            arr_refn.swap(arr_swap)
        coll_keys[0].swap(arr_keys)
        coll_refn[0].swap(arr_refn)
        # Print
        cdef int nc, nm
        nc = coll_refn[0].numberOfOnes()
        nm = coll_keys[0].numberOfOnes()
        cdef tuple nout = (nc, nm)
        if verbose == 1:
            print("{: 10d}/{: 10d} collisions at coarse refinement.  ({: 10.5f}%)".format(nc,nm,100.0*float(nc)/nm))
        return nout

    cdef tuple _find_collisions_refined(self, BoolArrayCollection coll, bint verbose = 0):
        cdef np.int32_t ifile
        cdef ewah_bool_array iarr, arr_two, arr_swap
        cdef ewah_bool_array* coll_refn
        cdef map[np.uint64_t, ewah_bool_array] map_keys, map_refn
        cdef map[np.uint64_t, ewah_bool_array]* coll_coll
        cdef map[np.uint64_t, ewah_bool_array]* map_bitmask
        coll_refn = <ewah_bool_array*> coll.ewah_refn
        if coll_refn[0].numberOfOnes() == 0:
            if verbose == 1:
                print("{: 10d}/{: 10d} collisions at refined refinement. ({: 10.5f}%)".format(0,0,0))
            return (0,0)
        coll_coll = <map[np.uint64_t, ewah_bool_array]*> coll.ewah_coll
        for ifile in range(self.nfiles):
            map_bitmask = (<map[np.uint64_t, ewah_bool_array]**> self.ewah_coll)[ifile]
            it_mi1 = map_bitmask[0].begin()
            while it_mi1 != map_bitmask[0].end():
                mi1 = dereference(it_mi1).first
                iarr = dereference(it_mi1).second
                map_keys[mi1].logicaland(iarr, arr_two)
                map_keys[mi1].logicalor(iarr, arr_swap)
                map_keys[mi1].swap(arr_swap)
                map_refn[mi1].logicalor(arr_two, arr_swap)
                map_refn[mi1].swap(arr_swap)
                preincrement(it_mi1)
        coll_coll[0] = map_refn
        # Count
        cdef int nc, nm
        nc = 0
        nm = 0
        it_mi1 = map_refn.begin()
        while it_mi1 != map_refn.end():
            mi1 = dereference(it_mi1).first
            iarr = dereference(it_mi1).second
            nc += iarr.numberOfOnes()
            iarr = map_keys[mi1]
            nm += iarr.numberOfOnes()
            preincrement(it_mi1)
        cdef tuple nout = (nc, nm)
        # Print
        if verbose == 1:
            if nm == 0:
                print("{: 10d}/{: 10d} collisions at refined refinement. ({: 10.5f}%)".format(nc,nm,0.0))
            else:
                print("{: 10d}/{: 10d} collisions at refined refinement. ({: 10.5f}%)".format(nc,nm,100.0*float(nc)/nm))
        return nout

    cdef void _set(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        ewah_keys[0].set(i1)
        if i2 != FLAG:
            ewah_refn[0].set(i1)
            ewah_coll[0][i1].set(i2)

    cdef void _set_coarse(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        ewah_keys[0].set(i1)

    cdef void _set_refined(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        ewah_refn[0].set(i1)
        ewah_coll[0][i1].set(i2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_owners(self, np.uint32_t[:,:] arr):
        cdef ewah_bool_array **ewah_owns = <ewah_bool_array **> self.ewah_owns
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1][0] != 0:
                ewah_owns[arr[i1][1]][0].set(i1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_coarse_array(self, np.uint32_t ifile, np.uint8_t[:] arr):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                ewah_keys[0].set(i1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_refined_array(self, np.uint32_t ifile, np.uint64_t i1, np.uint8_t[:] arr):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        cdef np.uint64_t i2
        for i2 in range(arr.shape[0]):
            if arr[i2] == 1:
                ewah_refn[0].set(i1)
                ewah_coll[0][i1].set(i2)

    cdef void _set_map(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        ewah_coll[0][i1].set(i2)

    cdef void _set_refn(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        ewah_refn[0].set(i1)

    cdef void _set_owns(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_owns = (<ewah_bool_array **> self.ewah_owns)[ifile]
        ewah_owns[0].set(i1)

    cdef bint _get(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        if (ewah_keys[0].get(i1) == 0): return 0
        if (ewah_refn[0].get(i1) == 0) or (i2 == FLAG): 
            return 1
        return ewah_coll[0][i1].get(i2)

    cdef bint _get_coarse(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        return ewah_keys[0].get(i1)

    cdef bint _isref(self, np.uint32_t ifile, np.uint64_t i):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        return ewah_refn[0].get(i)

    cdef int _count_coarse(self, np.uint32_t ifile):
        return self._count_total(ifile) - self._count_refined(ifile)

    cdef int _count_total(self, np.uint32_t ifile):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef int out
        out = ewah_keys.numberOfOnes()
        return out

    cdef int _count_refined(self, np.uint32_t ifile):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef int out
        out = ewah_refn.numberOfOnes()
        return out

    cdef int _count_owned(self, np.uint32_t ifile):
        cdef ewah_bool_array *ewah_owns = (<ewah_bool_array **> self.ewah_owns)[ifile]
        cdef int out
        out = ewah_owns.numberOfOnes()
        return out

    cdef void _append(self, np.uint32_t ifile, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll1 = (<ewah_map **> self.ewah_coll)[ifile]
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

    cdef bint _intersects(self, np.uint32_t ifile, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll1 = (<ewah_map **> self.ewah_coll)[ifile]
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # No intersection
        if ewah_keys1[0].intersects(ewah_keys2[0]) == 0:
            return 0
        # Intersection at coarse level
        ewah_keys1[0].logicalxor(ewah_refn1[0],ewah_coar1)
        ewah_keys2[0].logicalxor(ewah_refn2[0],ewah_coar2)
        if ewah_coar1.intersects(ewah_keys2[0]) == 1:
            return 1
        if ewah_coar2.intersects(ewah_keys1[0]) == 1:
            return 1
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
        return 0
    
    cdef void _logicalxor(self, np.uint32_t ifile, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll1 = (<ewah_map **> self.ewah_coll)[ifile]
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # Keys
        ewah_keys1[0].logicalxor(ewah_keys2[0],ewah_keys_out[0])
        # Refn
        ewah_refn1[0].logicalxor(ewah_refn2[0],ewah_refn_out[0])
        # Coll
        it_map1 = ewah_coll1[0].begin()
        while (it_map1 != ewah_coll1[0].end()):
            mi1 = dereference(it_map1).first
            mi1_ewah1 = dereference(it_map1).second
            it_map2 = ewah_coll2[0].find(mi1)
            if it_map2 == ewah_coll2[0].end():
                ewah_coll_out[0][mi1] = mi1_ewah1
            else:
                mi1_ewah2 = dereference(it_map2).second
                mi1_ewah1.logicalxor(mi1_ewah2, swap)
                ewah_coll_out[0][mi1] = swap
            preincrement(it_map1)
        it_map2 = ewah_coll2[0].begin()
        while (it_map2 != ewah_coll2[0].end()):
            mi1 = dereference(it_map2).first
            mi1_ewah2 = dereference(it_map2).second
            it_map1 = ewah_coll1[0].find(mi1)
            if it_map1 == ewah_coll1[0].end():
                ewah_coll_out[0][mi1] = mi1_ewah2
            preincrement(it_map2)

    def logicalxor(self, ifile, solf, out):
        return self._logicalxor(ifile, solf, out)

    cdef void _logicaland(self, np.uint32_t ifile, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll1 = (<ewah_map **> self.ewah_coll)[ifile]
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # Keys
        ewah_keys1[0].logicaland(ewah_keys2[0],ewah_keys_out[0])
        # Refn
        ewah_refn1[0].logicaland(ewah_refn2[0],ewah_refn_out[0])
        # Coll
        if ewah_refn_out[0].numberOfOnes() > 0:
            it_map1 = ewah_coll1[0].begin()
            while (it_map1 != ewah_coll1[0].end()):
                mi1 = dereference(it_map1).first
                it_map2 = ewah_coll2[0].find(mi1)
                if it_map2 != ewah_coll2[0].end():
                    mi1_ewah1 = dereference(it_map1).second
                    mi1_ewah2 = dereference(it_map2).second
                    mi1_ewah1.logicaland(mi1_ewah2, swap)
                    ewah_coll_out[0][mi1] = swap
                preincrement(it_map1)

    def logicaland(self, ifile, solf, out):
        return self._logicaland(ifile, solf, out)

    cdef bytes _dumps(self, np.uint32_t ifile):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_bool_array *ewah_owns = (<ewah_bool_array **> self.ewah_owns)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
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
        # Owners EWAH
        ewah_owns[0].write(ss,1)
        # Return type cast python bytes string
        return <bytes>ss.str()

    cdef bint _loads(self, np.uint32_t ifile, bytes s):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_bool_array *ewah_owns = (<ewah_bool_array **> self.ewah_owns)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
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
        # Owners
        if ss.eof() == 1: return 0
        ewah_owns[0].read(ss,1)
        return 1

    # def save(self, fname):
    #     cdef bytes serial_BAC
    #     cdef int ifile
    #     f = open(fname,'wb')
    #     f.write(struct.pack('Q',self.nfiles))
    #     for ifile in range(self.nfiles):
    #         serial_BAC = self._dumps(ifile)
    #         f.write(struct.pack('Q',len(serial_BAC)))
    #         f.write(serial_BAC)
    #     f.close()

    # def load(self, fname):
    #     cdef np.uint64_t nfiles
    #     cdef np.uint64_t size_serial
    #     cdef bint read_flag = 1
    #     cdef bint irflag
    #     f = open(fname,'rb')
    #     nfiles, = struct.unpack('Q',f.read(struct.calcsize('Q')))
    #     if nfiles != self.nfiles:
    #         raise Exception("Number of bitmasks ({}) conflicts with number of files ({})".format(nfiles,self.nfiles))
    #     for ifile in range(nfiles):
    #         size_serial, = struct.unpack('Q',f.read(struct.calcsize('Q')))
    #         irflag = self._loads(ifile, f.read(size_serial))
    #         if irflag == 0: read_flag = 0
    #     f.close()
    #     return read_flag

    def __dealloc__(self):
        cdef ewah_bool_array *ewah_keys
        cdef ewah_bool_array *ewah_refn
        cdef ewah_bool_array *ewah_owns
        cdef ewah_map *ewah_coll
        for ifile in range(self.nfiles):
            ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
            ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
            ewah_owns = (<ewah_bool_array **> self.ewah_owns)[ifile]
            ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
            del ewah_keys
            del ewah_refn
            del ewah_owns
            del ewah_coll

    def print_info(self, ifile, prefix=''):
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(prefix,
                                                                     self._count_coarse(ifile),
                                                                     self._count_refined(ifile),
                                                                     self._count_total(ifile)))

cdef class BoolArrayCollection:

    def __cinit__(self):
        cdef ewah_bool_array *ewah_keys = new ewah_bool_array()
        cdef ewah_bool_array *ewah_refn = new ewah_bool_array()
        cdef ewah_bool_array *ewah_coar = new ewah_bool_array()
        cdef ewah_bool_array *ewah_owns = new ewah_bool_array()
        cdef ewah_map *ewah_coll = new ewah_map()
        self.ewah_keys = <void *> ewah_keys
        self.ewah_refn = <void *> ewah_refn
        self.ewah_coar = <void *> ewah_coar
        self.ewah_owns = <void *> ewah_owns
        self.ewah_coll = <void *> ewah_coll

    def reset(self):
        self.__dealloc__()
        self.__init__()

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
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                ewah_keys[0].set(i1)
                # self._set_coarse(i1)

    cdef void _set_refined_array(self, np.uint64_t i1, np.uint8_t[:] arr):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        cdef np.uint64_t i2
        for i2 in range(arr.shape[0]):
            if arr[i2] == 1:
                ewah_refn[0].set(i1)
                ewah_coll[0][i1].set(i2)
                # self._set_refined(i1, i2)

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

    cdef void _set_owns(self, np.uint64_t i1):
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> self.ewah_owns
        ewah_owns[0].set(i1)

    def set_owns(self, i1):
        self._set_owns(i1)

    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        if (ewah_keys[0].get(i1) == 0): return 0
        if (ewah_refn[0].get(i1) == 0) or (i2 == FLAG): 
            return 1
        return ewah_coll[0][i1].get(i2)

    def get(self, i1, i2 = FLAG):
        return self._get(i1, i2)

    cdef bint _get_coarse(self, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        return ewah_keys[0].get(i1)

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

    cdef int _count_owned(self):
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> self.ewah_owns
        cdef int out
        out = ewah_owns.numberOfOnes()
        return out

    def count_owned(self):
        return self._count_owned()

    cdef void _append(self, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_owns1 = <ewah_bool_array *> self.ewah_owns
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewah_bool_array *ewah_owns2 = <ewah_bool_array *> solf.ewah_owns
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
        # Owners
        ewah_owns1[0].logicalor(ewah_owns2[0], swap)
        ewah_owns1[0].swap(swap)
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
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # No intersection
        if ewah_keys1[0].intersects(ewah_keys2[0]) == 0:
            return 0
        # Intersection at coarse level
        ewah_keys1[0].logicalxor(ewah_refn1[0],ewah_coar1)
        ewah_keys2[0].logicalxor(ewah_refn2[0],ewah_coar2)
        if ewah_coar1.intersects(ewah_keys2[0]) == 1:
            return 1
        if ewah_coar2.intersects(ewah_keys1[0]) == 1:
            return 1
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
        return 0

    cdef void _logicalxor(self, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_owns1 = <ewah_bool_array *> self.ewah_owns
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewah_bool_array *ewah_owns2 = <ewah_bool_array *> solf.ewah_owns
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_bool_array *ewah_owns_out = <ewah_bool_array *> out.ewah_owns
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # Keys
        ewah_keys1[0].logicalxor(ewah_keys2[0],ewah_keys_out[0])
        # Refn
        ewah_refn1[0].logicalxor(ewah_refn2[0],ewah_refn_out[0])
        # Owns
        ewah_owns1[0].logicalxor(ewah_owns2[0],ewah_owns_out[0])
        # Coll
        it_map1 = ewah_coll1[0].begin()
        while (it_map1 != ewah_coll1[0].end()):
            mi1 = dereference(it_map1).first
            mi1_ewah1 = dereference(it_map1).second
            it_map2 = ewah_coll2[0].find(mi1)
            if it_map2 == ewah_coll2[0].end():
                ewah_coll_out[0][mi1] = mi1_ewah1
            else:
                mi1_ewah2 = dereference(it_map2).second
                mi1_ewah1.logicalxor(mi1_ewah2, swap)
                ewah_coll_out[0][mi1] = swap
            preincrement(it_map1)
        it_map2 = ewah_coll2[0].begin()
        while (it_map2 != ewah_coll2[0].end()):
            mi1 = dereference(it_map2).first
            mi1_ewah2 = dereference(it_map2).second
            it_map1 = ewah_coll1[0].find(mi1)
            if it_map1 == ewah_coll1[0].end():
                ewah_coll_out[0][mi1] = mi1_ewah2
            preincrement(it_map2)

    def logicalxor(self, solf, out):
        return self._logicalxor(solf, out)

    cdef void _logicaland(self, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_owns1 = <ewah_bool_array *> self.ewah_owns
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewah_bool_array *ewah_owns2 = <ewah_bool_array *> solf.ewah_owns
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_bool_array *ewah_owns_out = <ewah_bool_array *> out.ewah_owns
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
        cdef ewah_bool_array ewah_coar1, ewah_coar2
        # Keys
        ewah_keys1[0].logicaland(ewah_keys2[0],ewah_keys_out[0])
        # Refn
        ewah_refn1[0].logicaland(ewah_refn2[0],ewah_refn_out[0])
        # Owns
        ewah_owns1[0].logicaland(ewah_owns2[0],ewah_owns_out[0])
        # Coll
        if ewah_refn_out[0].numberOfOnes() > 0:
            it_map1 = ewah_coll1[0].begin()
            while (it_map1 != ewah_coll1[0].end()):
                mi1 = dereference(it_map1).first
                mi1_ewah1 = dereference(it_map1).second
                it_map2 = ewah_coll2[0].find(mi1)
                if it_map2 != ewah_coll2[0].end():
                    mi1_ewah2 = dereference(it_map2).second
                    mi1_ewah1.logicaland(mi1_ewah2, swap)
                    ewah_coll_out[0][mi1] = swap
                preincrement(it_map1)

    def logicaland(self, solf, out):
        return self._logicaland(solf, out)

    cdef bytes _dumps(self):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> self.ewah_owns
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
        # Owners
        ewah_owns[0].write(ss,1)
        # Return type cast python bytes string
        return <bytes>ss.str()

    def dumps(self):
        return self._dumps()

    cdef bint _loads(self, bytes s):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> self.ewah_owns
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
        # Owners
        if ss.eof() == 1: return 0
        ewah_owns[0].read(ss,1)
        return 1

    def loads(self, s):
        return self._loads(s)

    def save(self, fname):
        cdef bytes serial_BAC
        f = open(fname,'wb')
        serial_BAC = self._dumps()
        f.write(struct.pack('Q',len(serial_BAC)))
        f.write(serial_BAC)
        f.close()

    def load(self, fname):
        cdef np.uint64_t size_serial
        cdef bint flag_read
        f = open(fname,'rb')
        size_serial, = struct.unpack('Q',f.read(struct.calcsize('Q')))
        flag_read = self._loads(f.read(size_serial))
        f.close()
        return flag_read

    def __dealloc__(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> self.ewah_owns
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_keys
        del ewah_refn
        del ewah_coar
        del ewah_owns
        del ewah_coll

    def print_info(self, prefix=''):
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(prefix,
                                                                     self._count_coarse(),
                                                                     self._count_refined(),
                                                                     self._count_total()))

cdef class BoolArrayCollectionUncompressed:

    def __cinit__(self, np.uint64_t nele1, np.uint64_t nele2):
        self.nele1 = <int>nele1
        self.nele2 = <int>nele2
        cdef ewah_map *ewah_coll = new ewah_map()
        self.ewah_coll = <void *> ewah_coll
        cdef np.uint64_t i
        IF UncompressedFormat == 'MemoryView':
            self.ewah_keys = malloc(sizeof(bitarrtype)*nele1)
            self.ewah_refn = malloc(sizeof(bitarrtype)*nele1)
            self.ewah_owns = malloc(sizeof(bitarrtype)*nele1)
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:nele1]>self.ewah_keys
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:nele1]>self.ewah_refn
            cdef bitarrtype[:] ewah_owns = <bitarrtype[:nele1]>self.ewah_owns
            for i in range(nele1):
                ewah_keys[i] = 0
                ewah_refn[i] = 0
                ewah_owns[i] = 0
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *>malloc(sizeof(bitarrtype)*nele1)
            cdef bitarrtype *ewah_refn = <bitarrtype *>malloc(sizeof(bitarrtype)*nele1)
            cdef bitarrtype *ewah_owns = <bitarrtype *>malloc(sizeof(bitarrtype)*nele1)
            for i in range(nele1):
                ewah_keys[i] = 0
                ewah_refn[i] = 0
                ewah_owns[i] = 0
            self.ewah_keys = <void *> ewah_keys
            self.ewah_refn = <void *> ewah_refn
            self.ewah_owns = <void *> ewah_owns

    def reset(self):
        self.__dealloc__()
        self.__init__(self.nele1,self.nele2)

    cdef void _compress(self, BoolArrayCollection solf):
        cdef np.uint64_t i
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> solf.ewah_refn
        cdef ewah_bool_array *ewah_owns = <ewah_bool_array *> solf.ewah_owns
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] bool_keys = <bitarrtype[:self.nele1]>self.ewah_keys
            cdef bitarrtype[:] bool_refn = <bitarrtype[:self.nele1]>self.ewah_refn
            cdef bitarrtype[:] bool_owns = <bitarrtype[:self.nele1]>self.ewah_owns
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *bool_keys = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *bool_refn = <bitarrtype *> self.ewah_refn
            cdef bitarrtype *bool_owns = <bitarrtype *> self.ewah_owns
        for i in range(self.nele1):
            if bool_keys[i] == 1:
                ewah_keys[0].set(i)
            if bool_refn[i] == 1:
                ewah_refn[0].set(i)
            if bool_owns[i] == 1:
                ewah_owns[0].set(i)
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_map *ewah_coll2 = <ewah_map *> solf.ewah_coll
        ewah_coll2[0] = ewah_coll1[0]

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_keys[i1] = 1
        # Note the 0 here, for dereferencing
        if i2 != FLAG:
            ewah_refn[i1] = 1
            ewah_coll[0][i1].set(i2)

    cdef void _set_coarse(self, np.uint64_t i1):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        ewah_keys[i1] = 1

    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        ewah_refn[i1] = 1
        ewah_coll[0][i1].set(i2)

    cdef void _set_coarse_array(self, np.uint8_t[:] arr):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                ewah_keys[i1] = 1

    cdef void _set_coarse_array_ptr(self, np.uint8_t *arr):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
        ELIF UncompressedFormat == 'Pointer':
            # TODO: memcpy?
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(self.nele1):
            if arr[i1] == 1:
                ewah_keys[i1] = 1

    cdef void _set_refined_array(self, np.uint64_t i1, np.uint8_t[:] arr):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        cdef np.uint64_t i2
        for i2 in range(arr.shape[0]):
            if arr[i2] == 1:
                ewah_refn[i1] = 1
                ewah_coll[0][i1].set(i2)

    cdef void _set_refined_array_ptr(self, np.uint64_t i1, np.uint8_t *arr):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        cdef np.uint64_t i2
        for i2 in range(self.nele2):
            if arr[i2] == 1:
                ewah_refn[i1] = 1
                ewah_coll[0][i1].set(i2)

    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_coll[0][i1].set(i2)

    cdef void _set_refn(self, np.uint64_t i1):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        ewah_refn[i1] = 1

    cdef void _set_owns(self, np.uint64_t i1):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_owns = <bitarrtype[:self.nele1]>self.ewah_owns
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_owns = <bitarrtype *> self.ewah_owns
        ewah_owns[i1] = 1

    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        if ewah_keys[i1] == 0: return 0
        if (ewah_refn[i1] == 0) or (i2 == FLAG): 
            return 1
        return ewah_coll[0][i1].get(i2)

    cdef bint _get_coarse(self, np.uint64_t i1):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        return <bint>ewah_keys[i1]
        # if (ewah_keys[i1] == 0): return 0
        # return 1

    cdef bint _isref(self, np.uint64_t i):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        return <bint>ewah_refn[i]

    cdef int _count_total(self):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys = <bitarrtype[:self.nele1]>self.ewah_keys
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i
        cdef int out = 0
        for i in range(self.nele1):
            out += ewah_keys[i]
        return out

    cdef int _count_refined(self):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_refn = <bitarrtype[:self.nele1]>self.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef np.uint64_t i
        cdef int out = 0
        for i in range(self.nele1):
            out += ewah_refn[i]
        return out

    cdef int _count_owned(self):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_owns = <bitarrtype[:self.nele1]>self.ewah_owns
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_owns = <bitarrtype *> self.ewah_owns
        cdef np.uint64_t i
        cdef int out = 0
        for i in range(self.nele1):
            out += ewah_owns[i]
        return out

    cdef void _append(self, BoolArrayCollectionUncompressed solf):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys1 = <bitarrtype[:self.nele1]>self.ewah_keys
            cdef bitarrtype[:] ewah_refn1 = <bitarrtype[:self.nele1]>self.ewah_refn
            cdef bitarrtype[:] ewah_owns1 = <bitarrtype[:self.nele1]>self.ewah_owns
            cdef bitarrtype[:] ewah_keys2 = <bitarrtype[:solf.nele1]>solf.ewah_keys
            cdef bitarrtype[:] ewah_refn2 = <bitarrtype[:solf.nele1]>solf.ewah_refn
            cdef bitarrtype[:] ewah_owns2 = <bitarrtype[:solf.nele1]>solf.ewah_owns
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys1 = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *ewah_refn1 = <bitarrtype *> self.ewah_refn
            cdef bitarrtype *ewah_owns1 = <bitarrtype *> self.ewah_owns
            cdef bitarrtype *ewah_keys2 = <bitarrtype *> solf.ewah_keys
            cdef bitarrtype *ewah_refn2 = <bitarrtype *> solf.ewah_refn
            cdef bitarrtype *ewah_owns2 = <bitarrtype *> solf.ewah_owns
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array swap, mi1_ewah1, mi1_ewah2
        cdef np.uint64_t nrefn, mi1
        # TODO: Check if nele1 is equal?
        # Keys
        for mi1 in range(solf.nele1):
            if ewah_keys2[mi1] == 1:
                ewah_keys1[mi1] = 1
        # Refined
        for mi1 in range(solf.nele1):
            if ewah_refn2[mi1] == 1:
                ewah_refn1[mi1] = 1
        # Owners
        for mi1 in range(solf.nele1):
            if ewah_owns2[mi1] == 1:
                ewah_owns1[mi1] = 1
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

    cdef bint _intersects(self, BoolArrayCollectionUncompressed solf):
        IF UncompressedFormat == 'MemoryView':
            cdef bitarrtype[:] ewah_keys1 = <bitarrtype[:self.nele1]>self.ewah_keys
            cdef bitarrtype[:] ewah_refn1 = <bitarrtype[:self.nele1]>self.ewah_refn
            cdef bitarrtype[:] ewah_keys2 = <bitarrtype[:solf.nele1]>solf.ewah_keys
            cdef bitarrtype[:] ewah_refn2 = <bitarrtype[:solf.nele1]>solf.ewah_refn
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys1 = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *ewah_refn1 = <bitarrtype *> self.ewah_refn
            cdef bitarrtype *ewah_keys2 = <bitarrtype *> solf.ewah_keys
            cdef bitarrtype *ewah_refn2 = <bitarrtype *> solf.ewah_refn
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll1 = <cmap[np.uint64_t, ewah_bool_array] *> self.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array] *ewah_coll2 = <cmap[np.uint64_t, ewah_bool_array] *> solf.ewah_coll
        cdef cmap[np.uint64_t, ewah_bool_array].iterator it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        # No intersection
        for mi1 in range(self.nele1):
            if (ewah_keys1[mi1] == 1) and (ewah_keys2[mi1] == 1):
                break
        if (mi1 < self.nele1):
            return 0
        # Intersection at refined level
        for mi1 in range(self.nele1):
            if (ewah_refn1[mi1] == 1) and (ewah_refn2[mi1] == 1):
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
                break
        # Intersection at coarse level or refined inside coarse
        if mi1 == self.nele1:
            return 1
        return 0

    def __dealloc__(self):
        IF UncompressedFormat == 'MemoryView':
            free(self.ewah_keys)
            free(self.ewah_refn)
            free(self.ewah_owns)
        ELIF UncompressedFormat == 'Pointer':
            cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
            cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
            cdef bitarrtype *ewah_owns = <bitarrtype *> self.ewah_owns
            free(ewah_keys)
            free(ewah_refn)
            free(ewah_owns)
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_coll

    def print_info(self, prefix=''):
        cdef int nrefn = self._count_refined()
        cdef int nkeys = self._count_total()
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(prefix,
                                                                     nkeys - nrefn,
                                                                     nrefn,
                                                                     nkeys))

    

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

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
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

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
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

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
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

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
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

