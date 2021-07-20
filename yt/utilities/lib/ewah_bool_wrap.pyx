# distutils: language = c++
# distutils: include_dirs = LIB_DIR_EWAH
# distutils: extra_compile_args = CPP14_FLAG
"""
Wrapper for EWAH Bool Array: https://github.com/lemire/EWAHBoolArray



"""


import struct

from cython.operator cimport dereference, preincrement
from libc.stdlib cimport free, malloc, qsort
from libcpp.algorithm cimport sort
from libcpp.map cimport map as cmap

import numpy as np

cimport cython
cimport numpy as np

from yt.utilities.lib.geometry_utils cimport (
    morton_neighbors_coarse,
    morton_neighbors_refined,
)


cdef extern from "<algorithm>" namespace "std" nogil:
    Iter unique[Iter](Iter first, Iter last)

cdef np.uint64_t FLAG = ~(<np.uint64_t>0)
cdef np.uint64_t MAX_VECTOR_SIZE = <np.uint64_t>1e7

ctypedef cmap[np.uint64_t, ewah_bool_array] ewahmap
ctypedef cmap[np.uint64_t, ewah_bool_array].iterator ewahmap_it
ctypedef pair[np.uint64_t, ewah_bool_array] ewahmap_p

cdef class FileBitmasks:

    def __cinit__(self, np.uint32_t nfiles):
        cdef int i
        self.nfiles = nfiles
        self.ewah_keys = <ewah_bool_array **>malloc(nfiles*sizeof(ewah_bool_array*))
        self.ewah_refn = <ewah_bool_array **>malloc(nfiles*sizeof(ewah_bool_array*))
        self.ewah_coll = <ewah_map **>malloc(nfiles*sizeof(ewah_map*))
        for i in range(nfiles):
            self.ewah_keys[i] = new ewah_bool_array()
            self.ewah_refn[i] = new ewah_bool_array()
            self.ewah_coll[i] = new ewah_map()

    cdef void _reset(self):
        cdef np.int32_t ifile
        for ifile in range(self.nfiles):
            self.ewah_keys[ifile].reset()
            self.ewah_refn[ifile].reset()
            self.ewah_coll[ifile].clear()

    cdef bint _iseq(self, FileBitmasks solf):
        cdef np.int32_t ifile
        cdef ewah_bool_array* arr1
        cdef ewah_bool_array* arr2
        cdef ewahmap *map1
        cdef ewahmap *map2
        cdef ewahmap_p pair1, pair2
        cdef ewahmap_it it_map1, it_map2
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
            map1 = (<ewahmap **> self.ewah_coll)[ifile]
            map2 = (<ewahmap **> solf.ewah_coll)[ifile]
            for pair1 in map1[0]:
                it_map2 = map2[0].find(pair1.first)
                if it_map2 == map2[0].end():
                    return 0
                if pair1.second != dereference(it_map2).second:
                    return 0
            for pair2 in map2[0]:
                it_map1 = map1[0].find(pair2.first)
                if it_map1 == map1[0].end():
                    return 0
                if pair2.second != dereference(it_map1).second:
                    return 0
            # Match
            return 1

    def iseq(self, solf):
        return self._iseq(solf)

    cdef BoolArrayCollection _get_bitmask(self, np.uint32_t ifile):
        cdef BoolArrayCollection out = BoolArrayCollection()
        cdef ewah_bool_array **ewah_keys = <ewah_bool_array **>self.ewah_keys
        cdef ewah_bool_array **ewah_refn = <ewah_bool_array **>self.ewah_refn
        cdef ewah_map **ewah_coll = <ewah_map **>self.ewah_coll
        # This version actually copies arrays, which can be costly
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *>out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *>out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *>out.ewah_coll
        ewah_keys_out[0] = ewah_keys[ifile][0]
        ewah_refn_out[0] = ewah_refn[ifile][0]
        ewah_coll_out[0] = ewah_coll[ifile][0]
        # This version only copies pointers which can lead to deallocation of
        # the source when the copy is deleted.
        # out.ewah_keys = <void *>ewah_keys[ifile]
        # out.ewah_refn = <void *>ewah_refn[ifile]
        # out.ewah_coll = <void *>ewah_coll[ifile]
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
        cdef cmap[np.uint64_t, ewah_bool_array] map_keys, map_refn
        cdef cmap[np.uint64_t, ewah_bool_array]* coll_coll
        cdef cmap[np.uint64_t, ewah_bool_array]* map_bitmask
        coll_refn = <ewah_bool_array*> coll.ewah_refn
        if coll_refn[0].numberOfOnes() == 0:
            if verbose == 1:
                print("{: 10d}/{: 10d} collisions at refined refinement. ({: 10.5f}%)".format(0,0,0))
            return (0,0)
        coll_coll = <cmap[np.uint64_t, ewah_bool_array]*> coll.ewah_coll
        for ifile in range(self.nfiles):
            map_bitmask = (<cmap[np.uint64_t, ewah_bool_array]**> self.ewah_coll)[ifile]
            for it_mi1 in map_bitmask[0]:
                mi1 = it_mi1.first
                iarr = it_mi1.second
                map_keys[mi1].logicaland(iarr, arr_two)
                map_keys[mi1].logicalor(iarr, arr_swap)
                map_keys[mi1].swap(arr_swap)
                map_refn[mi1].logicalor(arr_two, arr_swap)
                map_refn[mi1].swap(arr_swap)
        coll_coll[0] = map_refn
        # Count
        cdef int nc, nm
        nc = 0
        nm = 0
        for it_mi1 in map_refn:
            mi1 = it_mi1.first
            iarr = it_mi1.second
            nc += iarr.numberOfOnes()
            iarr = map_keys[mi1]
            nm += iarr.numberOfOnes()
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_refined_index_array(self, np.uint32_t ifile, np.int64_t nsub_mi,
                                       np.ndarray[np.uint64_t, ndim=1] sub_mi1,
                                       np.ndarray[np.uint64_t, ndim=1] sub_mi2):
        cdef np.ndarray[np.int64_t, ndim=1] ind = np.lexsort((sub_mi2[:nsub_mi],
                                                              sub_mi1[:nsub_mi]))
        cdef np.int64_t i, p
        cdef BoolArrayCollection temp
        if self._count_refined(ifile) == 0:
            # Add to file bitmask in order
            for i in range(nsub_mi):
                p = ind[i]
                self._set_refined(ifile, sub_mi1[p], sub_mi2[p])
        else:
            # Add to dummy bitmask in order, then combine
            temp = BoolArrayCollection()
            for i in range(nsub_mi):
                p = ind[i]
                temp._set_coarse(sub_mi1[p])
                temp._set_refined(sub_mi1[p], sub_mi2[p])
                self._append(ifile, temp)

    cdef void _set_map(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        ewah_coll[0][i1].set(i2)

    cdef void _set_refn(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        ewah_refn[0].set(i1)

    cdef bint _get(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        if (ewah_keys[0].get(i1) == 0): return 0
        if (i2 == FLAG) or (ewah_refn[0].get(i1) == 0):
            return 1
        return ewah_coll[0][i1].get(i2)

    cdef bint _get_coarse(self, np.uint32_t ifile, np.uint64_t i1):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        return ewah_keys[0].get(i1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _get_coarse_array(self, np.uint32_t ifile, np.uint64_t imax,
                                np.uint8_t[:] arr) except *:
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_keys[0].begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_keys[0].end())
        cdef np.uint64_t iset
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            if iset >= imax:
                raise IndexError("Index {} exceedes max {}.".format(iset, imax))
            arr[iset] = 1
            preincrement(iter_set[0])

    cdef bint _isref(self, np.uint32_t ifile, np.uint64_t i):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        return ewah_refn[0].get(i)

    def count_coarse(self, ifile):
        return self._count_coarse(ifile)

    def count_total(self, ifile):
        return self._count_total(ifile)

    def count_refined(self, ifile):
        return self._count_refined(ifile)

    cdef np.uint64_t _count_coarse(self, np.uint32_t ifile):
        return self._count_total(ifile) - self._count_refined(ifile)

    cdef np.uint64_t _count_total(self, np.uint32_t ifile):
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef np.uint64_t out = ewah_keys[0].numberOfOnes()
        return out

    cdef np.uint64_t _count_refined(self, np.uint32_t ifile):
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef np.uint64_t out = ewah_refn[0].numberOfOnes()
        return out

    def append(self, np.uint32_t ifile, BoolArrayCollection solf):
        if solf is None: return
        self._append(ifile, solf)

    cdef void _append(self, np.uint32_t ifile, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn1 = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll1 = (<ewah_map **> self.ewah_coll)[ifile]
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array swap, mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
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
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
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
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
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
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
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

    cdef void _select_contaminated(self, np.uint32_t ifile,
                                   BoolArrayCollection mask, np.uint8_t[:] out,
                                   np.uint8_t[:] secondary_files,
                                   BoolArrayCollection mask2 = None):
        # Fill mask at indices owned by this file that are also contaminated by
        # other files.
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_bool_array ewah_mask
        cdef ewah_bool_array *ewah_mask1
        cdef ewah_bool_array *ewah_mask2
        cdef ewah_bool_array ewah_slct
        cdef ewah_bool_array *ewah_file
        cdef np.uint64_t iset
        # Merge masks as necessary
        if mask2 is None:
            ewah_mask = (<ewah_bool_array *> mask.ewah_keys)[0]
        else:
            ewah_mask1 = <ewah_bool_array *> mask.ewah_keys
            ewah_mask2 = <ewah_bool_array *> mask2.ewah_keys
            ewah_mask1[0].logicalor(ewah_mask2[0],ewah_mask)
        # Get just refined cells owned by this file
        ewah_mask.logicaland(ewah_refn[0], ewah_slct)
        # Set array values
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct.begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct.end())
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            out[iset] = 1
            preincrement(iter_set[0])
        # Find files that intersect this one
        cdef np.uint32_t isfile
        for isfile in range(self.nfiles):
            if isfile == ifile: continue
            ewah_file = (<ewah_bool_array **> self.ewah_keys)[isfile]
            if ewah_slct.intersects(ewah_file[0]) == 1:
                secondary_files[isfile] = 1

    cdef void _select_uncontaminated(self, np.uint32_t ifile,
                                     BoolArrayCollection mask, np.uint8_t[:] out,
                                     BoolArrayCollection mask2 = None):
        # Fill mask at indices that are owned by this file and no other.
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_bool_array ewah_mask
        cdef ewah_bool_array *ewah_mask1
        cdef ewah_bool_array *ewah_mask2
        cdef ewah_bool_array ewah_slct
        cdef ewah_bool_array ewah_coar
        cdef np.uint64_t iset
        # Merge masks if necessary
        if mask2 is None:
            ewah_mask = (<ewah_bool_array *> mask.ewah_keys)[0]
        else:
            ewah_mask1 = <ewah_bool_array *> mask.ewah_keys
            ewah_mask2 = <ewah_bool_array *> mask2.ewah_keys
            ewah_mask1[0].logicalor(ewah_mask2[0],ewah_mask)
        # Get coarse cells owned by this file
        ewah_keys[0].logicalxor(ewah_refn[0],ewah_coar)
        ewah_coar.logicaland(ewah_mask,ewah_slct)
        # Set array elements
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct.begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct.end())
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            out[iset] = 1
            preincrement(iter_set[0])

    cdef bytes _dumps(self, np.uint32_t ifile):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        cdef ewahmap_it it_map
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

    cdef bint _loads(self, np.uint32_t ifile, bytes s):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
        cdef ewah_bool_array *ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
        cdef ewah_map *ewah_coll = (<ewah_map **> self.ewah_coll)[ifile]
        cdef np.uint64_t nrefn, mi1
        nrefn = mi1 = 0
        # Write string to string stream
        if len(s) == 0: return 1
        ss.write(s, len(s))
        # Read keys and refinement arrays
        ewah_keys[0].read(ss,1)
        if ss.eof(): return 1
        ewah_refn[0].read(ss,1)
        # Read and check number of refined cells
        ss.read(<char *> (&nrefn), sizeof(nrefn))
        if nrefn != ewah_refn[0].numberOfOnes():
            raise Exception("Error in read. File indicates {} refinements, but bool array has {}.".format(nrefn,ewah_refn[0].numberOfOnes()))
        # Loop over refined cells
        for _ in range(nrefn):
            ss.read(<char *> (&mi1), sizeof(mi1))
            if ss.eof(): return 1
            ewah_coll[0][mi1].read(ss,1)
            # or...
            #mi1_ewah.read(ss,1)
            #ewah_coll[0][mi1].swap(mi1_ewah)
        return 1

    cdef bint _check(self):
        cdef np.uint32_t ifile
        cdef ewah_bool_array *ewah_keys
        cdef ewah_bool_array *ewah_refn
        cdef ewah_bool_array tmp1, tmp2
        cdef np.uint64_t nchk
        cdef str msg
        # Check individual files
        for ifile in range(self.nfiles):
            ewah_keys = (<ewah_bool_array **> self.ewah_keys)[ifile]
            ewah_refn = (<ewah_bool_array **> self.ewah_refn)[ifile]
            # Check that there are not any refn that are not keys
            ewah_keys[0].logicalxor(ewah_refn[0], tmp1)
            ewah_refn[0].logicaland(tmp1, tmp2)
            nchk = tmp2.numberOfOnes()
            if nchk > 0:
                msg = "File {}: There are {} refined cells that are not set on coarse level.".format(ifile,nchk)
                print(msg)
                return 0
                # raise Exception(msg)
        return 1

    def check(self):
        return self._check()

    def __dealloc__(self):
        for ifile in range(self.nfiles):
            del self.ewah_keys[ifile]
            del self.ewah_refn[ifile]
            del self.ewah_coll[ifile]

    def print_info(self, ifile, prefix=''):
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(
            prefix,
            self._count_coarse(ifile),
            self._count_refined(ifile),
            self._count_total(ifile)))

cdef class BoolArrayCollection:

    def __cinit__(self):
        self.ewah_keys = new ewah_bool_array()
        self.ewah_refn = new ewah_bool_array()
        self.ewah_coar = new ewah_bool_array()
        self.ewah_coll = new ewah_map()

    cdef void _reset(self):
        self.ewah_keys[0].reset()
        self.ewah_refn[0].reset()
        self.ewah_coar[0].reset()
        self.ewah_coll[0].clear()

    cdef int _richcmp(self, BoolArrayCollection solf, int op) except -1:

        cdef ewah_bool_array *arr1
        cdef ewah_bool_array *arr2
        cdef ewahmap *map1
        cdef ewahmap *map2
        cdef ewahmap_it it_map1, it_map2
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
            map1 = <ewahmap *> self.ewah_coll
            map2 = <ewahmap *> solf.ewah_coll
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    def set_from(self, np.uint64_t[:] ids):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef np.uint64_t i
        cdef np.uint64_t last = 0
        for i in range(ids.shape[0]):
            if ids[i] < last:
                raise RuntimeError
            self._set(ids[i])
            last = ids[i]
        print("Set from %s array and ended up with %s bytes" % (
            ids.size, ewah_keys[0].sizeInBytes()))

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_coarse_array(self, np.uint8_t[:] arr):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                ewah_keys[0].set(i1)
                # self._set_coarse(i1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _get_coarse_array(self, np.uint64_t imax, np.uint8_t[:] arr) except *:
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_keys[0].begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_keys[0].end())
        cdef np.uint64_t iset
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            if iset >= imax:
                raise IndexError("Index {} exceedes max {}.".format(iset, imax))
            arr[iset] = 1
            preincrement(iter_set[0])

    def get_coarse_array(self, imax, arr):
        return self._get_coarse_array(imax, arr)

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

    cdef np.uint64_t _count_total(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef np.uint64_t out = ewah_keys.numberOfOnes()
        return out

    def count_total(self):
        return self._count_total()

    cdef np.uint64_t _count_refined(self):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef np.uint64_t out = ewah_refn.numberOfOnes()
        return out

    def count_refined(self):
        return self._count_refined()

    cdef np.uint64_t _count_coarse(self):
        self._ewah_coarse()
        cdef ewah_bool_array *ewah_coar = <ewah_bool_array *> self.ewah_coar
        cdef np.uint64_t out = ewah_coar.numberOfOnes()
        return out

    def count_coarse(self):
        return self._count_coarse()

    cdef void _logicalor(self, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = self.ewah_refn
        cdef ewahmap *ewah_coll1 = self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = solf.ewah_refn
        cdef ewahmap *ewah_coll2 = solf.ewah_coll
        cdef ewah_bool_array *ewah_keys3 = out.ewah_keys
        cdef ewah_bool_array *ewah_refn3 = out.ewah_refn
        cdef ewahmap *ewah_coll3 = out.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        # Keys
        ewah_keys1[0].logicalor(ewah_keys2[0], ewah_keys3[0])
        # Refined
        ewah_refn1[0].logicalor(ewah_refn2[0], ewah_refn3[0])
        # Map
        it_map1 = ewah_coll1[0].begin()
        while it_map1 != ewah_coll1[0].end():
            mi1 = dereference(it_map1).first
            mi1_ewah1 = dereference(it_map1).second
            ewah_coll3[0][mi1] = mi1_ewah1
            preincrement(it_map1)
        it_map2 = ewah_coll2[0].begin()
        while it_map2 != ewah_coll2[0].end():
            mi1 = dereference(it_map2).first
            mi1_ewah2 = dereference(it_map2).second
            it_map1 = ewah_coll1[0].find(mi1)
            if it_map1 != ewah_coll1[0].end():
                mi1_ewah1 = dereference(it_map1).second
                mi1_ewah1.logicalor(mi1_ewah2, ewah_coll3[0][mi1])
            else:
                ewah_coll3[0][mi1] = mi1_ewah2
            preincrement(it_map2)

    cdef void _append(self, BoolArrayCollection solf):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewahmap *ewah_coll1 = <ewahmap *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array swap, mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
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
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewahmap *ewah_coll1 = <ewahmap *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
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
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
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

    def logicalxor(self, solf, out):
        return self._logicalxor(solf, out)

    cdef void _logicaland(self, BoolArrayCollection solf, BoolArrayCollection out):
        cdef ewah_bool_array *ewah_keys1 = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn1 = <ewah_bool_array *> self.ewah_refn
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_bool_array *ewah_keys2 = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn2 = <ewah_bool_array *> solf.ewah_refn
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewah_bool_array *ewah_keys_out = <ewah_bool_array *> out.ewah_keys
        cdef ewah_bool_array *ewah_refn_out = <ewah_bool_array *> out.ewah_refn
        cdef ewah_map *ewah_coll_out = <ewah_map *> out.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2, swap
        cdef np.uint64_t mi1
        # Keys
        ewah_keys1[0].logicaland(ewah_keys2[0],ewah_keys_out[0])
        # Refn
        ewah_refn1[0].logicaland(ewah_refn2[0],ewah_refn_out[0])
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

    cdef void _select_contaminated(self, BoolArrayCollection mask, np.uint8_t[:] out,
                                   BoolArrayCollection mask2 = None):
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array ewah_mask
        cdef ewah_bool_array *ewah_mask1
        cdef ewah_bool_array *ewah_mask2
        if mask2 is None:
            ewah_mask = (<ewah_bool_array *> mask.ewah_keys)[0]
        else:
            ewah_mask1 = <ewah_bool_array *> mask.ewah_keys
            ewah_mask2 = <ewah_bool_array *> mask2.ewah_keys
            ewah_mask1[0].logicalor(ewah_mask2[0],ewah_mask)
        cdef ewah_bool_array ewah_slct
        ewah_refn[0].logicaland(ewah_mask,ewah_slct)
        cdef np.uint64_t iset
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct.begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct.end())
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            out[iset] = 1
            preincrement(iter_set[0])

    cdef void _select_uncontaminated(self, BoolArrayCollection mask, np.uint8_t[:] out,
                                     BoolArrayCollection mask2 = None):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array ewah_mask
        cdef ewah_bool_array *ewah_mask1
        cdef ewah_bool_array *ewah_mask2
        if mask2 is None:
            ewah_mask = (<ewah_bool_array *> mask.ewah_keys)[0]
        else:
            ewah_mask1 = <ewah_bool_array *> mask.ewah_keys
            ewah_mask2 = <ewah_bool_array *> mask2.ewah_keys
            ewah_mask1[0].logicalor(ewah_mask2[0],ewah_mask)
        cdef ewah_bool_array ewah_slct
        cdef ewah_bool_array ewah_coar
        ewah_keys[0].logicalxor(ewah_refn[0],ewah_coar)
        ewah_coar.logicaland(ewah_mask,ewah_slct)
        cdef np.uint64_t iset
        cdef ewah_bool_iterator *iter_set = new ewah_bool_iterator(ewah_slct.begin())
        cdef ewah_bool_iterator *iter_end = new ewah_bool_iterator(ewah_slct.end())
        while iter_set[0] != iter_end[0]:
            iset = dereference(iter_set[0])
            out[iset] = 1
            preincrement(iter_set[0])

    cdef void _get_ghost_zones(self, int ngz, int order1, int order2,
                               bint periodicity[3], BoolArrayCollection out_ewah,
                               bint coarse_ghosts = 0):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewahmap *ewah_coll = <ewahmap *> self.ewah_coll
        cdef ewah_bool_iterator *iter_set1 = new ewah_bool_iterator(ewah_keys.begin())
        cdef ewah_bool_iterator *iter_end1 = new ewah_bool_iterator(ewah_keys.end())
        cdef ewah_bool_iterator *iter_set2
        cdef ewah_bool_iterator *iter_end2
        cdef np.uint64_t max_index1 = <np.uint64_t>(1 << order1)
        cdef np.uint64_t max_index2 = <np.uint64_t>(1 << order2)
        cdef np.uint64_t nele1 = <np.uint64_t>(max_index1**3)
        cdef np.uint64_t nele2 = <np.uint64_t>(max_index2**3)
        cdef BoolArrayCollectionUncompressed temp_bool = BoolArrayCollectionUncompressed(nele1, nele2)
        cdef BoolArrayCollectionUncompressed out_bool = BoolArrayCollectionUncompressed(nele1, nele2)
        cdef np.uint64_t mi1, mi2, mi1_n, mi2_n
        cdef np.uint32_t ntot, i
        cdef void* pointers[7]
        pointers[0] = malloc( sizeof(np.int32_t) * (2*ngz+1)*3)
        pointers[1] = malloc( sizeof(np.uint64_t) * (2*ngz+1)*3)
        pointers[2] = malloc( sizeof(np.uint64_t) * (2*ngz+1)*3)
        pointers[3] = malloc( sizeof(np.uint64_t) * (2*ngz+1)**3)
        pointers[4] = malloc( sizeof(np.uint64_t) * (2*ngz+1)**3)
        pointers[5] = malloc( sizeof(np.uint8_t) * nele1)
        pointers[6] = malloc( sizeof(np.uint8_t) * nele2)
        cdef np.uint32_t[:,:] index = <np.uint32_t[:2*ngz+1,:3]> pointers[0]
        cdef np.uint64_t[:,:] ind1_n = <np.uint64_t[:2*ngz+1,:3]> pointers[1]
        cdef np.uint64_t[:,:] ind2_n = <np.uint64_t[:2*ngz+1,:3]> pointers[2]
        cdef np.uint64_t[:] neighbor_list1 = <np.uint64_t[:((2*ngz+1)**3)]> pointers[3]
        cdef np.uint64_t[:] neighbor_list2 = <np.uint64_t[:((2*ngz+1)**3)]> pointers[4]
        cdef np.uint8_t *bool_keys = <np.uint8_t *> pointers[5]
        cdef np.uint8_t *bool_coll = <np.uint8_t *> pointers[6]
        cdef SparseUnorderedRefinedBitmaskSet list_coll = SparseUnorderedRefinedBitmaskSet()
        for i in range(nele1):
            bool_keys[i] = 0
        while iter_set1[0] != iter_end1[0]:
            mi1 = dereference(iter_set1[0])
            if (coarse_ghosts == 1) or (ewah_refn[0].get(mi1) == 0):
                # Coarse neighbors
                ntot = morton_neighbors_coarse(mi1, max_index1, periodicity, ngz,
                                               index, ind1_n, neighbor_list1)
                for i in range(ntot):
                    mi1_n = neighbor_list1[i]
                    if ewah_keys[0].get(mi1_n) == 0:
                        bool_keys[mi1_n] = 1
            else:
                for i in range(nele2):
                    bool_coll[i] = 0
                # Refined neighbors
                iter_set2 = new ewah_bool_iterator(ewah_coll[0][mi1].begin())
                iter_end2 = new ewah_bool_iterator(ewah_coll[0][mi1].end())
                while iter_set2[0] != iter_end2[0]:
                    mi2 = dereference(iter_set2[0])
                    ntot = morton_neighbors_refined(mi1, mi2,
                                                    max_index1, max_index2,
                                                    periodicity, ngz, index,
                                                    ind1_n, ind2_n,
                                                    neighbor_list1,
                                                    neighbor_list2)
                    for i in range(ntot):
                        mi1_n = neighbor_list1[i]
                        mi2_n = neighbor_list2[i]
                        if mi1_n == mi1:
                            if ewah_coll[0][mi1].get(mi2_n) == 0:
                                bool_keys[mi1_n] = 1
                                bool_coll[mi2_n] = 1
                        else:
                            if ewah_refn[0].get(mi1_n) == 1:
                                if ewah_coll[0][mi1_n].get(mi2_n) == 0:
                                    bool_keys[mi1_n] = 1
                                    list_coll._set(mi1_n, mi2_n)
                            else:
                                if ewah_keys[0].get(mi1_n) == 0:
                                    bool_keys[mi1_n] = 1
                    preincrement(iter_set2[0])
                # Add to running list
                temp_bool._set_refined_array_ptr(mi1, bool_coll)
            preincrement(iter_set1[0])
        # Set keys
        out_bool._set_coarse_array_ptr(bool_keys)
        list_coll._fill_bool(out_bool)
        out_bool._append(temp_bool)
        out_bool._compress(out_ewah)
        # Free things
        for i in range(7):
            free(pointers[i])

    cdef bytes _dumps(self):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewahmap *ewah_coll = <ewahmap *> self.ewah_coll
        cdef ewahmap_it it_map
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

    cdef bint _loads(self, bytes s):
        # TODO: write word size
        cdef sstream ss
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewahmap *ewah_coll = <ewahmap *> self.ewah_coll
        cdef np.uint64_t nrefn, mi1
        nrefn = mi1 = 0
        # Write string to string stream
        if len(s) == 0: return 1
        ss.write(s, len(s))
        # Read keys and refinement arrays
        if ss.eof(): return 1
        ewah_keys[0].read(ss,1)
        if ss.eof(): return 1
        ewah_refn[0].read(ss,1)
        # Read and check number of refined cells
        if ss.eof(): return 1
        ss.read(<char *> (&nrefn), sizeof(nrefn))
        if nrefn != ewah_refn[0].numberOfOnes():
            raise Exception("Error in read. File indicates {} refinements, but bool array has {}.".format(nrefn,ewah_refn[0].numberOfOnes()))
        # Loop over refined cells
        for _ in range(nrefn):
            ss.read(<char *> (&mi1), sizeof(mi1))
            if ss.eof():
                # A brief note about why we do this!
                # In previous versions of the EWAH code, which were more
                # susceptible to issues with differences in sizes of size_t
                # etc, the ewah_coll.read would use instance variables as
                # destinations; these were initialized to zero.  In recent
                # versions, it uses (uninitialized) temporary variables.  We
                # were passing in streams that were already at EOF - so the
                # uninitialized memory would not be written to, and it would
                # retain the previous values, which would invariably be really
                # really big!  So we do a check for EOF here to make sure we're
                # not up to no good.
                break
            ewah_coll[0][mi1].read(ss,1)
            # or...
            #mi1_ewah.read(ss,1)
            #ewah_coll[0][mi1].swap(mi1_ewah)
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

    cdef bint _check(self):
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> self.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> self.ewah_refn
        cdef ewah_bool_array tmp1, tmp2
        cdef np.uint64_t nchk
        cdef str msg
        # Check that there are not any refn that are not keys
        ewah_keys[0].logicalxor(ewah_refn[0], tmp1)
        ewah_refn[0].logicaland(tmp1, tmp2)
        nchk = tmp2.numberOfOnes()
        if nchk > 0:
            msg = "There are {} refined cells that are not set on coarse level.".format(nchk)
            print(msg)
            return 0
            # raise Exception(msg)
        return 1

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

cdef class BoolArrayCollectionUncompressed:

    def __cinit__(self, np.uint64_t nele1, np.uint64_t nele2):
        self.nele1 = <int>nele1
        self.nele2 = <int>nele2
        self.ewah_coll = new ewah_map()
        cdef np.uint64_t i
        self.ewah_keys = <bitarrtype *>malloc(sizeof(bitarrtype)*nele1)
        self.ewah_refn = <bitarrtype *>malloc(sizeof(bitarrtype)*nele1)
        for i in range(nele1):
            self.ewah_keys[i] = 0
            self.ewah_refn[i] = 0

    def reset(self):
        self.__dealloc__()
        self.__init__(self.nele1,self.nele2)

    cdef void _compress(self, BoolArrayCollection solf):
        cdef np.uint64_t i
        cdef ewah_bool_array *ewah_keys = <ewah_bool_array *> solf.ewah_keys
        cdef ewah_bool_array *ewah_refn = <ewah_bool_array *> solf.ewah_refn
        cdef bitarrtype *bool_keys = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *bool_refn = <bitarrtype *> self.ewah_refn
        for i in range(self.nele1):
            if bool_keys[i] == 1:
                ewah_keys[0].set(i)
            if bool_refn[i] == 1:
                ewah_refn[0].set(i)
        cdef ewah_map *ewah_coll1 = <ewah_map *> self.ewah_coll
        cdef ewah_map *ewah_coll2 = <ewah_map *> solf.ewah_coll
        ewah_coll2[0] = ewah_coll1[0]

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_keys[i1] = 1
        # Note the 0 here, for dereferencing
        if i2 != FLAG:
            ewah_refn[i1] = 1
            ewah_coll[0][i1].set(i2)

    cdef void _set_coarse(self, np.uint64_t i1):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        ewah_keys[i1] = 1

    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        ewah_refn[i1] = 1
        ewah_coll[0][i1].set(i2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_coarse_array(self, np.uint8_t[:] arr):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(arr.shape[0]):
            if arr[i1] == 1:
                ewah_keys[i1] = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_coarse_array_ptr(self, np.uint8_t *arr):
        # TODO: memcpy?
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i1
        for i1 in range(self.nele1):
            if arr[i1] == 1:
                ewah_keys[i1] = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_refined_array(self, np.uint64_t i1, np.uint8_t[:] arr):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        cdef np.uint64_t i2
        for i2 in range(arr.shape[0]):
            if arr[i2] == 1:
                ewah_refn[i1] = 1
                ewah_coll[0][i1].set(i2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    cdef void _set_refined_array_ptr(self, np.uint64_t i1, np.uint8_t *arr):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        cdef np.uint64_t i2
        cdef ewah_bool_array *barr = &ewah_coll[0][i1]
        for i2 in range(self.nele2):
            if arr[i2] == 1:
                ewah_refn[i1] = 1
                barr.set(i2)

    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2):
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        ewah_coll[0][i1].set(i2)

    cdef void _set_refn(self, np.uint64_t i1):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        ewah_refn[i1] = 1

    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2 = FLAG):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        # Note the 0 here, for dereferencing
        if ewah_keys[i1] == 0: return 0
        if (ewah_refn[i1] == 0) or (i2 == FLAG):
            return 1
        return ewah_coll[0][i1].get(i2)

    cdef bint _get_coarse(self, np.uint64_t i1):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        return <bint>ewah_keys[i1]
        # if (ewah_keys[i1] == 0): return 0
        # return 1

    cdef bint _isref(self, np.uint64_t i):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        return <bint>ewah_refn[i]

    cdef np.uint64_t _count_total(self):
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef np.uint64_t i
        cdef np.uint64_t out = 0
        for i in range(self.nele1):
            out += ewah_keys[i]
        return out

    cdef np.uint64_t _count_refined(self):
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        cdef np.uint64_t i
        cdef np.uint64_t out = 0
        for i in range(self.nele1):
            out += ewah_refn[i]
        return out

    cdef void _append(self, BoolArrayCollectionUncompressed solf):
        cdef bitarrtype *ewah_keys1 = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *ewah_refn1 = <bitarrtype *> self.ewah_refn
        cdef bitarrtype *ewah_keys2 = <bitarrtype *> solf.ewah_keys
        cdef bitarrtype *ewah_refn2 = <bitarrtype *> solf.ewah_refn
        cdef ewahmap *ewah_coll1 = <ewahmap *> self.ewah_coll
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array swap, mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        # TODO: Check if nele1 is equal?
        # Keys
        for mi1 in range(solf.nele1):
            if ewah_keys2[mi1] == 1:
                ewah_keys1[mi1] = 1
        # Refined
        for mi1 in range(solf.nele1):
            if ewah_refn2[mi1] == 1:
                ewah_refn1[mi1] = 1
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
        cdef bitarrtype *ewah_keys1 = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *ewah_refn1 = <bitarrtype *> self.ewah_refn
        cdef bitarrtype *ewah_keys2 = <bitarrtype *> solf.ewah_keys
        cdef bitarrtype *ewah_refn2 = <bitarrtype *> solf.ewah_refn
        cdef ewahmap *ewah_coll1 = <ewahmap *> self.ewah_coll
        cdef ewahmap *ewah_coll2 = <ewahmap *> solf.ewah_coll
        cdef ewahmap_it it_map1, it_map2
        cdef ewah_bool_array mi1_ewah1, mi1_ewah2
        cdef np.uint64_t mi1
        # No intersection
        for mi1 in range(self.nele1):
            if (ewah_keys1[mi1] == 1) and (ewah_keys2[mi1] == 1):
                break
        if (mi1 < self.nele1):
            return 0
        mi1 = self.nele1 # This is to get rid of a warning
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
        cdef bitarrtype *ewah_keys = <bitarrtype *> self.ewah_keys
        cdef bitarrtype *ewah_refn = <bitarrtype *> self.ewah_refn
        free(ewah_keys)
        free(ewah_refn)
        cdef ewah_map *ewah_coll = <ewah_map *> self.ewah_coll
        del ewah_coll

    def print_info(self, prefix=''):
        cdef np.uint64_t nrefn = self._count_refined()
        cdef np.uint64_t nkeys = self._count_total()
        print("{}{: 8d} coarse, {: 8d} refined, {: 8d} total".format(prefix,
                                                                     nkeys - nrefn,
                                                                     nrefn,
                                                                     nkeys))



# Vector version
cdef class SparseUnorderedBitmaskVector:
    def __cinit__(self):
        self.total = 0

    cdef void _set(self, np.uint64_t ind):
        self.entries.push_back(ind)
        self.total += 1

    def set(self, ind):
        self._set(ind)

    cdef void _fill(self, np.uint8_t[:] mask):
        cdef np.uint64_t i, ind
        for i in range(self.entries.size()):
            ind = self.entries[i]
            mask[ind] = 1

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        self._remove_duplicates()
        cdef np.uint64_t i, ind
        for i in range(self.entries.size()):
            ind = self.entries[i]
            mm._set_coarse(ind)

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
        self._remove_duplicates()
        cdef np.uint64_t i, ind
        for i in range(self.entries.size()):
            ind = self.entries[i]
            mm._set_coarse(ind)

    cdef void _reset(self):
        self.entries.erase(self.entries.begin(), self.entries.end())
        self.total = 0

    cdef to_array(self):
        self._remove_duplicates()
        cdef np.ndarray[np.uint64_t, ndim=1] rv
        rv = np.empty(self.entries.size(), dtype='uint64')
        for i in range(self.entries.size()):
            rv[i] = self.entries[i]
        return rv

    cdef void _remove_duplicates(self):
        cdef vector[np.uint64_t].iterator last
        sort(self.entries.begin(), self.entries.end())
        last = unique(self.entries.begin(), self.entries.end())
        self.entries.erase(last, self.entries.end())

    cdef void _prune(self):
        if self.total > MAX_VECTOR_SIZE:
            self._remove_duplicates()
            self.total = 0

    def __dealloc__(self):
        self.entries.clear()

# Set version
cdef class SparseUnorderedBitmaskSet:
    cdef void _set(self, np.uint64_t ind):
        self.entries.insert(ind)

    def set(self, ind):
        self._set(ind)

    cdef void _fill(self, np.uint8_t[:] mask):
        for it in self.entries:
            mask[it] = 1

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        for it in self.entries:
            mm._set_coarse(it)

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
        for it in self.entries:
            mm._set_coarse(it)

    cdef void _reset(self):
        self.entries.clear()

    cdef to_array(self):
        cdef np.uint64_t ind
        cdef np.ndarray[np.uint64_t, ndim=1] rv
        cdef cset[np.uint64_t].iterator it
        rv = np.empty(self.entries.size(), dtype='uint64')
        it = self.entries.begin()
        i = 0
        while it != self.entries.end():
            ind = dereference(it)
            rv[i] = ind
            preincrement(it)
            i += 1
        return rv

    def __dealloc__(self):
        self.entries.clear()

# vector version
cdef class SparseUnorderedRefinedBitmaskVector:
    def __cinit__(self):
        self.total = 0

    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2):
        cdef ind_pair ind
        ind.first = ind1
        ind.second = ind2
        self.entries.push_back(ind)
        self.total += 1

    def set(self, ind1, ind2):
        self._set(ind1, ind2)

    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:] mask2):
        for it in self.entries:
            mask1[it.first] = mask2[it.second] = 1

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        self._remove_duplicates()
        for it in self.entries:
            mm._set_refined(it.first, it.second)

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
        self._remove_duplicates()
        for it in self.entries:
            mm._set_refined(it.first, it.second)

    cdef void _reset(self):
        self.entries.erase(self.entries.begin(), self.entries.end())
        self.total = 0

    cdef to_array(self):
        cdef np.uint64_t i
        cdef np.ndarray[np.uint64_t, ndim=2] rv
        self._remove_duplicates()
        rv = np.empty((self.entries.size(),2),dtype='uint64')
        i = 0
        for it in self.entries:
            rv[i,0] = it.first
            rv[i,1] = it.second
            i += 1
        return rv

    cdef void _remove_duplicates(self):
        cdef vector[ind_pair].iterator last
        sort(self.entries.begin(), self.entries.end())
        last = unique(self.entries.begin(), self.entries.end())
        self.entries.erase(last, self.entries.end())
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
        self.entries.clear()

# Set version
cdef class SparseUnorderedRefinedBitmaskSet:
    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2):
        cdef ind_pair ind
        ind.first = ind1
        ind.second = ind2
        self.entries.insert(ind)

    def set(self, ind1, ind2):
        self._set(ind1, ind2)

    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:] mask2):
        for p in self.entries:
            mask1[p.first] = mask2[p.second] = 1

    cdef void _fill_ewah(self, BoolArrayCollection mm):
        for it in self.entries:
            mm._set_refined(it.first, it.second)

    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm):
        for it in self.entries:
            mm._set_refined(it.first, it.second)

    cdef void _reset(self):
        self.entries.clear()

    cdef to_array(self):
        cdef np.uint64_t i
        cdef np.ndarray[np.uint64_t, ndim=2] rv
        rv = np.empty((self.entries.size(),2),dtype='uint64')
        i = 0
        for it in self.entries:
            rv[i,0] = it.first
            rv[i,1] = it.second
            i += 1
        return rv

    def __dealloc__(self):
        self.entries.clear()
