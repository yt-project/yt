cimport numpy as np
from libcpp.pair cimport pair
from libcpp.set cimport set as cset
from libcpp.vector cimport vector

from ewah_bool_utils.ewah_bool_array cimport (
    ewah_bool_array,
    ewah_bool_iterator,
    ewah_map,
    sstream,
)

ctypedef bint bitarrtype
ctypedef pair[np.uint64_t, np.uint64_t] ind_pair

cdef class FileBitmasks:
    cdef np.uint32_t nfiles
    cdef ewah_map** ewah_coll
    cdef ewah_bool_array** ewah_keys
    cdef ewah_bool_array** ewah_refn

    cdef void _reset(self)
    cdef bint _iseq(self, FileBitmasks solf)
    cdef BoolArrayCollection _get_bitmask(self, np.uint32_t ifile)
    cdef tuple _find_collisions(self, BoolArrayCollection coll, bint verbose=*)
    cdef tuple _find_collisions_coarse(self, BoolArrayCollection coll, bint
                verbose=*, file_list=*)
    cdef tuple _find_collisions_refined(self, BoolArrayCollection coll, bint verbose=*)
    cdef void _set(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2=*)
    cdef void _set_coarse(self, np.uint32_t ifile, np.uint64_t i1)
    cdef void _set_refined(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_coarse_array(self, np.uint32_t ifile, np.uint8_t[:] arr)
    cdef void _set_refined_array(self, np.uint32_t ifile, np.uint64_t mi1, np.uint8_t[:] arr)
    cdef void _set_refined_index_array(self, np.uint32_t ifile, np.int64_t nsub_mi,
                                       np.ndarray[np.uint64_t, ndim=1] sub_mi1,
                                       np.ndarray[np.uint64_t, ndim=1] sub_mi2)
    cdef void _set_map(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_refn(self, np.uint32_t ifile, np.uint64_t i1)
    cdef bint _get(self, np.uint32_t ifile, np.uint64_t i1, np.uint64_t i2=*)
    cdef bint _get_coarse(self, np.uint32_t ifile, np.uint64_t i1)
    cdef void _get_coarse_array(self, np.uint32_t ifile, np.uint64_t imax, np.uint8_t[:] arr) except *
    cdef bint _isref(self, np.uint32_t ifile, np.uint64_t i)
    cdef np.uint64_t _count_total(self, np.uint32_t ifile)
    cdef np.uint64_t _count_refined(self, np.uint32_t ifile)
    cdef np.uint64_t _count_coarse(self, np.uint32_t ifile)
    cdef void _append(self, np.uint32_t ifile, BoolArrayCollection solf)
    cdef bint _intersects(self, np.uint32_t ifile, BoolArrayCollection solf)
    cdef void _logicalxor(self, np.uint32_t ifile, BoolArrayCollection solf, BoolArrayCollection out)
    cdef void _logicaland(self, np.uint32_t ifile, BoolArrayCollection solf, BoolArrayCollection out)
    cdef void _select_contaminated(self, np.uint32_t ifile, BoolArrayCollection mask, np.uint8_t[:] out,
               np.uint8_t[:] secondary_files, BoolArrayCollection mask2=*)
    cdef void _select_uncontaminated(self, np.uint32_t ifile, BoolArrayCollection mask, np.uint8_t[:] out,
               BoolArrayCollection mask2=*)
    cdef bytes _dumps(self, np.uint32_t ifile)
    cdef bint _loads(self, np.uint32_t ifile, bytes s)
    cdef bint _check(self)

cdef class BoolArrayCollection:
    cdef ewah_map* ewah_coll
    cdef ewah_bool_array* ewah_keys
    cdef ewah_bool_array* ewah_refn
    cdef ewah_bool_array* ewah_coar

    cdef void _reset(self)
    cdef int _richcmp(self, BoolArrayCollection solf, int op) except -1
    cdef void _set(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef void _set_coarse(self, np.uint64_t i1)
    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_coarse_array(self, np.uint8_t[:] arr)
    cdef void _set_refined_array(self, np.uint64_t mi1, np.uint8_t[:] arr)
    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_refn(self, np.uint64_t i1)
    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef bint _get_coarse(self, np.uint64_t i1)
    cdef void _get_coarse_array(self, np.uint64_t imax, np.uint8_t[:] arr) except *
    cdef bint _contains(self, np.uint64_t i)
    cdef bint _isref(self, np.uint64_t i)
    cdef void _ewah_coarse(self)
    cdef np.uint64_t _count_total(self)
    cdef np.uint64_t _count_refined(self)
    cdef np.uint64_t _count_coarse(self)
    cdef void _append(self, BoolArrayCollection solf)
    cdef void _logicalor(self, BoolArrayCollection solf, BoolArrayCollection out)
    cdef bint _intersects(self, BoolArrayCollection solf)
    cdef void _logicalxor(self, BoolArrayCollection solf, BoolArrayCollection out)
    cdef void _logicaland(self, BoolArrayCollection solf, BoolArrayCollection out)
    cdef void _select_contaminated(self, BoolArrayCollection mask, np.uint8_t[:] out,
        BoolArrayCollection mask2=*)
    cdef void _select_uncontaminated(self, BoolArrayCollection mask, np.uint8_t[:] out,
        BoolArrayCollection mask2=*)
    cdef void _get_ghost_zones(self, int ngz, int order1, int order2,
                               bint periodicity[3], BoolArrayCollection out_ewah,
                               bint coarse_ghosts=*)
    cdef bytes _dumps(self)
    cdef bint _loads(self, bytes s)
    cdef bint _check(self)

cdef class BoolArrayCollectionUncompressed:
    cdef int nele1
    cdef int nele2
    cdef ewah_map* ewah_coll
    cdef bitarrtype* ewah_keys
    cdef bitarrtype* ewah_refn

    cdef void _set(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef void _set_coarse(self, np.uint64_t i1)
    cdef void _set_refined(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_coarse_array(self, np.uint8_t[:] arr)
    cdef void _set_coarse_array_ptr(self, np.uint8_t *arr)
    cdef void _set_refined_array(self, np.uint64_t mi1, np.uint8_t[:] arr)
    cdef void _set_refined_array_ptr(self, np.uint64_t mi1, np.uint8_t *arr)
    cdef void _set_map(self, np.uint64_t i1, np.uint64_t i2)
    cdef void _set_refn(self, np.uint64_t i1)
    cdef bint _get(self, np.uint64_t i1, np.uint64_t i2=*)
    cdef bint _get_coarse(self, np.uint64_t i1)
    cdef bint _isref(self, np.uint64_t i)
    cdef np.uint64_t _count_total(self)
    cdef np.uint64_t _count_refined(self)
    cdef void _append(self, BoolArrayCollectionUncompressed solf)
    cdef bint _intersects(self, BoolArrayCollectionUncompressed solf)
    cdef void _compress(self, BoolArrayCollection solf)

cdef class SparseUnorderedBitmaskSet:
    cdef cset[np.uint64_t] entries
    cdef void _set(self, np.uint64_t ind)
    cdef void _fill(self, np.uint8_t[:] mask)
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm)
    cdef void _reset(self)
    cdef to_array(self)

cdef class SparseUnorderedBitmaskVector:
    cdef int total
    cdef vector[np.uint64_t] entries
    cdef void _set(self, np.uint64_t ind)
    cdef void _fill(self, np.uint8_t[:] mask)
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm)
    cdef void _reset(self)
    cdef to_array(self)
    cdef void _remove_duplicates(self)
    cdef void _prune(self)

cdef class SparseUnorderedRefinedBitmaskSet:
    cdef cset[ind_pair] entries
    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2)
    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:])
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm)
    cdef void _reset(self)
    cdef to_array(self)

cdef class SparseUnorderedRefinedBitmaskVector:
    cdef int total
    cdef vector[ind_pair] entries
    cdef void _set(self, np.uint64_t ind1, np.uint64_t ind2)
    cdef void _fill(self, np.uint8_t[:] mask1, np.uint8_t[:])
    cdef void _fill_ewah(self, BoolArrayCollection mm)
    cdef void _fill_bool(self, BoolArrayCollectionUncompressed mm)
    cdef void _reset(self)
    cdef to_array(self)
    cdef void _remove_duplicates(self)
    cdef void _prune(self)
