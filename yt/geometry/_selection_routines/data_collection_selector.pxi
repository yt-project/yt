cdef class DataCollectionSelector(SelectorObject):
    cdef object obj_ids
    cdef np.int64_t nids

    def __init__(self, dobj):
        self.obj_ids = dobj._obj_ids
        self.nids = self.obj_ids.shape[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.ndarray[np.int64_t, ndim=1] oids = self.obj_ids
        with nogil:
            for n in range(self.nids):
                gridi[oids[n]] = 1
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        mask = np.ones(gobj.ActiveDimensions, dtype='uint8')
        return mask.astype("bool")

    def _hash_vals(self):
        return (hash(self.obj_ids.tobytes()), self.nids)

data_collection_selector = DataCollectionSelector
