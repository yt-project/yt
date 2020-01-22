"""
RAMSES Oct definitions file




"""
from oct_container cimport SparseOctreeContainer, OctInfo
from .oct_visitors cimport OctVisitor, Oct, cind
from yt.utilities.lib.fp_utils cimport *
cimport numpy as np

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    cdef Oct neighbour_in_direction(self, OctInfo *oinfo, np.int64_t *nneighbors,
                                    Oct *o, bint periodicity[3])
