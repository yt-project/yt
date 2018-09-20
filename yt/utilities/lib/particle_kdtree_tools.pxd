cimport numpy as np
cimport cython

from cykdtree.kdtree cimport PyKDTree, KDTree, Node, uint64_t, uint32_t
from yt.utilities.lib.bounded_priority_queue cimport BoundedPriorityQueue

cdef int knn_position(np.float64_t[:, ::1] tree_positions,
                      np.float64_t[::1] position,
                      BoundedPriorityQueue queue, KDTree * kdtree,
                      np.int64_t skipaxis)

cdef int knn_grid(np.float64_t[:, ::1] tree_positions,
                  np.float64_t[:, :, :, ::1] dists,
                  np.int64_t[:, :, :, ::1] pids,  KDTree * kdtree,
                  np.float64_t[:] bounds, np.int64_t[:] size,
                  np.int64_t skipaxis)
