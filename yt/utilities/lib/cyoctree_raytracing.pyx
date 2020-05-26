"""This is a wrapper around the C++ class to efficiently cast rays into an octree.
It relies on the seminal paper by  J. Revelles,, C.Ureña and M.Lastra.
"""


cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
cimport cython
from libc.stdlib cimport free

cdef extern from "octree_raytracing.cpp":
    cdef cppclass RayInfo[T]:
        vector[T] keys
        vector[double] t

    cdef cppclass Octree3D[T]:
        Octree3D(int depth, double* size)
        Octree3D(int depth, double* LE, double* RE)
        void insert_node_no_ret(const int* ipos, const int lvl, T key)
        RayInfo[T]** cast_rays(const double* origins, const double* directions, const int Nrays)
        
cdef class CythonOctreeRayTracing:
    cdef Octree3D[int]* oct
    cdef int depth

    def __init__(self, np.ndarray LE, np.ndarray RE, int depth):
        cdef double* LE_ptr = <double *>LE.data
        cdef double* RE_ptr = <double *>RE.data
        self.oct = new Octree3D[int](depth, LE_ptr, RE_ptr)
        self.depth = depth
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_nodes(self, int[:, :] ipos_view, int[:] lvl_view, int[:] key):
        cdef int i
        cdef int ii[3]
        
        for i in range(len(key)):
            ii[0] = ipos_view[i, 0]
            ii[1] = ipos_view[i, 1]
            ii[2] = ipos_view[i, 2]
            self.oct.insert_node_no_ret(ii, lvl_view[i], <int> key[i])

    @cython.boundscheck(False)
    @cython.wraparound(False)       
    def cast_rays(self, double[:, ::1] o, double[:, ::1] d):
        cdef RayInfo[int]** ret
        cdef int Nrays = len(o)
        cdef RayInfo[int]* ri
        
        if Nrays == 0:
            return
        
        # print('Casting rays')
        
        ret = self.oct.cast_rays(&o[0,0], &d[0,0], Nrays)

        # print('cast!')
        # Now pack all the rays in numpy arrays
        cdef int[:] key_view
        cdef int* key_ptr
        
        cdef double[:] t_view
        cdef double* t_ptr
        
        # print('Taking ownership of data')
        
        key_array, t_array = [], []
        for i in range(Nrays):
            ri = ret[i]
            if ri.keys.size() == 0:
                key_array.append(np.array([], dtype=int))
                t_array.append(np.array([], dtype=np.float64))
            else:
                key_ptr = &ri.keys[0]
                key_view = <int[:ri.keys.size()]> key_ptr
                key_array.append(np.asarray(key_view))
                
                t_ptr = &ri.t[0]
                t_view = <double[:ri.t.size()]> t_ptr
                t_array.append(np.asarray(t_view).reshape(-1, 2))
                
            free(ret[i])
            
        # We can now free the *list* of vectors, note that the memory is now managed by numpy!
        free(ret)
        return key_array, t_array

    def __dealloc__(self):
        del self.oct