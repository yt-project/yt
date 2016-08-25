import numpy as np
cimport numpy as np
cimport cython
from yt.utilities.physical_constants import kboltz
from libc.math cimport exp, fabs, sqrt

cdef double kb = kboltz.v
cdef double pi = np.pi

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def compute_weight(np.uint8_t thermal_broad,
    	           double dv,
                   double m_part,
	               np.ndarray[np.float64_t, ndim=1] v,
    		       np.ndarray[np.float64_t, ndim=1] T):

    cdef int i, n
    cdef double v2_th, x
    cdef np.ndarray[np.float64_t, ndim=1] w

    n = v.shape[0]
    w = np.zeros(n)

    for i in range(n):
        if thermal_broad:
            if T[i] > 0.0:
                v2_th = 2.*kb*T[i]/m_part
                w[i] = dv*exp(-v[i]*v[i]/v2_th)/sqrt(v2_th*pi)
        else:
            x = 1.-fabs(v[i])/dv
            if x > 0.0:
                w[i] = x
                
    return w
