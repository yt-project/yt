"""
Particle Deposition onto Octs

Author: Christopher Moody <chris.e.moody@gmail.com>
Affiliation: UC Santa Cruz
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np
cimport cython

from oct_container cimport Oct, OctAllocationContainer, OctreeContainer

# Mode functions
ctypedef np.float64_t (*type_opt)(np.float64_t, np.float64_t)
cdef np.float64_t opt_count(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += 1.0

cdef np.float64_t opt_sum(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += pdata 

cdef np.float64_t opt_diff(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += (data_in[index] - pdata) 

cdef np.float64_t opt_wcount(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += weight

cdef np.float64_t opt_wsum(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += pdata * weight

cdef np.float64_t opt_wdiff(np.float64_t pdata,
                            np.float64_t weight,
                            np.int64_t index,
                            np.ndarray[np.float64_t, ndim=2] data_out, 
                            np.ndarray[np.float64_t, ndim=2] data_in):
    data_out[index] += (data_in[index] - pdata) * weight

# Selection functions
ctypedef NOTSURE (*type_sel)(OctreeContainer, 
                                np.ndarray[np.float64_t, ndim=1],
                                np.float64_t)
cdef NOTSURE select_nearest(OctreeContainer oct_handler,
                            np.ndarray[np.float64_t, ndim=1] pos,
                            np.float64_t radius):
    #return only the nearest oct
    pass


cdef NOTSURE select_radius(OctreeContainer oct_handler,
                            np.ndarray[np.float64_t, ndim=1] pos,
                            np.float64_t radius):
    #return a list of octs within the radius
    pass
    

# Kernel functions
ctypedef np.float64_t (*type_ker)(np.float64_t)
cdef np.float64_t kernel_sph(np.float64_t x) nogil:
    cdef np.float64_t kernel
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel

cdef np.float64_t kernel_null(np.float64_t x) nogil: return 0.0

cdef deposit(OctreeContainer oct_handler, 
        np.ndarray[np.float64_t, ndim=2] ppos, #positions,columns are x,y,z
        np.ndarray[np.float64_t, ndim=2] pd, # particle fields
        np.ndarray[np.float64_t, ndim=1] pr, # particle radius
        np.ndarray[np.float64_t, ndim=2] data_in, #used to calc diff, same shape as data_out
        np.ndarray[np.float64_t, ndim=2] data_out, #write deposited here
        mode='count', selection='nearest', kernel='null'):
    cdef type_opt fopt
    cdef type_sel fsel
    cdef type_ker fker
    cdef long pi #particle index
    cdef long nocts #number of octs in selection
    cdef Oct oct 
    cdef np.float64_t w
    # Can we do this with dicts?
    # Setup the function pointers
    if mode == 'count':
        fopt = opt_count
    elif mode == 'sum':
        fopt = opt_sum
    elif mode == 'diff':
        fopt = opt_diff
    if mode == 'wcount':
        fopt = opt_count
    elif mode == 'wsum':
        fopt = opt_sum
    elif mode == 'wdiff':
        fopt = opt_diff
    if selection == 'nearest':
        fsel = select_nearest
    elif selection == 'radius':
        fsel = select_radius
    if kernel == 'null':
        fker = kernel_null
    if kernel == 'sph':
        fker = kernel_sph
    for pi in range(particles):
        octs = fsel(oct_handler, ppos[pi], pr[pi])
        for oct in octs:
            for cell in oct.cells:
                w = fker(pr[pi],cell) 
                weights.append(w)
        norm = weights.sum()
        for w, oct in zip(weights, octs):
            for cell in oct.cells:
                fopt(pd[pi], w/norm, oct.index, data_in, data_out)


