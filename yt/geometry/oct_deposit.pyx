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

from libc.stdlib cimport malloc, free, qsort
from libc.math cimport floor
cimport numpy as np
import numpy as np
from oct_container cimport Oct, OctAllocationContainer, OctreeContainer
cimport cython

cdef np.float64_t kernel_sph(np.float64_t x) nogil:
    cdef np.float64_t kernel
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel

#modes = count, sum, diff
modes = {'count': opt_count, 'sum': opt_sum, 'diff': opt_diff}
selections = {'direct': select_nearest, 'cic': select_radius}
kernels = {'unitary': kernel_unitary, 'sph': kernel_sph}
cdef deposit_direct(oct_handler, 
        np.ndarray[np.float64_t, ndim=2] ppos, #positions,columns are x,y,z
        np.ndarray[np.float64_t, ndim=2] pd, # particle fields
        np.ndarray[np.float64_t, ndim=1] pr, # particle radius
        np.ndarray[np.float64_t, ndim=2] data_out, #write deposited here
        np.ndarray[np.float64_t, ndim=2] data_in, #used to calc diff, same shape as data_out
        mode='count', selection='direct', kernel='sph'):
    fopt = modes[mode]
    fsel = selections[selection]
    fker = kernels[kernel]
    for pi in np.arange(particles):
        octs = fsel(oct_handler, pr[pi])
        for oct in octs:
            w = fker(pr[pi],oct) 
            weights.append(w)
        norm = weights.sum()
        for w, oct in zip(weights, octs):
            fopt(pd[pi], w/norm, oct.index, data_in, data_out)


