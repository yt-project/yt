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

cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
cimport cython
from libc.math cimport sqrt

from fp_utils cimport *
from oct_container cimport Oct, OctAllocationContainer, OctreeContainer
from .particle_deposit cimport sph_kernel, gind

cdef extern from "alloca.h":
    void *alloca(int)

cdef class ParticleSmoothOperation:
    # We assume each will allocate and define their own temporary storage
    cdef public object nvals
    cdef void process(self, int dim[3], np.float64_t left_edge[3],
                      np.float64_t dds[3], np.float64_t *ppos,
                      np.float64_t **fields, np.int64_t nneighbors,
                      np.int64_t *nind, np.int64_t *doffs,
                      np.int64_t *pinds, np.int64_t *pcounts,
                      np.int64_t offset)
