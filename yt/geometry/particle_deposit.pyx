"""
Particle Deposition onto Cells

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

from fp_utils cimport *
from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, OctInfo

cdef class ParticleDepositOperation:
    def __init__(self, nvals):
        self.nvals = nvals

    def initialize(self, *args):
        raise NotImplementedError

    def finalize(self, *args):
        raise NotImplementedError

    def process_octree(self, OctreeContainer octree,
                     np.ndarray[np.int64_t, ndim=1] dom_ind,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None, int domain_id = -1):
        cdef int nf, i, j
        if fields is None:
            fields = []
        nf = len(fields)
        cdef np.float64_t **field_pointers, *field_vals, pos[3]
        cdef np.ndarray[np.float64_t, ndim=1] tarr
        field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        field_vals = <np.float64_t*>alloca(sizeof(np.float64_t) * nf)
        for i in range(nf):
            tarr = fields[i]
            field_pointers[i] = <np.float64_t *> tarr.data
        cdef int dims[3]
        dims[0] = dims[1] = dims[2] = 2
        cdef OctInfo oi
        cdef np.int64_t offset
        cdef Oct *oct
        for i in range(positions.shape[0]):
            # We should check if particle remains inside the Oct here
            for j in range(nf):
                field_vals[j] = field_pointers[j][i]
            for j in range(3):
                pos[j] = positions[i, j]
            oct = octree.get(pos, &oi)
            #print oct.local_ind, oct.pos[0], oct.pos[1], oct.pos[2]
            offset = dom_ind[oct.ind] * 8
            # Check that we found the oct ...
            self.process(dims, oi.left_edge, oi.dds,
                         offset, pos, field_vals)
        
    def process_grid(self, gobj,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None):
        cdef int nf, i, j
        if fields is None:
            fields = []
        nf = len(fields)
        cdef np.float64_t **field_pointers, *field_vals, pos[3]
        cdef np.ndarray[np.float64_t, ndim=1] tarr
        field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        field_vals = <np.float64_t*>alloca(sizeof(np.float64_t) * nf)
        for i in range(nf):
            tarr = fields[i]
            field_pointers[i] = <np.float64_t *> tarr.data
        cdef np.float64_t dds[3], left_edge[3]
        cdef int dims[3]
        for i in range(3):
            dds[i] = gobj.dds[i]
            left_edge[i] = gobj.LeftEdge[i]
            dims[i] = gobj.ActiveDimensions[i]
        for i in range(positions.shape[0]):
            # Now we process
            for j in range(nf):
                field_vals[j] = field_pointers[j][i]
            for j in range(3):
                pos[j] = positions[i, j]
            self.process(dims, left_edge, dds, 0, pos, field_vals)

    cdef void process(self, int dim[3], np.float64_t left_edge[3],
                      np.float64_t dds[3], np.int64_t offset,
                      np.float64_t ppos[3], np.float64_t *fields):
        raise NotImplementedError

cdef class CountParticles(ParticleDepositOperation):
    cdef np.float64_t *count # float, for ease
    cdef public object ocount
    def initialize(self):
        self.ocount = np.zeros(self.nvals, dtype="float64")
        cdef np.ndarray arr = self.ocount
        self.count = <np.float64_t*> arr.data

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3], 
                      np.float64_t dds[3],
                      np.int64_t offset, # offset into IO field
                      np.float64_t ppos[3], # this particle's position
                      np.float64_t *fields # any other fields we need
                      ):
        # here we do our thing; this is the kernel
        cdef int ii[3], i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i])/dds[i])
        self.count[gind(ii[0], ii[1], ii[2], dim) + offset] += 1
        
    def finalize(self):
        return self.ocount

deposit_count = CountParticles

"""
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
"""
