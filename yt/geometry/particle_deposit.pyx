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
        cdef np.int64_t offset, moff
        cdef Oct *oct
        moff = octree.get_domain_offset(domain_id)
        for i in range(positions.shape[0]):
            # We should check if particle remains inside the Oct here
            for j in range(nf):
                field_vals[j] = field_pointers[j][i]
            for j in range(3):
                pos[j] = positions[i, j]
            oct = octree.get(pos, &oi)
            # This next line is unfortunate.  Basically it says, sometimes we
            # might have particles that belong to octs outside our domain.
            if oct.domain != domain_id: continue
            #print domain_id, oct.local_ind, oct.ind, oct.domain, oct.pos[0], oct.pos[1], oct.pos[2]
            # Note that this has to be our local index, not our in-file index.
            offset = dom_ind[oct.domain_ind - moff] * 8
            if offset < 0: continue
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
    cdef np.int64_t *count # float, for ease
    cdef public object ocount
    def initialize(self):
        # Create a numpy array accessible to python
        self.ocount = np.zeros(self.nvals, dtype="int64")
        cdef np.ndarray arr = self.ocount
        # alias the C-view for use in cython
        self.count = <np.int64_t*> arr.data

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
        return self.ocount.astype('f8')

deposit_count = CountParticles

cdef class SumParticleField(ParticleDepositOperation):
    cdef np.float64_t *sum
    cdef public object osum
    def initialize(self):
        self.osum = np.zeros(self.nvals, dtype="float64")
        cdef np.ndarray arr = self.osum
        self.sum = <np.float64_t*> arr.data

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3], 
                      np.float64_t dds[3],
                      np.int64_t offset, 
                      np.float64_t ppos[3],
                      np.float64_t *fields 
                      ):
        cdef int ii[3], i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i]) / dds[i])
        self.sum[gind(ii[0], ii[1], ii[2], dim) + offset] += fields[0]
        
    def finalize(self):
        return self.osum

deposit_sum = SumParticleField

cdef class StdParticleField(ParticleDepositOperation):
    # Thanks to Britton and MJ Turk for the link
    # to a single-pass STD
    # http://www.cs.berkeley.edu/~mhoemmen/cs194/Tutorials/variance.pdf
    cdef np.float64_t *mk
    cdef np.float64_t *qk
    cdef np.float64_t *i
    cdef public object omk
    cdef public object oqk
    cdef public object oi
    def initialize(self):
        # we do this in a single pass, but need two scalar
        # per cell, M_k, and Q_k and also the number of particles
        # deposited into each one
        # the M_k term
        self.omk= np.zeros(self.nvals, dtype="float64")
        cdef np.ndarray omkarr= self.omk
        self.mk= <np.float64_t*> omkarr.data
        # the Q_k term
        self.oqk= np.zeros(self.nvals, dtype="float64")
        cdef np.ndarray oqkarr= self.oqk
        self.qk= <np.float64_t*> oqkarr.data
        # particle count
        self.oi = np.zeros(self.nvals, dtype="float64")
        cdef np.ndarray oiarr = self.oi
        self.i = <np.float64_t*> oiarr.data

    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3], 
                      np.float64_t dds[3],
                      np.int64_t offset,
                      np.float64_t ppos[3],
                      np.float64_t *fields
                      ):
        cdef int ii[3], i, cell_index
        cdef float k, mk, qk
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i])/dds[i])
        cell_index = gind(ii[0], ii[1], ii[2], dim) + offset
        k = self.i[cell_index] 
        mk = self.mk[cell_index]
        qk = self.qk[cell_index] 
        #print k, mk, qk, cell_index
        if k == 0.0:
            # Initialize cell values
            self.mk[cell_index] = fields[0]
        else:
            self.mk[cell_index] = mk + (fields[0] - mk) / k
            self.qk[cell_index] = qk + (k - 1.0) * (fields[0] - mk)**2.0 / k
        self.i[cell_index] += 1
        
    def finalize(self):
        # This is the standard variance
        # if we want sample variance divide by (self.oi - 1.0)
        std2 = self.oqk / self.oi
        std2[self.oi == 0.0] = 0.0
        return np.sqrt(std2)

deposit_std = StdParticleField


cdef class WeightedMeanParticleField(ParticleDepositOperation):
    # Deposit both mass * field and mass into two scalars
    # then in finalize divide mass * field / mass
    cdef np.float64_t *wf
    cdef public object owf
    cdef np.float64_t *w
    cdef public object ow
    def initialize(self):
        self.owf = np.zeros(self.nvals, dtype='float64')
        cdef np.ndarray wfarr = self.owf
        self.wf = <np.float64_t*> wfarr.data
        
        self.ow = np.zeros(self.nvals, dtype='float64')
        cdef np.ndarray warr = self.ow
        self.w = <np.float64_t*> warr.data
    
    @cython.cdivision(True)
    cdef void process(self, int dim[3],
                      np.float64_t left_edge[3], 
                      np.float64_t dds[3],
                      np.int64_t offset, 
                      np.float64_t ppos[3],
                      np.float64_t *fields 
                      ):
        cdef int ii[3], i
        for i in range(3):
            ii[i] = <int>((ppos[i] - left_edge[i]) / dds[i])
        self.w[ gind(ii[0], ii[1], ii[2], dim) + offset] += fields[1]
        self.wf[gind(ii[0], ii[1], ii[2], dim) + offset] += fields[0] * fields[1]
        
    def finalize(self):
        return self.owf / self.ow

deposit_weighted_mean= WeightedMeanParticleField

