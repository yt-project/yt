"""
Particle smoothing in cells

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
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
from libc.stdlib cimport malloc, free, realloc
cimport cython
from libc.math cimport sqrt

from fp_utils cimport *
from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, OctInfo

cdef class ParticleSmoothOperation:
    def __init__(self, nvals, nfields, max_neighbors):
        # This is the set of cells, in grids, blocks or octs, we are handling.
        cdef int i
        self.nvals = nvals 
        self.nfields = nfields
        self.maxn = max_neighbors
        self.neighbors = <NeighborList *> malloc(
            sizeof(NeighborList) * self.maxn)
        self.neighbor_reset()

    def initialize(self, *args):
        raise NotImplementedError

    def finalize(self, *args):
        raise NotImplementedError

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_octree(self, OctreeContainer octree,
                     np.ndarray[np.int64_t, ndim=1] dom_ind,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None, int domain_id = -1,
                     int domain_offset = 0,
                     int test_neighbors = 0):
        # This will be a several-step operation.
        #
        # We first take all of our particles and assign them to Octs.  If they
        # are not in an Oct, we will assume they are out of bounds.  Note that
        # this means that if we have loaded neighbor particles for which an Oct
        # does not exist, we are going to be discarding them -- so sparse
        # octrees will need to ensure that neighbor octs *exist*.  Particles
        # will be assigned in a new NumPy array.  Note that this incurs
        # overhead, but reduces complexity as we will now be able to use
        # argsort.
        #
        # After the particles have been assigned to Octs, we process each Oct
        # individually.  We will do this by calling "get" for the *first*
        # particle in each set of Octs in the sorted list.  After this, we get
        # neighbors for each Oct.
        #
        # Now, with the set of neighbors (and thus their indices) we allocate
        # an array of particles and their fields, fill these in, and call our
        # process function.
        #
        # This is not terribly efficient -- for starters, the neighbor function
        # is not the most efficient yet.  We will also need to handle some
        # mechanism of an expandable array for holding pointers to Octs, so
        # that we can deal with >27 neighbors.  As I write this comment,
        # neighbors() only returns 27 neighbors.
        cdef int nf, i, j, dims[3]
        cdef np.float64_t **field_pointers, *field_vals, pos[3], *ppos
        cdef int nsize = 0
        cdef np.int64_t *nind = NULL
        cdef OctInfo oi
        cdef Oct *oct, **neighbors = NULL
        cdef np.int64_t nneighbors, numpart, offset, moff, local_ind
        cdef np.int64_t *doffs, *pinds, *pcounts
        cdef np.ndarray[np.int64_t, ndim=1] pind, doff, pdoms, pcount
        cdef np.ndarray[np.float64_t, ndim=1] tarr
        dims[0] = dims[1] = dims[2] = 2
        numpart = positions.shape[0]
        pcount = np.zeros_like(dom_ind)
        doff = np.zeros_like(dom_ind) - 1
        moff = octree.get_domain_offset(domain_id + domain_offset)
        pdoms = np.zeros(positions.shape[0], dtype="int64") - 1
        nf = len(fields)
        if fields is None:
            fields = []
        field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        for i in range(nf):
            tarr = fields[i]
            field_pointers[i] = <np.float64_t *> tarr.data
        for i in range(3):
            self.DW[i] = (octree.DRE[i] - octree.DLE[i])
        for i in range(positions.shape[0]):
            for j in range(3):
                pos[j] = positions[i, j]
            oct = octree.get(pos)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            # Note that this has to be our local index, not our in-file index.
            # This is the particle count, which we'll use once we have sorted
            # the particles to calculate the offsets into each oct's particles.
            pcount[oct.domain_ind - moff] += 1
            pdoms[i] = oct.domain_ind - moff # We store the *actual* offset.
        # Now we have oct assignments.  Let's sort them.
        # Note that what we will be providing to our processing functions will
        # actually be indirectly-sorted fields.  This preserves memory at the
        # expense of additional pointer lookups.
        pind = np.argsort(pdoms)
        # This now gives us the indices to the particles for each domain.
        for i in range(positions.shape[0]):
            # This is the domain_ind (minus moff) for this particle
            offset = pdoms[pind[i]] 
            if doff[offset] < 0: doff[offset] = i
        # Now doff is full of offsets to the first entry in the pind that
        # refers to that oct's particles.
        ppos = <np.float64_t *> positions.data
        doffs = <np.int64_t*> doff.data
        pinds = <np.int64_t*> pind.data
        pcounts = <np.int64_t*> pcount.data
        cdef np.int64_t pn
        nsize = 27
        nind = <np.int64_t *> malloc(sizeof(np.int64_t)*nsize)
        for i in range(doff.shape[0]):
            # Nothing assigned.
            if doff[i] < 0: continue
            for j in range(3):
                pos[j] = positions[pind[doff[i]], j]
            oct = octree.get(pos, &oi)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            offset = dom_ind[oct.domain_ind - moff] * 8
            neighbors = octree.neighbors(&oi, &nneighbors)
            # Now we have all our neighbors.  And, we should be set for what
            # else we need to do.
            if nneighbors > nsize:
                nind = <np.int64_t *> realloc(
                    nind, sizeof(np.int64_t)*nneighbors)
                nsize = nneighbors
            for j in range(nneighbors):
                nind[j] = neighbors[j].domain_ind - moff
            free(neighbors)
            self.neighbor_process(dims, oi.left_edge, oi.dds,
                         ppos, field_pointers, nneighbors, nind, doffs,
                         pinds, pcounts, offset)
            # This is allocated by the neighbors function, so we deallocate it.
        if nind != NULL:
            free(nind)
        
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_grid(self, gobj,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None):
        raise NotImplementedError

    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3]):
        raise NotImplementedError

    cdef void neighbor_reset(self):
        self.curn = 0
        for i in range(self.maxn):
            self.neighbors[i].pn = -1
            self.neighbors[i].r2 = 1e300

    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3]):
        cdef NeighborList *cur
        cdef int i
        # _c means candidate (what we're evaluating)
        # _o means other (the item in the list)
        cdef np.float64_t r2_c, r2_o
        cdef np.int64_t pn_c, pn_o
        if self.curn < self.maxn:
            cur = &self.neighbors[self.curn]
            cur.pn = pn
            cur.r2 = r2dist(ppos, cpos, self.DW)
            self.curn += 1
            return
        # This will go (curn - 1) through 0.
        r2_c = r2dist(ppos, cpos, self.DW)
        pn_c = pn
        for i in range((self.curn - 1), -1, -1):
            # First we evaluate against i.  If our candidate radius is greater
            # than the one we're inspecting, we quit early.
            cur = &self.neighbors[i]
            r2_o = cur.r2
            pn_o = cur.pn
            if r2_c >= r2_o:
                break
            # Now we know we need to swap them.  First we assign our candidate
            # values to cur.
            cur.r2 = r2_c
            cur.pn = pn_c
            if i + 1 >= self.maxn:
                continue # No swapping
            cur = &self.neighbors[i + 1]
            cur.r2 = r2_o
            cur.pn = pn_o
        # At this point, we've evaluated all the particles and we should have a
        # sorted set of values.  So, we're done.

    cdef void neighbor_find(self,
                            np.int64_t nneighbors,
                            np.int64_t *nind,
                            np.int64_t *doffs,
                            np.int64_t *pcounts,
                            np.int64_t *pinds,
                            np.float64_t *ppos,
                            np.float64_t cpos[3]
                            ):
        # We are now given the number of neighbors, the indices into the
        # domains for them, and the number of particles for each.
        cdef int ni, i, j
        cdef np.int64_t offset, pn, pc
        cdef np.float64_t pos[3]
        self.neighbor_reset()
        for ni in range(nneighbors):
            offset = doffs[nind[ni]]
            pc = pcounts[nind[ni]]
            for i in range(pc):
                pn = pinds[offset + i]
                for j in range(3):
                    pos[j] = ppos[pn * 3 + j]
                self.neighbor_eval(pn, pos, cpos)

    cdef void neighbor_process(self, int dim[3], np.float64_t left_edge[3],
                               np.float64_t dds[3], np.float64_t *ppos,
                               np.float64_t **fields, np.int64_t nneighbors,
                               np.int64_t *nind, np.int64_t *doffs,
                               np.int64_t *pinds, np.int64_t *pcounts,
                               np.int64_t offset):
        # Note that we assume that fields[0] == smoothing length in the native
        # units supplied.  We can now iterate over every cell in the block and
        # every particle to find the nearest.  We will use a priority heap.
        cdef int i, j, k
        cdef np.float64_t cpos[3]
        cpos[0] = left_edge[0] + 0.5*dds[0]
        for i in range(dim[0]):
            cpos[1] = left_edge[1] + 0.5*dds[1]
            for j in range(dim[1]):
                cpos[2] = left_edge[2] + 0.5*dds[2]
                for k in range(dim[2]):
                    self.neighbor_find(nneighbors, nind, doffs, pcounts,
                        pinds, ppos, cpos)
                    # Now we have all our neighbors in our neighbor list.
                    self.process(offset, i, j, k, dim, cpos)
                    cpos[2] += dds[2]
                cpos[1] += dds[1]
            cpos[0] += dds[0]


cdef class SimpleNeighborSmooth(ParticleSmoothOperation):
    cdef np.float64_t **fp
    cdef public object vals
    def initialize(self):
        if self.nvals < 2:
            # We need at least two fields, the smoothing length and the 
            # field to smooth, to operate.
            raise RuntimeError
        self.vals = [np.zeros(self.nvals, dtype="float64")
                     for i in range(self.nfields)]
        self.fp = <np.float64_t **> malloc(
            sizeof(np.float64_t *) * self.nfields)

    def finalize(self):
        free(self.fp)
        return self.vals

    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3]):
        # We have our i, j, k for our cell 
        #print "Offset", offset, i, j, k, self.curn
        return

simple_neighbor_smooth = SimpleNeighborSmooth
