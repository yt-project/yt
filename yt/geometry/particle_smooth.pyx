"""
Particle smoothing in cells




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free, realloc
cimport cython
from libc.math cimport sqrt, fabs, sin, cos

from fp_utils cimport *
from oct_container cimport Oct, OctAllocationContainer, \
    OctreeContainer, OctInfo

cdef int Neighbor_compare(void *on1, void *on2) nogil:
    cdef NeighborList *n1, *n2
    n1 = <NeighborList *> on1
    n2 = <NeighborList *> on2
    # Note that we set this up so that "greatest" evaluates to the *end* of the
    # list, so we can do standard radius comparisons.
    if n1.r2 < n2.r2:
        return -1
    elif n1.r2 == n2.r2:
        return 0
    else:
        return 1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.float64_t r2dist(np.float64_t ppos[3],
                         np.float64_t cpos[3],
                         np.float64_t DW[3],
                         bint periodicity[3],
                         np.float64_t max_dist2):
    cdef int i
    cdef np.float64_t r2, DR
    r2 = 0.0
    for i in range(3):
        DR = (ppos[i] - cpos[i])
        if not periodicity[i]:
            pass
        elif (DR > DW[i]/2.0):
            DR -= DW[i]
        elif (DR < -DW[i]/2.0):
            DR += DW[i]
        r2 += DR * DR
        if max_dist2 >= 0.0 and r2 > max_dist2:
            return -1.0
    return r2

cdef void spherical_coord_setup(np.float64_t ipos[3], np.float64_t opos[3]):
    opos[0] = ipos[0] * sin(ipos[1]) * cos(ipos[2])
    opos[1] = ipos[0] * sin(ipos[1]) * sin(ipos[2])
    opos[2] = ipos[0] * cos(ipos[1])

cdef void cart_coord_setup(np.float64_t ipos[3], np.float64_t opos[3]):
    opos[0] = ipos[0]
    opos[1] = ipos[1]
    opos[2] = ipos[2]

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
    def process_octree(self, OctreeContainer mesh_octree,
                     np.ndarray[np.int64_t, ndim=1] mdom_ind,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     np.ndarray[np.float64_t, ndim=2] oct_positions,
                     fields = None, int domain_id = -1,
                     int domain_offset = 0,
                     periodicity = (True, True, True),
                     index_fields = None,
                     OctreeContainer particle_octree = None,
                     np.ndarray[np.int64_t, ndim=1] pdom_ind = None,
                     geometry = "cartesian"):
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
        # that we can deal with >27 neighbors.  
        if particle_octree is None:
            particle_octree = mesh_octree
            pdom_ind = mdom_ind
        cdef int nf, i, j, dims[3], n
        cdef np.float64_t **field_pointers, *field_vals, pos[3], *ppos, dds[3]
        cdef np.float64_t **octree_field_pointers
        cdef int nsize = 0
        cdef np.int64_t *nind = NULL
        cdef OctInfo moi, poi
        cdef Oct *oct
        cdef np.int64_t numpart, offset, local_ind
        cdef np.int64_t moff_p, moff_m
        cdef np.int64_t *doffs, *pinds, *pcounts, poff
        cdef np.ndarray[np.int64_t, ndim=1] pind, doff, pdoms, pcount
        cdef np.ndarray[np.int64_t, ndim=2] doff_m
        cdef np.ndarray[np.float64_t, ndim=1] tarr
        cdef np.ndarray[np.float64_t, ndim=4] iarr
        cdef np.ndarray[np.float64_t, ndim=2] cart_positions
        if geometry == "cartesian":
            self.pos_setup = cart_coord_setup
            cart_positions = positions
        elif geometry == "spherical":
            self.pos_setup = spherical_coord_setup
            cart_positions = np.empty((positions.shape[0], 3), dtype="float64")

            cart_positions[:,0] = positions[:,0] * \
                                  np.sin(positions[:,1]) * \
                                  np.cos(positions[:,2])
            cart_positions[:,1] = positions[:,0] * \
                                  np.sin(positions[:,1]) * \
                                  np.sin(positions[:,2])
            cart_positions[:,2] = positions[:,0] * \
                                  np.cos(positions[:,1])
            periodicity = (False, False, False)
        else:
            raise NotImplementedError
        dims[0] = dims[1] = dims[2] = (1 << mesh_octree.oref)
        cdef int nz = dims[0] * dims[1] * dims[2]
        numpart = positions.shape[0]
        # pcount is the number of particles per oct.
        pcount = np.zeros_like(pdom_ind)
        # doff is the offset to a given oct in the sorted particles.
        doff = np.zeros_like(pdom_ind) - 1
        doff_m = np.zeros((mdom_ind.shape[0], 2), dtype="int64")
        moff_p = particle_octree.get_domain_offset(domain_id + domain_offset)
        moff_m = mesh_octree.get_domain_offset(domain_id + domain_offset)
        # pdoms points particles at their octs.  So the value in this array, for
        # a given index, is the local oct index.
        pdoms = np.zeros(positions.shape[0], dtype="int64") - 1
        nf = len(fields)
        if fields is None:
            fields = []
        field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        for i in range(nf):
            tarr = fields[i]
            field_pointers[i] = <np.float64_t *> tarr.data
        if index_fields is None:
            index_fields = []
        nf = len(index_fields)
        index_field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        for i in range(nf):
            iarr = index_fields[i]
            index_field_pointers[i] = <np.float64_t *> iarr.data
        for i in range(3):
            self.DW[i] = (mesh_octree.DRE[i] - mesh_octree.DLE[i])
            self.periodicity[i] = periodicity[i]
        for i in range(positions.shape[0]):
            for j in range(3):
                pos[j] = positions[i, j]
            oct = particle_octree.get(pos)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            # Note that this has to be our local index, not our in-file index.
            # This is the particle count, which we'll use once we have sorted
            # the particles to calculate the offsets into each oct's particles.
            offset = oct.domain_ind - moff_p
            pcount[offset] += 1
            pdoms[i] = offset # We store the *actual* offset.
            oct = mesh_octree.get(pos)
            offset = oct.domain_ind - moff_m
        # Now we have oct assignments.  Let's sort them.
        # Note that what we will be providing to our processing functions will
        # actually be indirectly-sorted fields.  This preserves memory at the
        # expense of additional pointer lookups.
        pind = np.argsort(pdoms)
        pind = np.asarray(pind, dtype='int64', order='C')
        # So what this means is that we now have all the oct-0 particle indices
        # in order, then the oct-1, etc etc.
        # This now gives us the indices to the particles for each domain.
        for i in range(positions.shape[0]):
            # This value, poff, is the index of the particle in the *unsorted*
            # arrays.
            poff = pind[i] 
            offset = pdoms[poff] 
            # If we have yet to assign the starting index to this oct, we do so
            # now.
            if doff[offset] < 0: doff[offset] = i
        #print domain_id, domain_offset, moff_p, moff_m
        #raise RuntimeError
        # Now doff is full of offsets to the first entry in the pind that
        # refers to that oct's particles.
        ppos = <np.float64_t *> positions.data
        cart_pos = <np.float64_t *> cart_positions.data
        doffs = <np.int64_t*> doff.data
        pinds = <np.int64_t*> pind.data
        pcounts = <np.int64_t*> pcount.data
        cdef np.ndarray[np.uint8_t, ndim=1] visited
        visited = np.zeros(mdom_ind.shape[0], dtype="uint8")
        cdef int nproc = 0
        for i in range(oct_positions.shape[0]):
            for j in range(3):
                pos[j] = oct_positions[i, j]
            oct = mesh_octree.get(pos, &moi)
            offset = mdom_ind[oct.domain_ind - moff_m] * nz
            if visited[oct.domain_ind - moff_m] == 1: continue
            visited[oct.domain_ind - moff_m] = 1
            if offset < 0: continue
            nproc += 1
            self.neighbor_process(dims, moi.left_edge, moi.dds,
                         cart_pos, field_pointers, doffs, &nind,
                         pinds, pcounts, offset, index_field_pointers,
                         particle_octree, domain_id, &nsize)
        #print "VISITED", visited.sum(), visited.size,
        #print 100.0*float(visited.sum())/visited.size
        if nind != NULL:
            free(nind)

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_particles(self, OctreeContainer particle_octree,
                     np.ndarray[np.int64_t, ndim=1] pdom_ind,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None, int domain_id = -1,
                     int domain_offset = 0,
                     periodicity = (True, True, True),
                     geometry = "cartesian"):
        # The other functions in this base class process particles in a way
        # that results in a modification to the *mesh*.  This function is
        # designed to process neighboring particles in such a way that a new
        # *particle* field is defined -- this means that new particle
        # attributes (*not* mesh attributes) can be created that rely on the
        # values of nearby particles.  For instance, a smoothing kernel, or a
        # nearest-neighbor field.
        cdef int nf, i, j, k, dims[3], n
        cdef np.float64_t **field_pointers, *field_vals, pos[3], *ppos, dds[3]
        cdef np.float64_t **octree_field_pointers
        cdef int nsize = 0
        cdef np.int64_t *nind = NULL
        cdef OctInfo moi, poi
        cdef Oct *oct, **neighbors = NULL
        cdef np.int64_t nneighbors, numpart, offset, local_ind
        cdef np.int64_t moff_p, moff_m, pind0
        cdef np.int64_t *doffs, *pinds, *pcounts, poff
        cdef np.ndarray[np.int64_t, ndim=1] pind, doff, pdoms, pcount
        cdef np.ndarray[np.float64_t, ndim=1] tarr
        cdef np.ndarray[np.float64_t, ndim=2] cart_positions
        if geometry == "cartesian":
            self.pos_setup = cart_coord_setup
            cart_positions = positions
        elif geometry == "spherical":
            self.pos_setup = spherical_coord_setup
            cart_positions = np.empty((positions.shape[0], 3), dtype="float64")

            cart_positions[:,0] = positions[:,0] * \
                                  np.sin(positions[:,1]) * \
                                  np.cos(positions[:,2])
            cart_positions[:,1] = positions[:,0] * \
                                  np.sin(positions[:,1]) * \
                                  np.sin(positions[:,2])
            cart_positions[:,2] = positions[:,0] * \
                                  np.cos(positions[:,1])
            periodicity = (False, False, False)
        else:
            raise NotImplementedError
        numpart = positions.shape[0]
        pcount = np.zeros_like(pdom_ind)
        doff = np.zeros_like(pdom_ind) - 1
        moff_p = particle_octree.get_domain_offset(domain_id + domain_offset)
        pdoms = np.zeros(positions.shape[0], dtype="int64") - 1
        nf = len(fields)
        if fields is None:
            fields = []
        field_pointers = <np.float64_t**> alloca(sizeof(np.float64_t *) * nf)
        for i in range(nf):
            tarr = fields[i]
            field_pointers[i] = <np.float64_t *> tarr.data
        for i in range(3):
            self.DW[i] = (particle_octree.DRE[i] - particle_octree.DLE[i])
            self.periodicity[i] = periodicity[i]
        for i in range(positions.shape[0]):
            for j in range(3):
                pos[j] = positions[i, j]
            oct = particle_octree.get(pos)
            if oct == NULL or (domain_id > 0 and oct.domain != domain_id):
                continue
            # Note that this has to be our local index, not our in-file index.
            # This is the particle count, which we'll use once we have sorted
            # the particles to calculate the offsets into each oct's particles.
            offset = oct.domain_ind - moff_p
            pcount[offset] += 1
            pdoms[i] = offset # We store the *actual* offset.
        # Now we have oct assignments.  Let's sort them.
        # Note that what we will be providing to our processing functions will
        # actually be indirectly-sorted fields.  This preserves memory at the
        # expense of additional pointer lookups.
        pind = np.argsort(pdoms)
        pind = np.asarray(pind, dtype='int64', order='C')
        # So what this means is that we now have all the oct-0 particle indices
        # in order, then the oct-1, etc etc.
        # This now gives us the indices to the particles for each domain.
        for i in range(positions.shape[0]):
            # This value, poff, is the index of the particle in the *unsorted*
            # arrays.
            poff = pind[i] 
            offset = pdoms[poff] 
            # If we have yet to assign the starting index to this oct, we do so
            # now.
            if doff[offset] < 0: doff[offset] = i
        #print domain_id, domain_offset, moff_p, moff_m
        #raise RuntimeError
        # Now doff is full of offsets to the first entry in the pind that
        # refers to that oct's particles.
        ppos = <np.float64_t *> positions.data
        cart_pos = <np.float64_t *> cart_positions.data
        doffs = <np.int64_t*> doff.data
        pinds = <np.int64_t*> pind.data
        pcounts = <np.int64_t*> pcount.data
        cdef int maxnei = 0
        cdef int nproc = 0
        for i in range(doff.shape[0]):
            if doff[i] < 0: continue
            offset = pind[doff[i]]
            for j in range(3):
                pos[j] = positions[offset, j]
            for j in range(pcount[i]):
                pind0 = pind[doff[i] + j]
                for k in range(3):
                    pos[k] = positions[pind0, k]
                self.neighbor_process_particle(pos, cart_pos, field_pointers,
                            doffs, &nind, pinds, pcounts, pind0,
                            NULL, particle_octree, domain_id, &nsize)
        #print "VISITED", visited.sum(), visited.size,
        #print 100.0*float(visited.sum())/visited.size
        if nind != NULL:
            free(nind)

    cdef int neighbor_search(self, np.float64_t pos[3], OctreeContainer octree,
                             np.int64_t **nind, int *nsize, 
                             np.int64_t nneighbors, np.int64_t domain_id,
                             Oct **oct = NULL):
        cdef OctInfo oi
        cdef Oct *ooct, **neighbors
        cdef int j
        cdef np.int64_t moff = octree.get_domain_offset(domain_id)
        ooct = octree.get(pos, &oi)
        if oct != NULL and ooct == oct[0]:
            return nneighbors
        oct[0] = ooct
        if nind[0] == NULL:
            nsize[0] = 27
            nind[0] = <np.int64_t *> malloc(sizeof(np.int64_t)*nsize[0])
        neighbors = octree.neighbors(&oi, &nneighbors, ooct, self.periodicity)
        # Now we have all our neighbors.  And, we should be set for what
        # else we need to do.
        if nneighbors > nsize[0]:
            nind[0] = <np.int64_t *> realloc(
                nind[0], sizeof(np.int64_t)*nneighbors)
            nsize[0] = nneighbors
        
        for j in range(nneighbors):
            # Particle octree neighbor indices
            nind[0][j] = neighbors[j].domain_ind - moff
            for n in range(j):
                if nind[0][j] == nind[0][n]:
                    nind[0][j] = -1
                break
        # This is allocated by the neighbors function, so we deallocate it.
        free(neighbors)
        return nneighbors
        
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def process_grid(self, gobj,
                     np.ndarray[np.float64_t, ndim=2] positions,
                     fields = None):
        raise NotImplementedError

    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **ifields):
        raise NotImplementedError

    cdef void neighbor_reset(self):
        self.curn = 0
        for i in range(self.maxn):
            self.neighbors[i].pn = -1
            self.neighbors[i].r2 = 1e300

    cdef void neighbor_eval(self, np.int64_t pn, np.float64_t ppos[3],
                            np.float64_t cpos[3]):
        cdef NeighborList *cur
        cdef int i, j
        # _c means candidate (what we're evaluating)
        # _o means other (the item in the list)
        cdef np.float64_t r2_c, r2_o
        cdef np.int64_t pn_c, pn_o
        # If we're less than the maximum number of neighbors, we simply append.
        # After that, we will sort, and then only compare against the rightmost
        # entries.
        if self.curn < self.maxn:
            cur = &self.neighbors[self.curn]
            cur.pn = pn
            cur.r2 = r2dist(ppos, cpos, self.DW, self.periodicity, -1)
            self.curn += 1
            if self.curn == self.maxn:
                # This time we sort it, so that future insertions will be able
                # to be done in order.
                qsort(self.neighbors, self.curn, sizeof(NeighborList), 
                      Neighbor_compare)
            return
        # This will go (curn - 1) through 0.
        r2_o = self.neighbors[self.curn - 1].r2
        r2_c = r2dist(ppos, cpos, self.DW, self.periodicity, r2_o)
        # Early terminate
        if r2_c < 0: return
        pn_c = pn
        for j in range(1, self.maxn + 1):
            i = self.maxn - j
            # First we evaluate against i.  If our candidate radius is greater
            # than the one we're inspecting, we quit.
            cur = &self.neighbors[i]
            r2_o = cur.r2
            pn_o = cur.pn
            if r2_c > r2_o:
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
            if nind[ni] == -1: continue
            offset = doffs[nind[ni]]
            pc = pcounts[nind[ni]]
            for i in range(pc):
                pn = pinds[offset + i]
                for j in range(3):
                    pos[j] = ppos[pn * 3 + j]
                self.neighbor_eval(pn, pos, cpos)

    cdef void neighbor_process(self, int dim[3], np.float64_t left_edge[3],
                               np.float64_t dds[3], np.float64_t *ppos,
                               np.float64_t **fields, 
                               np.int64_t *doffs, np.int64_t **nind,
                               np.int64_t *pinds, np.int64_t *pcounts,
                               np.int64_t offset,
                               np.float64_t **index_fields,
                               OctreeContainer octree, np.int64_t domain_id,
                               int *nsize):
        # Note that we assume that fields[0] == smoothing length in the native
        # units supplied.  We can now iterate over every cell in the block and
        # every particle to find the nearest.  We will use a priority heap.
        cdef int i, j, k, ntot, nntot, m, nneighbors
        cdef np.float64_t cpos[3], opos[3]
        cdef Oct* oct = NULL
        cpos[0] = left_edge[0] + 0.5*dds[0]
        for i in range(dim[0]):
            cpos[1] = left_edge[1] + 0.5*dds[1]
            for j in range(dim[1]):
                cpos[2] = left_edge[2] + 0.5*dds[2]
                for k in range(dim[2]):
                    self.pos_setup(cpos, opos)
                    nneighbors = self.neighbor_search(opos, octree,
                                    nind, nsize, nneighbors, domain_id, &oct)
                    self.neighbor_find(nneighbors, nind[0], doffs, pcounts,
                        pinds, ppos, opos)
                    # Now we have all our neighbors in our neighbor list.
                    if self.curn <-1*self.maxn:
                        ntot = nntot = 0
                        for m in range(nneighbors):
                            if nind[0][m] < 0: continue
                            nntot += 1
                            ntot += pcounts[nind[0][m]]
                        print "SOMETHING WRONG", self.curn, nneighbors, ntot, nntot
                    self.process(offset, i, j, k, dim, opos, fields,
                                 index_fields)
                    cpos[2] += dds[2]
                cpos[1] += dds[1]
            cpos[0] += dds[0]

    cdef void neighbor_process_particle(self, np.float64_t cpos[3],
                               np.float64_t *ppos,
                               np.float64_t **fields,
                               np.int64_t *doffs, np.int64_t **nind,
                               np.int64_t *pinds, np.int64_t *pcounts,
                               np.int64_t offset,
                               np.float64_t **index_fields,
                               OctreeContainer octree, 
                               np.int64_t domain_id, int *nsize):
        # Note that we assume that fields[0] == smoothing length in the native
        # units supplied.  We can now iterate over every cell in the block and
        # every particle to find the nearest.  We will use a priority heap.
        cdef int i, j, k, ntot, nntot, m, dim[3]
        cdef Oct *oct = NULL
        cdef np.int64_t nneighbors = 0
        i = j = k = 0
        dim[0] = dim[1] = dim[2] = 1
        cdef np.float64_t opos[3]
        self.pos_setup(cpos, opos)
        nneighbors = self.neighbor_search(opos, octree,
                        nind, nsize, nneighbors, domain_id, &oct)
        self.neighbor_find(nneighbors, nind[0], doffs, pcounts, pinds, ppos, opos)
        self.process(offset, i, j, k, dim, opos, fields, index_fields)

cdef class VolumeWeightedSmooth(ParticleSmoothOperation):
    cdef np.float64_t **fp
    cdef public object vals
    def initialize(self):
        cdef int i
        if self.nfields < 3:
            # We need four fields -- the mass should be the first, then the
            # smoothing length for particles, the normalization factor to
            # ensure mass conservation, then the field we're smoothing.
            raise RuntimeError
        cdef np.ndarray tarr
        self.fp = <np.float64_t **> malloc(
            sizeof(np.float64_t *) * (self.nfields - 3))
        self.vals = []
        for i in range(self.nfields - 2):
            tarr = np.zeros(self.nvals, dtype="float64", order="F")
            self.vals.append(tarr)
            self.fp[i] = <np.float64_t *> tarr.data

    def finalize(self):
        free(self.fp)
        vv = self.vals.pop(-1)
        for v in self.vals:
            v /= vv
        return self.vals

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields):
        # We have our i, j, k for our cell, as well as the cell position.
        # We also have a list of neighboring particles with particle numbers.
        cdef int n, fi
        cdef np.float64_t weight, r2, val, hsml, dens, mass, coeff, max_r
        coeff = 0.0
        cdef np.int64_t pn
        # We get back our mass 
        # rho_i = sum(j = 1 .. n) m_j * W_ij
        max_r = sqrt(self.neighbors[self.curn-1].r2)
        for n in range(self.curn):
            # No normalization for the moment.
            # fields[0] is the smoothing length.
            r2 = self.neighbors[n].r2
            pn = self.neighbors[n].pn
            # Smoothing kernel weight function
            mass = fields[0][pn]
            hsml = fields[1][pn]
            if hsml < 0:
                hsml = max_r
            if hsml == 0: continue
            # Usually this density has been computed
            dens = fields[2][pn]
            if dens == 0.0: continue
            weight = mass * sph_kernel(sqrt(r2) / hsml) / dens
            # Mass of the particle times the value 
            for fi in range(self.nfields - 3):
                val = fields[fi + 3][pn]
                self.fp[fi][gind(i,j,k,dim) + offset] += val * weight
            self.fp[self.nfields - 3][gind(i,j,k,dim) + offset] += weight
        return

volume_weighted_smooth = VolumeWeightedSmooth

cdef class NearestNeighborSmooth(ParticleSmoothOperation):
    cdef np.float64_t *fp
    cdef public object vals
    def initialize(self):
        cdef np.ndarray tarr
        assert(self.nfields == 1)
        tarr = np.zeros(self.nvals, dtype="float64", order="F")
        self.vals = tarr
        self.fp = <np.float64_t *> tarr.data

    def finalize(self):
        return self.vals

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields):
        # We have our i, j, k for our cell, as well as the cell position.
        # We also have a list of neighboring particles with particle numbers.
        cdef np.int64_t pn
        # We get back our mass 
        # rho_i = sum(j = 1 .. n) m_j * W_ij
        pn = self.neighbors[0].pn
        self.fp[gind(i,j,k,dim) + offset] = fields[0][pn]
        #self.fp[gind(i,j,k,dim) + offset] = self.neighbors[0].r2
        return

nearest_smooth = NearestNeighborSmooth

cdef class IDWInterpolationSmooth(ParticleSmoothOperation):
    cdef np.float64_t *fp
    cdef public int p2
    cdef public object vals
    def initialize(self):
        cdef np.ndarray tarr
        assert(self.nfields == 1)
        tarr = np.zeros(self.nvals, dtype="float64", order="F")
        self.vals = tarr
        self.fp = <np.float64_t *> tarr.data
        self.p2 = 2 # Power, for IDW, in units of 2.  So we only do even p's.

    def finalize(self):
        return self.vals

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields):
        # We have our i, j, k for our cell, as well as the cell position.
        # We also have a list of neighboring particles with particle numbers.
        cdef np.int64_t pn, ni, di
        cdef np.float64_t total_weight = 0.0, total_value = 0.0, r2, val, w
        # We're going to do a very simple IDW average
        if self.neighbors[0].r2 == 0.0:
            pn = self.neighbors[0].pn
            self.fp[gind(i,j,k,dim) + offset] = fields[0][pn]
        for ni in range(self.curn):
            r2 = self.neighbors[ni].r2
            val = fields[0][self.neighbors[ni].pn]
            w = r2
            for di in range(self.p2 - 1):
                w *= r2
            total_value += w * val
            total_weight += w
        self.fp[gind(i,j,k,dim) + offset] = total_value / total_weight
        return

idw_smooth = IDWInterpolationSmooth

cdef class NthNeighborDistanceSmooth(ParticleSmoothOperation):

    def initialize(self):
        return

    def finalize(self):
        return

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields):
        cdef np.float64_t max_r
        # We assume "offset" here is the particle index.
        max_r = sqrt(self.neighbors[self.curn-1].r2)
        fields[0][offset] = max_r

nth_neighbor_smooth = NthNeighborDistanceSmooth

cdef class SmoothedDensityEstimate(ParticleSmoothOperation):
    def initialize(self):
        return

    def finalize(self):
        return

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void process(self, np.int64_t offset, int i, int j, int k,
                      int dim[3], np.float64_t cpos[3], np.float64_t **fields,
                      np.float64_t **index_fields):
        cdef np.float64_t r2, hsml, dens, mass, weight, lw
        cdef int pn
        # We assume "offset" here is the particle index.
        hsml = sqrt(self.neighbors[self.curn-1].r2)
        dens = 0.0
        weight = 0.0
        for pn in range(self.curn):
            mass = fields[0][self.neighbors[pn].pn]
            r2 = self.neighbors[pn].r2
            lw = sph_kernel(sqrt(r2) / hsml)
            dens += mass * lw
        weight = (4.0/3.0) * 3.1415926 * hsml**3
        fields[1][offset] = dens/weight

density_smooth = SmoothedDensityEstimate
