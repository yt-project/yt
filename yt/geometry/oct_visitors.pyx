"""
Oct visitor functions




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy
import numpy
from yt.utilities.lib.fp_utils cimport *
from libc.stdlib cimport malloc, free
from yt.geometry.oct_container cimport OctreeContainer

# Now some visitor functions

cdef class OctVisitor:
    def __init__(self, OctreeContainer octree, int domain_id = -1):
        cdef int i
        self.index = 0
        self.last = -1
        self.global_index = -1
        for i in range(3):
            self.pos[i] = -1
            self.ind[i] = -1
        self.dims = 0
        self.domain = domain_id
        self.level = -1
        self.oref = octree.oref
        self.nz = (1 << (self.oref*3))

    cdef void visit(self, Oct* o, np.uint8_t selected):
        raise NotImplementedError

# This copies an integer array from the source to the destination, based on the
# selection criteria.
cdef class CopyArrayI64(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We should always have global_index less than our source.
        # "last" here tells us the dimensionality of the array.
        if selected == 0: return
        # There are this many records between "octs"
        self.dest[self.index, :] = self.source[
                self.ind[2], self.ind[1], self.ind[0],
                self.global_index, :]
        self.index += 1

# This copies a floating point array from the source to the destination, based
# on the selection criteria.
cdef class CopyArrayF64(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We should always have global_index less than our source.
        # "last" here tells us the dimensionality of the array.
        if selected == 0: return
        # There are this many records between "octs"
        self.dest[self.index, :] = self.source[
                self.ind[2], self.ind[1], self.ind[0],
                self.global_index, :]
        self.index += 1

# This counts the number of octs, selected or not, that the selector hits.
# Note that the selector will not recursively visit unselected octs, so this is
# still useful.
cdef class CountTotalOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Count even if not selected.
        # Number of *octs* visited.
        if self.last != o.domain_ind:
            self.index += 1
            self.last = o.domain_ind

# This counts the number of selected cells.
cdef class CountTotalCells(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Number of *cells* visited and selected.
        self.index += selected

# Every time a cell is visited, mark it.  This will be for all visited octs.
cdef class MarkOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We mark them even if they are not selected
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            self.index += 1
        self.mark[self.index, self.ind[2], self.ind[1], self.ind[0]] = 1

# Mask all the selected cells.
cdef class MaskOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        self.mask[self.global_index, self.ind[2], self.ind[1], self.ind[0]] = 1

# Compute a mapping from domain_ind to flattened index.
cdef class IndexOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that we provide an index even if the cell is not selected.
        if self.last != o.domain_ind:
            self.last = o.domain_ind
            self.oct_index[o.domain_ind] = self.index
            self.index += 1

# Integer coordinates
cdef class ICoordsOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        cdef int i
        for i in range(3):
            self.icoords[self.index,i] = (self.pos[i] << self.oref) + self.ind[i]
        self.index += 1

# Level
cdef class IResOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        self.ires[self.index] = self.level
        self.index += 1

# Floating point coordinates
cdef class FCoordsOcts(OctVisitor):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that this does not actually give the correct floating point
        # coordinates.  It gives them in some unit system where the domain is 1.0
        # in all directions, and assumes that they will be scaled later.
        if selected == 0: return
        cdef int i
        cdef np.float64_t c, dx
        dx = 1.0 / ((1 << self.oref) << self.level)
        for i in range(3):
            c = <np.float64_t> ((self.pos[i] << self.oref ) + self.ind[i])
            self.fcoords[self.index,i] = (c + 0.5) * dx
        self.index += 1

# Floating point widths; domain modifications are done later.
cdef class FWidthOcts(OctVisitor):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # Note that this does not actually give the correct floating point
        # coordinates.  It gives them in some unit system where the domain is 1.0
        # in all directions, and assumes that they will be scaled later.
        if selected == 0: return
        cdef int i
        cdef np.float64_t dx
        dx = 1.0 / ((1 << self.oref) << self.level)
        for i in range(3):
            self.fwidth[self.index,i] = dx
        self.index += 1

# Mark which domains are touched by a selector.
cdef class IdentifyOcts(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We assume that our domain has *already* been selected by, which means
        # we'll get all cells within the domain for a by-domain selector and all
        # cells within the domain *and* selector for the selector itself.
        if selected == 0: return
        self.domain_mask[o.domain - 1] = 1

# Assign domain indices to octs
cdef class AssignDomainInd(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        o.domain_ind = self.global_index
        self.index += 1

# From the file, fill in C order
cdef class FillFileIndicesO(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We fill these arrays, then inside the level filler we use these as
        # indices as we fill a second array from the self.
        if selected == 0: return
        self.levels[self.index] = self.level
        self.file_inds[self.index] = o.file_ind
        self.cell_inds[self.index] = self.oind()
        self.index +=1

# From the file, fill in F order
cdef class FillFileIndicesR(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        # We fill these arrays, then inside the level filler we use these as
        # indices as we fill a second array from the self.
        if selected == 0: return
        self.levels[self.index] = self.level
        self.file_inds[self.index] = o.file_ind
        self.cell_inds[self.index] = self.rind()
        self.index +=1

# Count octs by domain
cdef class CountByDomain(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        if selected == 0: return
        # NOTE: We do this for every *cell*.
        self.domain_counts[o.domain - 1] += 1

# Store the refinement mapping of the octree to be loaded later
cdef class StoreOctree(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef np.uint8_t res, ii
        ii = cind(self.ind[0], self.ind[1], self.ind[2])
        if o.children == NULL:
            # Not refined.
            res = 0
        else:
            res = 1
        self.ref_mask[self.index] = res
        self.index += 1

# Go from a refinement mapping to a new octree
cdef class LoadOctree(OctVisitor):
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    cdef void visit(self, Oct* o, np.uint8_t selected):
        cdef int i, ii
        ii = cind(self.ind[0], self.ind[1], self.ind[2])
        if self.ref_mask[self.index] == 0:
            # We only want to do this once.  Otherwise we end up with way too many
            # nfinest for our tastes.
            if o.file_ind == -1:
                o.children = NULL
                o.file_ind = self.nfinest[0]
                o.domain = 1
                self.nfinest[0] += 1
        elif self.ref_mask[self.index] > 0:
            if self.ref_mask[self.index] != 1 and self.ref_mask[self.index] != 8:
                print "ARRAY CLUE: ", self.ref_mask[self.index], "UNKNOWN"
                raise RuntimeError
            if o.children == NULL:
                o.children = <Oct **> malloc(sizeof(Oct *) * 8)
                for i in range(8):
                    o.children[i] = NULL
            for i in range(8):
                o.children[ii + i] = &self.octs[self.nocts[0]]
                o.children[ii + i].domain_ind = self.nocts[0]
                o.children[ii + i].file_ind = -1
                o.children[ii + i].domain = -1
                o.children[ii + i].children = NULL
                self.nocts[0] += 1
        else:
            print "SOMETHING IS AMISS", self.index
            raise RuntimeError
        self.index += 1
