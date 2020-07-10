# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Geometry selection routines.




"""


import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
from cython cimport floating
from libc.stdlib cimport malloc, free
from yt.utilities.lib.fnv_hash cimport c_fnv_hash as fnv_hash
from yt.utilities.lib.fp_utils cimport fclip, iclip, fmax, fmin, imin, imax
from .oct_container cimport OctreeContainer, Oct
cimport oct_visitors
from .oct_visitors cimport cind
from yt.utilities.lib.volume_container cimport \
    VolumeContainer
from yt.utilities.lib.grid_traversal cimport \
    sampler_function, walk_volume
from yt.utilities.lib.bitarray cimport ba_get_value, ba_set_value
from yt.utilities.lib.geometry_utils cimport encode_morton_64bit, decode_morton_64bit, \
    bounded_morton_dds, morton_neighbors_coarse, morton_neighbors_refined

cdef extern from "math.h":
    double exp(double x) nogil
    float expf(float x) nogil
    long double expl(long double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil
    double fmod(double x, double y) nogil
    double log2(double x) nogil
    long int lrint(double x) nogil
    double fabs(double x) nogil

# use this as an epsilon test for grids aligned with selector
# define here to avoid the gil later
cdef np.float64_t grid_eps = np.finfo(np.float64).eps
grid_eps = 0.0

cdef inline np.float64_t dot(np.float64_t* v1,
                             np.float64_t* v2) nogil:
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

cdef inline np.float64_t norm(np.float64_t* v) nogil:
    return sqrt(dot(v, v))

# These routines are separated into a couple different categories:
#
#   * Routines for identifying intersections of an object with a bounding box
#   * Routines for identifying cells/points inside a bounding box that
#     intersect with an object
#   * Routines that speed up some type of geometric calculation

# First, bounding box / object intersection routines.
# These all respect the interface "dobj" and a set of left_edges, right_edges,
# sometimes also accepting level and mask information.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def convert_mask_to_indices(np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
            int count, int transpose = 0):
    cdef int i, j, k, cpos
    cdef np.ndarray[np.int64_t, ndim=2] indices
    indices = np.zeros((count, 3), dtype='int64')
    cpos = 0
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i, j, k] == 1:
                    if transpose == 1:
                        indices[cpos, 0] = k
                        indices[cpos, 1] = j
                        indices[cpos, 2] = i
                    else:
                        indices[cpos, 0] = i
                        indices[cpos, 1] = j
                        indices[cpos, 2] = k
                    cpos += 1
    return indices


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _mask_fill(np.ndarray[np.float64_t, ndim=1] out,
                np.int64_t offset,
                np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
                np.ndarray[floating, ndim=3] vals):
    cdef np.int64_t count = 0
    cdef int i, j, k
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            for k in range(mask.shape[2]):
                if mask[i, j, k] == 1:
                    out[offset + count] = vals[i, j, k]
                    count += 1
    return count

def mask_fill(np.ndarray[np.float64_t, ndim=1] out,
              np.int64_t offset,
              np.ndarray[np.uint8_t, ndim=3, cast=True] mask,
              np.ndarray vals):
    if vals.dtype == np.float32:
        return _mask_fill[np.float32_t](out, offset, mask, vals)
    elif vals.dtype == np.float64:
        return _mask_fill[np.float64_t](out, offset, mask, vals)
    else:
        raise RuntimeError

cdef class SelectorObject:

    def __cinit__(self, dobj, *args):
        self._hash_initialized = 0
        cdef np.float64_t [:] DLE
        cdef np.float64_t [:] DRE
        min_level = getattr(dobj, "min_level", None)
        max_level = getattr(dobj, "max_level", None)
        if min_level is None:
            min_level = 0
        if max_level is None:
            max_level = 99
        self.min_level = min_level
        self.max_level = max_level
        self.overlap_cells = 0

        ds = getattr(dobj, 'ds', None)
        if ds is None:
            for i in range(3):
                # NOTE that this is not universal.
                self.domain_width[i] = 1.0
                self.periodicity[i] = False
        else:
            DLE = _ensure_code(ds.domain_left_edge)
            DRE = _ensure_code(ds.domain_right_edge)
            for i in range(3):
                self.domain_width[i] = DRE[i] - DLE[i]
                self.periodicity[i] = ds.periodicity[i]

    def get_periodicity(self):
        cdef int i
        cdef np.ndarray[np.uint8_t, ndim=1] periodicity
        periodicity = np.zeros(3, dtype='uint8')
        for i in range(3):
            periodicity[i] = self.periodicity[i]
        return periodicity

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int i, n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        _ensure_code(left_edges)
        _ensure_code(right_edges)
        with nogil:
            for n in range(ng):
                # Call our selector function
                # Check if the sphere is inside the grid
                for i in range(3):
                    LE[i] = left_edges[n, i]
                    RE[i] = right_edges[n, i]
                gridi[n] = self.select_grid(LE, RE, levels[n, 0])
        return gridi.astype("bool")

    def count_octs(self, OctreeContainer octree, int domain_id = -1):
        cdef oct_visitors.CountTotalOcts visitor
        visitor = oct_visitors.CountTotalOcts(octree, domain_id)
        octree.visit_all_octs(self, visitor)
        return visitor.index

    def count_oct_cells(self, OctreeContainer octree, int domain_id = -1):
        cdef oct_visitors.CountTotalCells visitor
        visitor = oct_visitors.CountTotalCells(octree, domain_id)
        octree.visit_all_octs(self, visitor)
        return visitor.index

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void recursively_visit_octs(self, Oct *root,
                        np.float64_t pos[3], np.float64_t dds[3],
                        int level,
                        OctVisitor visitor,
                        int visit_covered = 0):
        # visit_covered tells us whether this octree supports partial
        # refinement.  If it does, we need to handle this specially -- first
        # we visit *this* oct, then we make a second pass to check any child
        # octs.
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t sdds[3]
        cdef np.float64_t spos[3]
        cdef int i, j, k, res
        cdef Oct *ch
        # Remember that pos is the *center* of the oct, and dds is the oct
        # width.  So to get to the edges, we add/subtract half of dds.
        for i in range(3):
            # sdds is the cell width
            sdds[i] = dds[i]/2.0
            LE[i] = pos[i] - dds[i]/2.0
            RE[i] = pos[i] + dds[i]/2.0
        #print(LE[0], RE[0], LE[1], RE[1], LE[2], RE[2])
        res = self.select_grid(LE, RE, level, root)
        if res == 1 and visitor.domain > 0 and root.domain != visitor.domain:
            res = -1
        cdef int increment = 1
        cdef int next_level, this_level
        # next_level: an int that says whether or not we can progress to children
        # this_level: an int that says whether or not we can select from this
        # level
        next_level = this_level = 1
        if res == -1:
            # This happens when we do domain selection but the oct has
            # children.  This would allow an oct to pass to its children but
            # not get accessed itself.
            next_level = 1
            this_level = 0
        elif level == self.max_level:
            next_level = 0
        elif level < self.min_level or level > self.max_level:
            this_level = 0
        if res == 0 and this_level == 1:
            return
        # Now we visit all our children.  We subtract off sdds for the first
        # pass because we center it on the first cell.
        cdef int iter = 1 - visit_covered # 2 if 1, 1 if 0.
        # So the order here goes like so.  If visit_covered is 1, which usually
        # comes from "partial_coverage", we visit the components of a zone even
        # if it has children.  But in general, the first iteration through, we
        # visit each cell.  This means that only if visit_covered is true do we
        # visit potentially covered cells.  The next time through, we visit
        # child cells.
        while iter < 2:
            spos[0] = pos[0] - sdds[0]/2.0
            for i in range(2):
                spos[1] = pos[1] - sdds[1]/2.0
                for j in range(2):
                    spos[2] = pos[2] - sdds[2]/2.0
                    for k in range(2):
                        ch = NULL
                        # We only supply a child if we are actually going to
                        # look at the next level.
                        if root.children != NULL and next_level == 1:
                            ch = root.children[cind(i, j, k)]
                        if iter == 1 and next_level == 1 and ch != NULL:
                            # Note that visitor.pos is always going to be the
                            # position of the Oct -- it is *not* always going
                            # to be the same as the position of the cell under
                            # investigation.
                            visitor.pos[0] = (visitor.pos[0] << 1) + i
                            visitor.pos[1] = (visitor.pos[1] << 1) + j
                            visitor.pos[2] = (visitor.pos[2] << 1) + k
                            visitor.level += 1
                            self.recursively_visit_octs(
                                ch, spos, sdds, level + 1, visitor,
                                visit_covered)
                            visitor.pos[0] = (visitor.pos[0] >> 1)
                            visitor.pos[1] = (visitor.pos[1] >> 1)
                            visitor.pos[2] = (visitor.pos[2] >> 1)
                            visitor.level -= 1
                        elif this_level == 1 and visitor.oref > 0:
                            visitor.global_index += increment
                            increment = 0
                            self.visit_oct_cells(root, ch, spos, sdds,
                                                 visitor, i, j, k)
                        elif this_level == 1 and increment == 1:
                            visitor.global_index += increment
                            increment = 0
                            visitor.ind[0] = visitor.ind[1] = visitor.ind[2] = 0
                            visitor.visit(root, 1)
                        spos[2] += sdds[2]
                    spos[1] += sdds[1]
                spos[0] += sdds[0]
            this_level = 0 # We turn this off for the second pass.
            iter += 1

    cdef void visit_oct_cells(self, Oct *root, Oct *ch,
                              np.float64_t spos[3], np.float64_t sdds[3],
                              OctVisitor visitor, int i, int j, int k):
        # We can short-circuit the whole process if data.oref == 1.
        # This saves us some funny-business.
        cdef int selected
        if visitor.oref == 1:
            selected = self.select_cell(spos, sdds)
            if ch != NULL:
                selected *= self.overlap_cells
            # visitor.ind refers to the cell, not to the oct.
            visitor.ind[0] = i
            visitor.ind[1] = j
            visitor.ind[2] = k
            visitor.visit(root, selected)
            return
        # Okay, now that we've got that out of the way, we have to do some
        # other checks here.  In this case, spos[] is the position of the
        # center of a *possible* oct child, which means it is the center of a
        # cluster of cells.  That cluster might have 1, 8, 64, ... cells in it.
        # But, we can figure it out by calculating the cell dds.
        cdef np.float64_t dds[3]
        cdef np.float64_t pos[3]
        cdef int ci, cj, ck
        cdef int nr = (1 << (visitor.oref - 1))
        for ci in range(3):
            dds[ci] = sdds[ci] / nr
        # Boot strap at the first index.
        pos[0] = (spos[0] - sdds[0]/2.0) + dds[0] * 0.5
        for ci in range(nr):
            pos[1] = (spos[1] - sdds[1]/2.0) + dds[1] * 0.5
            for cj in range(nr):
                pos[2] = (spos[2] - sdds[2]/2.0) + dds[2] * 0.5
                for ck in range(nr):
                    selected = self.select_cell(pos, dds)
                    if ch != NULL:
                        selected *= self.overlap_cells
                    visitor.ind[0] = ci + i * nr
                    visitor.ind[1] = cj + j * nr
                    visitor.ind[2] = ck + k * nr
                    visitor.visit(root, selected)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3],
                               np.int32_t level, Oct *o = NULL) nogil:
        if level < self.min_level or level > self.max_level: return 0
        return self.select_bbox(left_edge, right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid_edge(self, np.float64_t left_edge[3],
                                    np.float64_t right_edge[3],
                                    np.int32_t level, Oct *o = NULL) nogil:
        if level < self.min_level or level > self.max_level: return 0
        return self.select_bbox_edge(left_edge, right_edge)

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 0

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        """
        Returns:
          0: If the selector does not touch the bounding box.
          1: If the selector overlaps the bounding box anywhere.
        """
        return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        """
        Returns:
          0: If the selector does not touch the bounding box.
          1: If the selector contains the entire bounding box.
          2: If the selector contains part of the bounding box.
        """
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.float64_t periodic_difference(self, np.float64_t x1, np.float64_t x2, int d) nogil:
        # domain_width is already in code units, and we assume what is fed in
        # is too.
        cdef np.float64_t rel = x1 - x2
        if self.periodicity[d]:
            if rel > self.domain_width[d] * 0.5:
                rel -= self.domain_width[d]
            elif rel < -self.domain_width[d] * 0.5:
                rel += self.domain_width[d]
        return rel

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mesh_mask(self, mesh):
        cdef np.float64_t pos[3]
        cdef np.ndarray[np.int64_t, ndim=2] indices
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef int i, j, k, selected
        cdef int npoints, nv = mesh._connectivity_length
        cdef int total = 0
        cdef int offset = mesh._index_offset
        coords = _ensure_code(mesh.connectivity_coords)
        indices = mesh.connectivity_indices
        npoints = indices.shape[0]
        mask = np.zeros(npoints, dtype='uint8')
        for i in range(npoints):
            selected = 0
            for j in range(nv):
                for k in range(3):
                    pos[k] = coords[indices[i, j] - offset, k]
                selected = self.select_point(pos)
                if selected == 1: break
            total += selected
            mask[i] = selected
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mesh_cell_mask(self, mesh):
        cdef np.float64_t pos
        cdef np.float64_t le[3]
        cdef np.float64_t re[3]
        cdef np.ndarray[np.int64_t, ndim=2] indices
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef int i, j, k, selected
        cdef int npoints, nv = mesh._connectivity_length
        cdef int ndim = mesh.connectivity_coords.shape[1]
        cdef int total = 0
        cdef int offset = mesh._index_offset
        coords = _ensure_code(mesh.connectivity_coords)
        indices = mesh.connectivity_indices
        npoints = indices.shape[0]
        mask = np.zeros(npoints, dtype='uint8')
        for i in range(npoints):
            selected = 0
            for k in range(ndim):
                le[k] = 1e60
                re[k] = -1e60
            for j in range(nv):
                for k in range(ndim):
                    pos = coords[indices[i, j] - offset, k]
                    le[k] = fmin(pos, le[k])
                    re[k] = fmax(pos, re[k])
            selected = self.select_bbox(le, re)
            total += selected
            mask[i] = selected
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        child_mask = gobj.child_mask
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        cdef int dim[3]
        _ensure_code(gobj.dds)
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        cdef np.ndarray[np.float64_t, ndim=1] odds = gobj.dds.d
        cdef np.ndarray[np.float64_t, ndim=1] oleft_edge = gobj.LeftEdge.d
        cdef np.ndarray[np.float64_t, ndim=1] oright_edge = gobj.RightEdge.d
        cdef int i
        cdef np.float64_t dds[3]
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        for i in range(3):
            dds[i] = odds[i]
            dim[i] = gobj.ActiveDimensions[i]
            left_edge[i] = oleft_edge[i]
            right_edge[i] = oright_edge[i]
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        # Check for the level bounds
        cdef np.int32_t level = gobj.Level
        # We set this to 1 if we ignore child_mask
        cdef int total
        total = self.fill_mask_selector(left_edge, right_edge, dds, dim,
                                        child_mask, mask, level)
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_mask_selector(self, np.float64_t left_edge[3],
                                np.float64_t right_edge[3],
                                np.float64_t dds[3], int dim[3],
                                np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask,
                                np.ndarray[np.uint8_t, ndim=3] mask,
                                int level):
        cdef int i, j, k
        cdef int total = 0, this_level = 0
        cdef np.float64_t pos[3]
        if level < self.min_level or level > self.max_level:
            return 0
        if level == self.max_level:
            this_level = 1
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            for i in range(dim[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                for j in range(dim[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    for k in range(dim[2]):
                        if child_mask[i, j, k] == 1 or this_level == 1:
                            mask[i, j, k] = self.select_cell(pos, dds)
                            total += mask[i, j, k]
                        pos[2] += dds[2]
                    pos[1] += dds[1]
                pos[0] += dds[0]
        return total

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void visit_grid_cells(self, GridVisitorData *data,
                              grid_visitor_function *func,
                              np.uint8_t *cached_mask = NULL):
        # This function accepts a grid visitor function, the data that
        # corresponds to the current grid being examined (the most important
        # aspect of which is the .grid attribute, along with index values and
        # void* pointers to arrays) and a possibly-pre-generated cached mask.
        # Each cell is visited with the grid visitor function.
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef np.float64_t dds[3]
        cdef int dim[3]
        cdef int this_level = 0, level, i
        cdef np.float64_t pos[3]
        level = data.grid.level
        if level < self.min_level or level > self.max_level:
            return
        if level == self.max_level:
            this_level = 1
        cdef np.uint8_t child_masked, selected
        for i in range(3):
            left_edge[i] = data.grid.left_edge[i]
            right_edge[i] = data.grid.right_edge[i]
            dds[i] = (right_edge[i] - left_edge[i])/data.grid.dims[i]
            dim[i] = data.grid.dims[i]
        with nogil:
            pos[0] = left_edge[0] + dds[0] * 0.5
            data.pos[0] = 0
            for i in range(dim[0]):
                pos[1] = left_edge[1] + dds[1] * 0.5
                data.pos[1] = 0
                for j in range(dim[1]):
                    pos[2] = left_edge[2] + dds[2] * 0.5
                    data.pos[2] = 0
                    for k in range(dim[2]):
                        # We short-circuit if we have a cache; if we don't, we
                        # only set selected to true if it's *not* masked by a
                        # child and it *is* selected.
                        if cached_mask != NULL:
                            selected = ba_get_value(cached_mask,
                                                    data.global_index)
                        else:
                            if this_level == 1:
                                child_masked = 0
                            else:
                                child_masked = check_child_masked(data)
                            if child_masked == 0:
                                selected = self.select_cell(pos, dds)
                            else:
                                selected = 0
                        func(data, selected)
                        data.global_index += 1
                        pos[2] += dds[2]
                        data.pos[2] += 1
                    pos[1] += dds[1]
                    data.pos[1] += 1
                pos[0] += dds[0]
                data.pos[0] += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def count_points(self, np.ndarray[floating, ndim=1] x,
                           np.ndarray[floating, ndim=1] y,
                           np.ndarray[floating, ndim=1] z,
                           radii):
        cdef int count = 0
        cdef int i
        cdef np.float64_t pos[3]
        cdef np.float64_t radius
        cdef np.float64_t[:] _radii
        if radii is not None:
            _radii = np.atleast_1d(np.array(radii, dtype='float64'))
        else:
            _radii = np.array([0.0], dtype='float64')
        _ensure_code(x)
        _ensure_code(y)
        _ensure_code(z)
        with nogil:
            for i in range(x.shape[0]):
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                if _radii.shape[0] == 1:
                    radius = _radii[0]
                else:
                    radius = _radii[i]
                if radius == 0:
                    count += self.select_point(pos)
                else:
                    count += self.select_sphere(pos, radius)
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_points(self,
                      np.ndarray[floating, ndim=1] x,
                      np.ndarray[floating, ndim=1] y,
                      np.ndarray[floating, ndim=1] z,
                      radii):
        cdef int count = 0
        cdef int i
        cdef np.float64_t pos[3]
        cdef np.float64_t radius
        cdef np.ndarray[np.uint8_t, ndim=1] mask
        cdef np.float64_t[:] _radii
        if radii is not None:
            _radii = np.atleast_1d(np.array(radii, dtype='float64'))
        else:
            _radii = np.array([0.0], dtype='float64')
        mask = np.empty(x.shape[0], dtype='uint8')
        _ensure_code(x)
        _ensure_code(y)
        _ensure_code(z)


        # this is to allow selectors to optimize the point vs
        # 0-radius sphere case.  These two may have different
        # effects for 0-volume selectors, however (collision
        # between a ray and a point is null, while ray and a
        # sphere is allowed)
        with nogil:
            for i in range(x.shape[0]) :
                pos[0] = x[i]
                pos[1] = y[i]
                pos[2] = z[i]
                if _radii.shape[0] == 1:
                    radius = 0
                else:
                    radius = _radii[i]
                if radius == 0:
                    mask[i] = self.select_point(pos)
                else:
                    mask[i] = self.select_sphere(pos, radius)
                count += mask[i]
        if count == 0: return None
        return mask.view("bool")

    def __hash__(self):
        # convert data to be hashed to a byte array, which FNV algorithm expects
        if self._hash_initialized == 1:
            return self._hash
        hash_data = bytearray()
        for v in self._hash_vals() + self._base_hash():
            if isinstance(v, tuple):
                hash_data.extend(v[0].encode('ascii'))
                hash_data.extend(repr(v[1]).encode('ascii'))
            else:
                hash_data.extend(repr(v).encode('ascii'))
        cdef np.int64_t hash_value = fnv_hash(hash_data)
        self._hash = hash_value
        self._hash_initialized = 1
        return hash_value

    def _hash_vals(self):
        raise NotImplementedError

    def _base_hash(self):
        return (("min_level", self.min_level),
                ("max_level", self.max_level),
                ("overlap_cells", self.overlap_cells),
                ("periodicity[0]", self.periodicity[0]),
                ("periodicity[1]", self.periodicity[1]),
                ("periodicity[2]", self.periodicity[2]),
                ("domain_width[0]", self.domain_width[0]),
                ("domain_width[1]", self.domain_width[1]),
                ("domain_width[2]", self.domain_width[2]))


cdef class PointSelector(SelectorObject):
    cdef np.float64_t p[3]

    def __init__(self, dobj):
        cdef np.float64_t[:] DLE = _ensure_code(dobj.ds.domain_left_edge)
        cdef np.float64_t[:] DRE = _ensure_code(dobj.ds.domain_right_edge)
        for i in range(3):
            self.p[i] = _ensure_code(dobj.p[i])

            # ensure the point lies in the domain
            if self.periodicity[i]:
                self.p[i] = np.fmod(self.p[i], self.domain_width[i])
                if self.p[i] < DLE[i]:
                    self.p[i] += self.domain_width[i]
                elif self.p[i] >= DRE[i]:
                    self.p[i] -= self.domain_width[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if (pos[0] - 0.5*dds[0] <= self.p[0] < pos[0]+0.5*dds[0] and
            pos[1] - 0.5*dds[1] <= self.p[1] < pos[1]+0.5*dds[1] and
            pos[2] - 0.5*dds[2] <= self.p[2] < pos[2]+0.5*dds[2]):
            return 1
        else:
            return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.p[i], i)
            dist2 += dist*dist
        if dist2 <= radius*radius: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # point definitely can only be in one cell
        if (left_edge[0] <= self.p[0] < right_edge[0] and
            left_edge[1] <= self.p[1] < right_edge[1] and
            left_edge[2] <= self.p[2] < right_edge[2]):
            return 1
        else:
            return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # point definitely can only be in one cell
        # Return 2 in all cases to indicate that the point only overlaps
        # portion of box
        if (left_edge[0] <= self.p[0] <= right_edge[0] and
            left_edge[1] <= self.p[1] <= right_edge[1] and
            left_edge[2] <= self.p[2] <= right_edge[2]):
            return 2
        else:
            return 0

    def _hash_vals(self):
        return (("p[0]", self.p[0]),
                ("p[1]", self.p[1]),
                ("p[2]", self.p[2]))

point_selector = PointSelector


cdef class SphereSelector(SelectorObject):
    cdef np.float64_t radius
    cdef np.float64_t radius2
    cdef np.float64_t center[3]
    cdef np.float64_t bbox[3][2]
    cdef bint check_box[3]

    def __init__(self, dobj):
        for i in range(3):
            self.center[i] = _ensure_code(dobj.center[i])
        self.radius = _ensure_code(dobj.radius)
        self.radius2 = self.radius * self.radius
        center = _ensure_code(dobj.center)
        for i in range(3):
            self.center[i] = center[i]
            self.bbox[i][0] = self.center[i] - self.radius
            self.bbox[i][1] = self.center[i] + self.radius
            if self.bbox[i][0] < dobj.ds.domain_left_edge[i]:
                self.check_box[i] = False
            elif self.bbox[i][1] > dobj.ds.domain_right_edge[i]:
                self.check_box[i] = False
            else:
                self.check_box[i] = True

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        # sphere center either inside cell or center of cell lies inside sphere
        if (pos[0] - 0.5*dds[0] <= self.center[0] <= pos[0]+0.5*dds[0] and
            pos[1] - 0.5*dds[1] <= self.center[1] <= pos[1]+0.5*dds[1] and
            pos[2] - 0.5*dds[2] <= self.center[2] <= pos[2]+0.5*dds[2]):
            return 1
        return self.select_point(pos)
        # # langmm: added to allow sphere to interesect edge/corner of cell
        # cdef np.float64_t LE[3]
        # cdef np.float64_t RE[3]
        # cdef int i
        # for i in range(3):
        #     LE[i] = pos[i] - 0.5*dds[i]
        #     RE[i] = pos[i] + 0.5*dds[i]
        # return self.select_bbox(LE, RE)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            if self.check_box[i] and \
              (pos[i] < self.bbox[i][0] or
               pos[i] > self.bbox[i][1]):
                return 0
            dist = _periodic_dist(pos[i], self.center[i], self.domain_width[i],
                                  self.periodicity[i])
            dist2 += dist*dist
            if dist2 > self.radius2: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t dist, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.center[i], i)
            dist2 += dist*dist
        dist = self.radius+radius
        if dist2 <= dist*dist: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, dist, edge
        cdef int i
        if (left_edge[0] <= self.center[0] < right_edge[0] and
            left_edge[1] <= self.center[1] < right_edge[1] and
            left_edge[2] <= self.center[2] < right_edge[2]):
            return 1
        for i in range(3):
            if not self.check_box[i]: continue
            if right_edge[i] < self.bbox[i][0] or \
               left_edge[i] > self.bbox[i][1]:
                return 0
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            # Early terminate
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest*closest
            if dist > self.radius2: return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef np.float64_t box_center, relcenter, closest, farthest, cdist, fdist, edge
        cdef int i
        if (left_edge[0] <= self.center[0] <= right_edge[0] and
            left_edge[1] <= self.center[1] <= right_edge[1] and
            left_edge[2] <= self.center[2] <= right_edge[2]):
            fdist = 0
            for i in range(3):
                edge = right_edge[i] - left_edge[i]
                box_center = (right_edge[i] + left_edge[i])/2.0
                relcenter = self.periodic_difference(
                    box_center, self.center[i], i)
                if relcenter >= 0:
                    farthest = relcenter + edge/2.0
                else:
                    farthest = relcenter - edge/2.0
                # farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
                fdist += farthest*farthest
                if fdist >= self.radius2:
                    return 2  # Box extends outside sphere
            return 1  # Box entirely inside sphere
        for i in range(3):
            if not self.check_box[i]: continue
            if right_edge[i] < self.bbox[i][0] or \
               left_edge[i] > self.bbox[i][1]:
                return 0  # Box outside sphere bounding box
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        cdist = 0
        fdist = 0
        for i in range(3):
            # Early terminate
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            if relcenter >= 0:
                farthest = relcenter + edge/2.0
            else:
                farthest = relcenter - edge/2.0
            #farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
            cdist += closest*closest
            fdist += farthest*farthest
            if cdist > self.radius2:
                return 0  # Box does not overlap sphere
        if fdist < self.radius2:
            return 1  # Sphere extends to entirely contain box
        else:
            return 2  # Sphere only partially overlaps box

    def _hash_vals(self):
        return (("radius", self.radius),
                ("radius2", self.radius2),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]))

sphere_selector = SphereSelector

cdef class RegionSelector(SelectorObject):
    cdef np.float64_t left_edge[3]
    cdef np.float64_t right_edge[3]
    cdef np.float64_t right_edge_shift[3]
    cdef public bint is_all_data
    cdef bint loose_selection
    cdef bint check_period[3]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, dobj):
        cdef int i
        # We are modifying dobj.left_edge and dobj.right_edge , so here we will
        # do an in-place conversion of those arrays.
        cdef np.float64_t[:] RE = _ensure_code(dobj.right_edge)
        cdef np.float64_t[:] LE = _ensure_code(dobj.left_edge)
        cdef np.float64_t[:] DW = _ensure_code(dobj.ds.domain_width)
        cdef np.float64_t[:] DLE = _ensure_code(dobj.ds.domain_left_edge)
        cdef np.float64_t[:] DRE = _ensure_code(dobj.ds.domain_right_edge)
        le_all = (np.array(LE) == _ensure_code(dobj.ds.domain_left_edge)).all()
        re_all = (np.array(RE) == _ensure_code(dobj.ds.domain_right_edge)).all()
        # If we have a bounding box, then we should *not* revert to all data
        domain_override = getattr(dobj.ds, '_domain_override', False)
        if le_all and re_all and not domain_override:
            self.is_all_data = True
        else:
            self.is_all_data = False
        cdef np.float64_t region_width[3]
        cdef bint p[3]
        # This is for if we want to include zones that overlap and whose
        # centers are not strictly included.
        self.loose_selection = getattr(dobj, "loose_selection", False)

        for i in range(3):
            self.check_period[i] = False
            region_width[i] = RE[i] - LE[i]
            p[i] = dobj.ds.periodicity[i]
            if region_width[i] <= 0:
                raise RuntimeError(
                    "Region right edge[%s] < left edge: width = %s" % (
                        i, region_width[i]))

        for i in range(3):

            if p[i]:
                # First, we check if any criteria requires a period check,
                # without any adjustments.  This is for short-circuiting the
                # short-circuit of the loop down below in mask filling.
                if LE[i] < DLE[i] or LE[i] > DRE[i] or RE[i] > DRE[i]:
                    self.check_period[i] = True
                # shift so left_edge guaranteed in domain
                if LE[i] < DLE[i]:
                    LE[i] += DW[i]
                    RE[i] += DW[i]
                elif LE[i] > DRE[i]:
                    LE[i] -= DW[i]
                    RE[i] -= DW[i]
            else:
                if LE[i] < DLE[i] or RE[i] > DRE[i]:
                    raise RuntimeError(
                        "Error: yt attempted to read outside the boundaries of "
                        "a non-periodic domain along dimension %s.\n"
                        "Region left edge = %s, Region right edge = %s\n"
                        "Dataset left edge = %s, Dataset right edge = %s\n\n"
                        "This commonly happens when trying to compute ghost cells "
                        "up to the domain boundary. Two possible solutions are to "
                        "load a smaller region that does not border the edge or "
                        "override the periodicity for this dataset." % \
                        (i, dobj.left_edge[i], dobj.right_edge[i],
                         dobj.ds.domain_left_edge[i], dobj.ds.domain_right_edge[i])
                    )
            # Already ensured in code
            self.left_edge[i] = LE[i]
            self.right_edge[i] = RE[i]
            self.right_edge_shift[i] = RE[i] - DW[i]
            if not self.periodicity[i]:
                self.right_edge_shift[i] = -np.inf

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i
        for i in range(3):
            if (right_edge[i] < self.left_edge[i] and \
                left_edge[i] >= self.right_edge_shift[i]) or \
                left_edge[i] >= self.right_edge[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i
        for i in range(3):
            if (right_edge[i] < self.left_edge[i] and \
                left_edge[i] >= self.right_edge_shift[i]) or \
                left_edge[i] >= self.right_edge[i]:
                return 0
        for i in range(3):
            if left_edge[i] < self.right_edge_shift[i]:
                if right_edge[i] >= self.right_edge_shift[i]:
                    return 2
            elif left_edge[i] < self.left_edge[i] or \
                 right_edge[i] >= self.right_edge[i]:
                return 2
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int i
        if self.loose_selection:
            for i in range(3):
                left_edge[i] = pos[i] - dds[i]*0.5
                right_edge[i] = pos[i] + dds[i]*0.5
            return self.select_bbox(left_edge, right_edge)
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        for i in range(3):
            if (self.right_edge_shift[i] <= pos[i] < self.left_edge[i]) or \
               pos[i] >= self.right_edge[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        # adapted from http://stackoverflow.com/a/4579192/1382869
        cdef int i
        cdef np.float64_t p
        cdef np.float64_t r2 = radius**2
        cdef np.float64_t dmin = 0
        for i in range(3):
            if self.periodicity[i] and self.check_period[i]:
                p = pos[i] + self.right_edge_shift[i]
            else:
                p = pos[i]
            if p < self.left_edge[i]:
                dmin += (p - self.left_edge[i])**2
            elif pos[i] > self.right_edge[i]:
                dmin += (p - self.right_edge[i])**2
        return int(dmin <= r2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int fill_mask_selector(self, np.float64_t left_edge[3],
                                np.float64_t right_edge[3],
                                np.float64_t dds[3], int dim[3],
                                np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask,
                                np.ndarray[np.uint8_t, ndim=3] mask,
                                int level):
        cdef int i, j, k
        cdef int total = 0, this_level = 0
        cdef np.float64_t pos[3]
        if level < self.min_level or level > self.max_level:
            return 0
        if level == self.max_level:
            this_level = 1
        cdef np.int64_t si[3]
        cdef np.int64_t ei[3]
        for i in range(3):
            if not self.check_period[i]:
                si[i] = <np.int64_t> ((self.left_edge[i] - left_edge[i])/dds[i])
                ei[i] = <np.int64_t> ((self.right_edge[i] - left_edge[i])/dds[i])
                si[i] = iclip(si[i] - 1, 0, dim[i])
                ei[i] = iclip(ei[i] + 1, 0, dim[i])
            else:
                si[i] = 0
                ei[i] = dim[i]
        with nogil:
            pos[0] = left_edge[0] + (si[0] + 0.5) * dds[0]
            for i in range(si[0], ei[0]):
                pos[1] = left_edge[1] + (si[1] + 0.5) * dds[1]
                for j in range(si[1], ei[1]):
                    pos[2] = left_edge[2] + (si[2] + 0.5) * dds[2]
                    for k in range(si[2], ei[2]):
                        if child_mask[i, j, k] == 1 or this_level == 1:
                            mask[i, j, k] = self.select_cell(pos, dds)
                            total += mask[i, j, k]
                        pos[2] += dds[2]
                    pos[1] += dds[1]
                pos[0] += dds[0]
        return total


    def _hash_vals(self):
        return (("left_edge[0]", self.left_edge[0]),
                ("left_edge[1]", self.left_edge[1]),
                ("left_edge[2]", self.left_edge[2]),
                ("right_edge[0]", self.right_edge[0]),
                ("right_edge[1]", self.right_edge[1]),
                ("right_edge[2]", self.right_edge[2]))

region_selector = RegionSelector

cdef class CutRegionSelector(SelectorObject):
    cdef set _positions
    cdef tuple _conditionals

    def __init__(self, dobj):
        axis_name = dobj.ds.coordinates.axis_name
        positions = np.array([dobj['index', axis_name[0]],
                              dobj['index', axis_name[1]],
                              dobj['index', axis_name[2]]]).T
        self._conditionals = tuple(dobj.conditionals)
        self._positions = set(tuple(position) for position in positions)

    cdef int select_bbox(self,  np.float64_t left_edge[3],
                     np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_bbox_dge(self,  np.float64_t left_edge[3],
                     np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            if (pos[0], pos[1], pos[2]) in self._positions:
                return 1
            else:
                return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    def _hash_vals(self):
        t = ()
        for i, c in enumerate(self._conditionals):
            t += ("conditional[%s]" % i, c)
        return ("conditionals", t)

cut_region_selector = CutRegionSelector

cdef class DiskSelector(SelectorObject):
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t center[3]
    cdef np.float64_t radius, radius2
    cdef np.float64_t height

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
            self.center[i] = _ensure_code(dobj.center[i])
        self.radius = _ensure_code(dobj.radius)
        self.radius2 = self.radius * self.radius
        self.height = _ensure_code(dobj.height)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.periodic_difference(pos[i], self.center[i], i)
            h += temp * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        if fabs(h) <= self.height and r2 <= self.radius2: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef np.float64_t h, d, r2, temp
        cdef int i
        h = d = 0
        for i in range(3):
            temp = self.periodic_difference(pos[i], self.center[i], i)
            h += temp * self.norm_vec[i]
            d += temp*temp
        r2 = (d - h*h)
        d = self.radius+radius
        if fabs(h) <= self.height+radius and r2 <= d*d: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # Until we can get our OBB/OBB intersection correct, disable this.
        return 1
        # cdef np.float64_t *arr[2]
        # cdef np.float64_t pos[3]
        # cdef np.float64_t H, D, R2, temp
        # cdef int i, j, k, n
        # cdef int all_under = 1
        # cdef int all_over = 1
        # cdef int any_radius = 0
        # # A moment of explanation (revised):
        # #    The disk and bounding box collide if any of the following are true:
        # #    1) the center of the disk is inside the bounding box
        # #    2) any corner of the box lies inside the disk
        # #    3) the box spans the plane (!all_under and !all_over) and at least
        # #       one corner is within the cylindrical radius

        # # check if disk center lies inside bbox
        # if left_edge[0] <= self.center[0] <= right_edge[0] and \
        #    left_edge[1] <= self.center[1] <= right_edge[1] and \
        #    left_edge[2] <= self.center[2] <= right_edge[2] :
        #     return 1

        # # check all corners
        # arr[0] = left_edge
        # arr[1] = right_edge
        # for i in range(2):
        #     pos[0] = arr[i][0]
        #     for j in range(2):
        #         pos[1] = arr[j][1]
        #         for k in range(2):
        #             pos[2] = arr[k][2]
        #             H = D = 0
        #             for n in range(3):
        #                 temp = self.difference(pos[n], self.center[n], n)
        #                 H += (temp * self.norm_vec[n])
        #                 D += temp*temp
        #             R2 = (D - H*H)
        #             if R2 < self.radius2 :
        #                 any_radius = 1
        #                 if fabs(H) < self.height: return 1
        #             if H < 0: all_over = 0
        #             if H > 0: all_under = 0
        # if all_over == 0 and all_under == 0 and any_radius == 1: return 1
        # return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # Until we can get our OBB/OBB intersection correct, disable this.
        return 2
        # cdef np.float64_t *arr[2]
        # cdef np.float64_t pos[3], H, D, R2, temp
        # cdef int i, j, k, n
        # cdef int all_under = 1
        # cdef int all_over = 1
        # cdef int any_radius = 0
        # # A moment of explanation (revised):
        # #    The disk and bounding box collide if any of the following are true:
        # #    1) the center of the disk is inside the bounding box
        # #    2) any corner of the box lies inside the disk
        # #    3) the box spans the plane (!all_under and !all_over) and at least
        # #       one corner is within the cylindrical radius

        # # check if disk center lies inside bbox
        # if left_edge[0] <= self.center[0] <= right_edge[0] and \
        #    left_edge[1] <= self.center[1] <= right_edge[1] and \
        #    left_edge[2] <= self.center[2] <= right_edge[2] :
        #     return 1

        # # check all corners
        # arr[0] = left_edge
        # arr[1] = right_edge
        # for i in range(2):
        #     pos[0] = arr[i][0]
        #     for j in range(2):
        #         pos[1] = arr[j][1]
        #         for k in range(2):
        #             pos[2] = arr[k][2]
        #             H = D = 0
        #             for n in range(3):
        #                 temp = self.periodic_difference(
        #                     pos[n], self.center[n], n)
        #                 H += (temp * self.norm_vec[n])
        #                 D += temp*temp
        #             R2 = (D - H*H)
        #             if R2 < self.radius2 :
        #                 any_radius = 1
        #                 if fabs(H) < self.height: return 1
        #             if H < 0: all_over = 0
        #             if H > 0: all_under = 0
        # if all_over == 0 and all_under == 0 and any_radius == 1: return 1
        # return 0

    def _hash_vals(self):
        return (("norm_vec[0]", self.norm_vec[0]),
                ("norm_vec[1]", self.norm_vec[1]),
                ("norm_vec[2]", self.norm_vec[2]),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]),
                ("radius", self.radius),
                ("radius2", self.radius2),
                ("height", self.height))

disk_selector = DiskSelector

cdef class CuttingPlaneSelector(SelectorObject):
    cdef np.float64_t norm_vec[3]
    cdef np.float64_t d

    def __init__(self, dobj):
        cdef int i
        for i in range(3):
            self.norm_vec[i] = dobj._norm_vec[i]
        self.d = _ensure_code(dobj._d)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        cdef int i
        for i in range(3):
            left_edge[i] = pos[i] - 0.5*dds[i]
            right_edge[i] = pos[i] + 0.5*dds[i]
        return self.select_bbox(left_edge, right_edge)

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int i
        cdef np.float64_t height = self.d
        for i in range(3) :
            height += pos[i] * self.norm_vec[i]
        if height*height <= radius*radius : return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, j, k, n
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3]
        cdef np.float64_t gd
        arr[0] = left_edge
        arr[1] = right_edge
        all_under = 1
        all_over = 1
        # Check each corner
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    gd = self.d
                    for n in range(3):
                        gd += pos[n] * self.norm_vec[n]
                    # this allows corners and faces on the low-end to
                    # collide, while not selecting cells on the high-side
                    if i == 0 and j == 0 and k == 0 :
                        if gd <= 0: all_over = 0
                        if gd >= 0: all_under = 0
                    else :
                        if gd < 0: all_over = 0
                        if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, j, k, n
        cdef np.float64_t *arr[2]
        cdef np.float64_t pos[3]
        cdef np.float64_t gd
        arr[0] = left_edge
        arr[1] = right_edge
        all_under = 1
        all_over = 1
        # Check each corner
        for i in range(2):
            pos[0] = arr[i][0]
            for j in range(2):
                pos[1] = arr[j][1]
                for k in range(2):
                    pos[2] = arr[k][2]
                    gd = self.d
                    for n in range(3):
                        gd += pos[n] * self.norm_vec[n]
                    # this allows corners and faces on the low-end to
                    # collide, while not selecting cells on the high-side
                    if i == 0 and j == 0 and k == 0 :
                        if gd <= 0: all_over = 0
                        if gd >= 0: all_under = 0
                    else :
                        if gd < 0: all_over = 0
                        if gd > 0: all_under = 0
        if all_over == 1 or all_under == 1:
            return 0
        return 2 # a box of non-zeros volume can't be inside a plane

    def _hash_vals(self):
        return (("norm_vec[0]", self.norm_vec[0]),
                ("norm_vec[1]", self.norm_vec[1]),
                ("norm_vec[2]", self.norm_vec[2]),
                ("d", self.d))

cutting_selector = CuttingPlaneSelector

cdef class SliceSelector(SelectorObject):
    cdef int axis
    cdef np.float64_t coord
    cdef int ax, ay

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.coord = _ensure_code(dobj.coord)

        self.ax = (self.axis+1) % 3
        self.ay = (self.axis+2) % 3

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k
        cdef int total = 0
        cdef int this_level = 0
        cdef int ind[3][2]
        cdef np.uint64_t icoord
        cdef np.int32_t level = gobj.Level
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.dds)

        if level < self.min_level or level > self.max_level:
            return None
        else:
            child_mask = gobj.child_mask
            mask = np.zeros(gobj.ActiveDimensions, dtype=np.uint8)
            if level == self.max_level:
                this_level = 1
            for i in range(3):
                if i == self.axis:
                    icoord = <np.uint64_t>(
                        (self.coord - gobj.LeftEdge.d[i])/gobj.dds[i])
                    # clip coordinate to avoid seg fault below if we're
                    # exactly at a grid boundary
                    ind[i][0] = iclip(
                        icoord, 0, gobj.ActiveDimensions[i]-1)
                    ind[i][1] = ind[i][0] + 1
                else:
                    ind[i][0] = 0
                    ind[i][1] = gobj.ActiveDimensions[i]
            with nogil:
                for i in range(ind[0][0], ind[0][1]):
                    for j in range(ind[1][0], ind[1][1]):
                        for k in range(ind[2][0], ind[2][1]):
                            if this_level == 1 or child_mask[i, j, k]:
                                mask[i, j, k] = 1
                                total += 1
            if total == 0: return None
            return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if pos[self.axis] + 0.5*dds[self.axis] > self.coord \
           and pos[self.axis] - 0.5*dds[self.axis] - grid_eps <= self.coord:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef np.float64_t dist = self.periodic_difference(
            pos[self.axis], self.coord, self.axis)
        if dist*dist < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.axis] - grid_eps <= self.coord < right_edge[self.axis]:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.axis] - grid_eps <= self.coord < right_edge[self.axis]:
            return 2 # a box with non-zero volume can't be inside a plane
        return 0

    def _hash_vals(self):
        return (("axis", self.axis),
                ("coord", self.coord))

slice_selector = SliceSelector

cdef class OrthoRaySelector(SelectorObject):

    cdef np.uint8_t px_ax
    cdef np.uint8_t py_ax
    cdef np.float64_t px
    cdef np.float64_t py
    cdef int axis

    def __init__(self, dobj):
        self.axis = dobj.axis
        self.px_ax = dobj.px_ax
        self.py_ax = dobj.py_ax
        self.px = dobj.px
        self.py = dobj.py

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k
        cdef int total = 0
        cdef int this_level = 0
        cdef int ind[3][2]
        cdef np.int32_t level = gobj.Level
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)

        if level < self.min_level or level > self.max_level:
            return None
        else:
            child_mask = gobj.child_mask
            mask = np.zeros(gobj.ActiveDimensions, dtype=np.uint8)
            if level == self.max_level:
                this_level = 1
            ind[self.axis][0] = 0
            ind[self.axis][1] = gobj.ActiveDimensions[self.axis]
            ind[self.px_ax][0] = \
                <int> ((self.px - (gobj.LeftEdge).to_ndarray()[self.px_ax]) /
                       gobj.dds[self.px_ax])
            ind[self.px_ax][1] = ind[self.px_ax][0] + 1
            ind[self.py_ax][0] = \
                <int> ((self.py - (gobj.LeftEdge).to_ndarray()[self.py_ax]) /
                       gobj.dds[self.py_ax])
            ind[self.py_ax][1] = ind[self.py_ax][0] + 1

            with nogil:
                for i in range(ind[0][0], ind[0][1]):
                    for j in range(ind[1][0], ind[1][1]):
                        for k in range(ind[2][0], ind[2][1]):
                            if this_level == 1 or child_mask[i, j, k]:
                                mask[i, j, k] = 1
                                total += 1
            if total == 0: return None
            return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if self.px >= pos[self.px_ax] - 0.5*dds[self.px_ax] and \
           self.px <  pos[self.px_ax] + 0.5*dds[self.px_ax] and \
           self.py >= pos[self.py_ax] - 0.5*dds[self.py_ax] and \
           self.py <  pos[self.py_ax] + 0.5*dds[self.py_ax]:
            return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef np.float64_t dx = self.periodic_difference(
            pos[self.px_ax], self.px, self.px_ax)
        cdef np.float64_t dy = self.periodic_difference(
            pos[self.py_ax], self.py, self.py_ax)
        if dx*dx + dy*dy < radius*radius:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.px_ax] <= self.px < right_edge[self.px_ax] and \
           left_edge[self.py_ax] <= self.py < right_edge[self.py_ax] :
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if left_edge[self.px_ax] <= self.px < right_edge[self.px_ax] and \
           left_edge[self.py_ax] <= self.py < right_edge[self.py_ax] :
            return 2 # a box of non-zero volume can't be inside a ray
        return 0

    def _hash_vals(self):
        return (("px_ax", self.px_ax),
                ("py_ax", self.py_ax),
                ("px", self.px),
                ("py", self.py),
                ("axis", self.axis))

ortho_ray_selector = OrthoRaySelector

cdef struct IntegrationAccumulator:
    np.float64_t *t
    np.float64_t *dt
    np.uint8_t *child_mask
    int hits

cdef void dt_sampler(
             VolumeContainer *vc,
             np.float64_t v_pos[3],
             np.float64_t v_dir[3],
             np.float64_t enter_t,
             np.float64_t exit_t,
             int index[3],
             void *data) nogil:
    cdef IntegrationAccumulator *am = <IntegrationAccumulator *> data
    cdef int di = (index[0]*vc.dims[1]+index[1])*vc.dims[2]+index[2]
    if am.child_mask[di] == 0 or enter_t == exit_t:
        return
    am.hits += 1
    am.t[di] = enter_t
    am.dt[di] = (exit_t - enter_t)

cdef class RaySelector(SelectorObject):

    cdef np.float64_t p1[3]
    cdef np.float64_t p2[3]
    cdef np.float64_t vec[3]

    def __init__(self, dobj):
        cdef int i
        _ensure_code(dobj.start_point)
        _ensure_code(dobj.end_point)
        for i in range(3):
            self.vec[i] = dobj.vec[i]
            self.p1[i] = dobj.start_point[i]
            self.p2[i] = dobj.end_point[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i
        cdef int total = 0
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef VolumeContainer vc
        mask = np.zeros(gobj.ActiveDimensions, dtype='uint8')
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i, j, k] >= 0:
                        mask[i, j, k] = 1
                        total += 1
        free(ia)
        if total == 0: return None
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_dt(self, gobj):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.float64_t, ndim=1] tr, dtr
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k, ni
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef VolumeContainer vc
        t = np.zeros(gobj.ActiveDimensions, dtype="float64")
        dt = np.zeros(gobj.ActiveDimensions, dtype="float64") - 1
        child_mask = gobj.child_mask
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        _ensure_code(gobj.LeftEdge)
        _ensure_code(gobj.RightEdge)
        _ensure_code(gobj.dds)
        for i in range(3):
            vc.left_edge[i] = gobj.LeftEdge[i]
            vc.right_edge[i] = gobj.RightEdge[i]
            vc.dds[i] = gobj.dds[i]
            vc.idds[i] = 1.0/gobj.dds[i]
            vc.dims[i] = dt.shape[i]
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        tr = np.zeros(ia.hits, dtype="float64")
        dtr = np.zeros(ia.hits, dtype="float64")
        ni = 0
        for i in range(dt.shape[0]):
            for j in range(dt.shape[1]):
                for k in range(dt.shape[2]):
                    if dt[i, j, k] >= 0:
                        tr[ni] = t[i, j, k]
                        dtr[ni] = dt[i, j, k]
                        ni += 1
        if not (ni == ia.hits):
            print(ni, ia.hits)
        free(ia)
        return dtr, tr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def get_dt_mesh(self, mesh, nz, int offset):
        cdef np.ndarray[np.float64_t, ndim=3] t, dt
        cdef np.ndarray[np.float64_t, ndim=1] tr, dtr
        cdef np.ndarray[np.uint8_t, ndim=3, cast=True] child_mask
        cdef int i, j, k, ni
        cdef np.float64_t LE[3]
        cdef np.float64_t RE[3]
        cdef np.float64_t pos
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef np.ndarray[np.float64_t, ndim=2] coords
        cdef np.ndarray[np.int64_t, ndim=2] indices
        indices = mesh.connectivity_indices
        coords = _ensure_code(mesh.connectivity_coords)
        cdef int nc = indices.shape[0]
        cdef int nv = indices.shape[1]
        if nv != 8:
            raise NotImplementedError
        cdef VolumeContainer vc
        child_mask = np.ones((1,1,1), dtype="uint8")
        t = np.zeros((1,1,1), dtype="float64")
        dt = np.zeros((1,1,1), dtype="float64") - 1
        tr = np.zeros(nz, dtype="float64")
        dtr = np.zeros(nz, dtype="float64")
        ia.t = <np.float64_t *> t.data
        ia.dt = <np.float64_t *> dt.data
        ia.child_mask = <np.uint8_t *> child_mask.data
        ia.hits = 0
        ni = 0
        for i in range(nc):
            for j in range(3):
                LE[j] = 1e60
                RE[j] = -1e60
            for j in range(nv):
                for k in range(3):
                    pos = coords[indices[i, j] - offset, k]
                    LE[k] = fmin(pos, LE[k])
                    RE[k] = fmax(pos, RE[k])
            for j in range(3):
                vc.left_edge[j] = LE[j]
                vc.right_edge[j] = RE[j]
                vc.dds[j] = RE[j] - LE[j]
                vc.idds[j] = 1.0/vc.dds[j]
                vc.dims[j] = 1
            t[0,0,0] = dt[0,0,0] = -1
            walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
            if dt[0,0,0] >= 0:
                tr[ni] = t[0,0,0]
                dtr[ni] = dt[0,0,0]
                ni += 1
        free(ia)
        return dtr, tr

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # two 0-volume constructs don't intersect
        return 0
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:

        cdef int i
        cdef np.float64_t length = norm(self.vec)
        cdef np.float64_t r[3]
        for i in range(3):
            r[i] = pos[i] - self.p1[i]
        # the projected position of the sphere along the ray
        cdef np.float64_t l = dot(r, self.vec) / length
        # the square of the impact parameter
        cdef np.float64_t b_sqr = dot(r, r) - l*l

        # only accept spheres with radii larger than the impact parameter and
        # with a projected position along the ray no more than a radius away
        # from the ray
        if -radius < l and l < (length+radius) and b_sqr < radius*radius:
            return 1

        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i, rv
        cdef VolumeContainer vc
        cdef IntegrationAccumulator *ia
        ia = <IntegrationAccumulator *> malloc(sizeof(IntegrationAccumulator))
        cdef np.float64_t dt[1]
        cdef np.float64_t t[1]
        cdef np.uint8_t cm[1]
        for i in range(3):
            vc.left_edge[i] = left_edge[i]
            vc.right_edge[i] = right_edge[i]
            vc.dds[i] = right_edge[i] - left_edge[i]
            vc.idds[i] = 1.0/vc.dds[i]
            vc.dims[i] = 1
        t[0] = dt[0] = 0.0
        cm[0] = 1
        ia.t = t
        ia.dt = dt
        ia.child_mask = cm
        ia.hits = 0
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> ia)
        rv = 0
        if ia.hits > 0:
            rv = 1
        free(ia)
        return rv

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int i
        cdef np.uint8_t cm = 1
        cdef VolumeContainer vc
        cdef IntegrationAccumulator ia
        cdef np.float64_t dt, t
        for i in range(3):
            vc.left_edge[i] = left_edge[i]
            vc.right_edge[i] = right_edge[i]
            vc.dds[i] = right_edge[i] - left_edge[i]
            vc.idds[i] = 1.0/vc.dds[i]
            vc.dims[i] = 1
        t = dt = 0.0
        ia.t = &t
        ia.dt = &dt
        ia.child_mask = &cm
        ia.hits = 0
        walk_volume(&vc, self.p1, self.vec, dt_sampler, <void*> &ia)
        if ia.hits > 0:
            return 2 # a box of non-zero volume cannot be inside a ray
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3],
                               np.float64_t dds[3]) nogil:
        # This is terribly inefficient for Octrees.  For grids, it will never
        # get called.
        cdef int i
        cdef np.float64_t left_edge[3]
        cdef np.float64_t right_edge[3]
        for i in range(3):
            left_edge[i] = pos[i] - dds[i]/2.0
            right_edge[i] = pos[i] + dds[i]/2.0
        return self.select_bbox(left_edge, right_edge)

    def _hash_vals(self):
        return (("p1[0]", self.p1[0]),
                ("p1[1]", self.p1[1]),
                ("p1[2]", self.p1[2]),
                ("p2[0]", self.p2[0]),
                ("p2[1]", self.p2[1]),
                ("p2[2]", self.p2[2]),
                ("vec[0]", self.vec[0]),
                ("vec[1]", self.vec[1]),
                ("vec[2]", self.vec[2]))

ray_selector = RaySelector

cdef class DataCollectionSelector(SelectorObject):
    cdef object obj_ids
    cdef np.int64_t nids

    def __init__(self, dobj):
        self.obj_ids = dobj._obj_ids
        self.nids = self.obj_ids.shape[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int n
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        cdef np.ndarray[np.int64_t, ndim=1] oids = self.obj_ids
        with nogil:
            for n in range(self.nids):
                gridi[oids[n]] = 1
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        cdef np.ndarray[np.uint8_t, ndim=3] mask
        mask = np.ones(gobj.ActiveDimensions, dtype='uint8')
        return mask.astype("bool")

    def _hash_vals(self):
        return (hash(self.obj_ids.tostring()), self.nids)

data_collection_selector = DataCollectionSelector

cdef class EllipsoidSelector(SelectorObject):
    cdef np.float64_t vec[3][3]
    cdef np.float64_t mag[3]
    cdef np.float64_t center[3]

    def __init__(self, dobj):
        cdef int i
        _ensure_code(dobj.center)
        _ensure_code(dobj._e0)
        _ensure_code(dobj._e1)
        _ensure_code(dobj._e2)
        _ensure_code(dobj._A)
        _ensure_code(dobj._B)
        _ensure_code(dobj._C)
        for i in range(3):
            self.center[i] = dobj.center[i]
            self.vec[0][i] = dobj._e0[i]
            self.vec[1][i] = dobj._e1[i]
            self.vec[2][i] = dobj._e2[i]
        self.mag[0] = dobj._A
        self.mag[1] = dobj._B
        self.mag[2] = dobj._C

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return self.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef np.float64_t dot_evec[3]
        cdef np.float64_t dist
        cdef int i, j
        dot_evec[0] = dot_evec[1] = dot_evec[2] = 0
        # Calculate the rotated dot product
        for i in range(3): # axis
            dist = self.periodic_difference(pos[i], self.center[i], i)
            for j in range(3):
                dot_evec[j] += dist * self.vec[j][i]
        dist = 0.0
        for i in range(3):
            dist += (dot_evec[i] * dot_evec[i])/(self.mag[i] * self.mag[i])
        if dist <= 1.0: return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        # this is the sphere selection
        cdef int i
        cdef np.float64_t dist, dist2_max, dist2 = 0
        for i in range(3):
            dist = self.periodic_difference(pos[i], self.center[i], i)
            dist2 += dist * dist
        dist2_max = (self.mag[0] + radius) * (self.mag[0] + radius)
        if dist2 <= dist2_max:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # This is the sphere selection
        cdef int i
        cdef np.float64_t box_center, relcenter, closest, dist, edge, dist_max
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2]:
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        dist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            dist += closest * closest
        dist_max = self.mag[0] * self.mag[0]
        if dist <= dist_max:
            return 1
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # This is the sphere selection
        cdef int i
        cdef np.float64_t box_center, relcenter, closest, farthest, cdist, fdist, edge
        if left_edge[0] <= self.center[0] <= right_edge[0] and \
           left_edge[1] <= self.center[1] <= right_edge[1] and \
           left_edge[2] <= self.center[2] <= right_edge[2]:
            fdist = 0
            for i in range(3):
                edge = right_edge[i] - left_edge[i]
                box_center = (right_edge[i] + left_edge[i])/2.0
                relcenter = self.periodic_difference(
                    box_center, self.center[i], i)
                farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
                fdist += farthest*farthest
                if fdist >= self.mag[0]**2: return 2
            return 1
        # http://www.gamedev.net/topic/335465-is-this-the-simplest-sphere-aabb-collision-test/
        cdist = 0
        fdist = 0
        for i in range(3):
            box_center = (right_edge[i] + left_edge[i])/2.0
            relcenter = self.periodic_difference(box_center, self.center[i], i)
            edge = right_edge[i] - left_edge[i]
            closest = relcenter - fclip(relcenter, -edge/2.0, edge/2.0)
            farthest = relcenter + fclip(relcenter, -edge/2.0, edge/2.0)
            cdist += closest * closest
            fdist += farthest * farthest
            if cdist > self.mag[0]**2: return 0
        if fdist < self.mag[0]**2:
            return 1
        else:
            return 2

    def _hash_vals(self):
        return (("vec[0][0]", self.vec[0][0]),
                ("vec[0][1]", self.vec[0][1]),
                ("vec[0][2]", self.vec[0][2]),
                ("vec[1][0]", self.vec[1][0]),
                ("vec[1][1]", self.vec[1][1]),
                ("vec[1][2]", self.vec[1][2]),
                ("vec[2][0]", self.vec[2][0]),
                ("vec[2][1]", self.vec[2][1]),
                ("vec[2][2]", self.vec[2][2]),
                ("mag[0]", self.mag[0]),
                ("mag[1]", self.mag[1]),
                ("mag[2]", self.mag[2]),
                ("center[0]", self.center[0]),
                ("center[1]", self.center[1]),
                ("center[2]", self.center[2]))

ellipsoid_selector = EllipsoidSelector

cdef class GridSelector(SelectorObject):
    cdef object ind

    def __init__(self, dobj):
        self.ind = dobj.id - dobj._id_offset

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.zeros(ng, dtype='uint8')
        gridi[self.ind] = 1
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_mask(self, gobj):
        return np.ones(gobj.ActiveDimensions, dtype='bool')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        # we apparently don't check if the point actually lies in the grid..
        return 1

    def _hash_vals(self):
        return (self.ind,)

grid_selector = GridSelector

cdef class OctreeSubsetSelector(SelectorObject):

    def __init__(self, dobj):
        self.base_selector = dobj.base_selector
        self.min_level = self.base_selector.min_level
        self.max_level = self.base_selector.max_level
        self.domain_id = dobj.domain_id
        self.overlap_cells = getattr(dobj.oct_handler, 'overlap_cells', 1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return self.base_selector.select_point(pos)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # return 1
        return self.base_selector.select_bbox(left_edge, right_edge)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        cdef int res
        res = self.base_selector.select_grid(left_edge, right_edge, level, o)
        if self.domain_id == -1:
            return res
        elif res == 1 and o != NULL and o.domain != self.domain_id:
            return -1
        return res

    def _hash_vals(self):
        return (hash(self.base_selector), self.domain_id)

octree_subset_selector = OctreeSubsetSelector

cdef class IndexedOctreeSubsetSelector(SelectorObject):
    # This is a numpy array, which will be a bool of ndim 1
    cdef np.uint64_t min_ind
    cdef np.uint64_t max_ind
    cdef public SelectorObject base_selector
    cdef int filter_bbox
    cdef np.float64_t DLE[3]
    cdef np.float64_t DRE[3]

    def __init__(self, dobj):
        self.min_ind = dobj.min_ind
        self.max_ind = dobj.max_ind
        self.base_selector = dobj.base_selector
        self.min_level = self.base_selector.min_level
        self.max_level = self.base_selector.max_level
        self.filter_bbox = 0
        if getattr(dobj.ds, "filter_bbox", False):
            self.filter_bbox = 1
        for i in range(3):
            self.DLE[i] = dobj.ds.domain_left_edge[i]
            self.DRE[i] = dobj.ds.domain_right_edge[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        raise RuntimeError

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int i
        if self.filter_bbox == 0:
            return 1
        for i in range(3):
            if pos[i] < self.DLE[i] or pos[i] > self.DRE[i]:
                return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return self.base_selector.select_bbox(left_edge, right_edge)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        # Because visitors now use select_grid, we should be explicitly
        # checking this.
        return self.base_selector.select_grid(left_edge, right_edge, level, o)

    def _hash_vals(self):
        return (hash(self.base_selector), self.min_ind, self.max_ind)

indexed_octree_subset_selector = IndexedOctreeSubsetSelector

cdef class AlwaysSelector(SelectorObject):

    def __init__(self, dobj):
        self.overlap_cells = 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        cdef int ng = left_edges.shape[0]
        cdef np.ndarray[np.uint8_t, ndim=1] gridi = np.ones(ng, dtype='uint8')
        return gridi.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        return 1

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        return 1

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        return 1

    def _hash_vals(self):
        return ("always", 1,)

always_selector = AlwaysSelector

cdef class ComposeSelector(SelectorObject):
    cdef SelectorObject selector1
    cdef SelectorObject selector2

    def __init__(self, dobj, selector1, selector2):
        self.selector1 = selector1
        self.selector2 = selector2
        self.min_level = max(selector1.min_level, selector2.min_level)
        self.max_level = min(selector1.max_level, selector2.max_level)

    def select_grids(self,
                     np.ndarray[np.float64_t, ndim=2] left_edges,
                     np.ndarray[np.float64_t, ndim=2] right_edges,
                     np.ndarray[np.int32_t, ndim=2] levels):
        return np.logical_or(
                    self.selector1.select_grids(left_edges, right_edges, levels),
                    self.selector2.select_grids(left_edges, right_edges, levels))

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        if self.selector1.select_cell(pos, dds) and \
                self.selector2.select_cell(pos, dds):
            return 1
        else:
            return 0

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        if self.selector1.select_grid(left_edge, right_edge, level, o) or \
                self.selector2.select_grid(left_edge, right_edge, level, o):
            return 1
        else:
            return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        if self.selector1.select_point(pos) and \
                self.selector2.select_point(pos):
            return 1
        else:
            return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        if self.selector1.select_sphere(pos, radius) and \
                self.selector2.select_sphere(pos, radius):
            return 1
        else:
            return 0

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        if self.selector1.select_bbox(left_edge, right_edge) and \
                self.selector2.select_bbox(left_edge, right_edge):
            return 1
        else:
            return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.selector1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.selector2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 0: return 0
        return max(rv1, rv2)

    def _hash_vals(self):
        return (hash(self.selector1), hash(self.selector2))

compose_selector = ComposeSelector

cdef class HaloParticlesSelector(SelectorObject):
    cdef public object base_source
    cdef SelectorObject base_selector
    cdef object pind
    cdef public np.int64_t halo_id
    def __init__(self, dobj):
        self.base_source = dobj.base_source
        self.base_selector = self.base_source.selector
        self.pind = dobj.particle_indices

    def _hash_vals(self):
        return ("halo_particles", self.halo_id)

halo_particles_selector = HaloParticlesSelector

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def points_in_cells(
        np.float64_t[:] cx,
        np.float64_t[:] cy,
        np.float64_t[:] cz,
        np.float64_t[:] dx,
        np.float64_t[:] dy,
        np.float64_t[:] dz,
        np.float64_t[:] px,
        np.float64_t[:] py,
        np.float64_t[:] pz):
    # Take a list of cells and particles and calculate which particles
    # are enclosed within one of the cells.  This is used for querying
    # particle fields on clump/contour objects.
    # We use brute force since the cells are a relatively unordered collection.

    cdef int p, c, n_p, n_c
    cdef np.ndarray[np.uint8_t, ndim=1, cast=True] mask

    n_p = px.size
    n_c = cx.size
    mask = np.zeros(n_p, dtype="bool")

    for p in range(n_p):
        for c in range(n_c):
            if (fabs(px[p] - cx[c]) <= 0.5 * dx[c] and
                fabs(py[p] - cy[c]) <= 0.5 * dy[c] and
                fabs(pz[p] - cz[c]) <= 0.5 * dz[c]):
                mask[p] = True
                break

    return mask

cdef class BooleanSelector(SelectorObject):

    def __init__(self, dobj):
        # Note that this has a different API than the other selector objects,
        # so will not work as a traditional data selector.
        if not hasattr(dobj.dobj1, "selector"):
            self.sel1 = dobj.dobj1
        else:
            self.sel1 = dobj.dobj1.selector
        if not hasattr(dobj.dobj2, "selector"):
            self.sel2 = dobj.dobj2
        else:
            self.sel2 = dobj.dobj2.selector

cdef class BooleanANDSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox(left_edge, right_edge)
        if rv2 == 0: return 0
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 0: return 0
        return max(rv1, rv2)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        cdef int rv1 = self.sel1.select_grid(left_edge, right_edge, level, o)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_grid(left_edge, right_edge, level, o)
        if rv2 == 0: return 0
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 0: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 0: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 0: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("and",) +
                self.sel2._hash_vals())

cdef class BooleanORSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox(left_edge, right_edge)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_bbox(left_edge, right_edge)
        if rv2 == 1: return 1
        return 0

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 1: return 1
        return max(rv1, rv2)

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        cdef int rv1 = self.sel1.select_grid(left_edge, right_edge, level, o)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_grid(left_edge, right_edge, level, o)
        if rv2 == 1: return 1
        if (rv1 == 2) or (rv2 == 2): return 2
        return 0

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 1: return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 1: return 1
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 1: return 1
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 1: return 1
        return 0

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("or",) +
                self.sel2._hash_vals())

cdef class BooleanNOTSelector(BooleanSelector):
    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 1
        elif rv1 == 1: return 0
        return 2

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 1
        return 0

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 1
        return 0

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 1
        return 0

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("not",))

cdef class BooleanXORSelector(BooleanSelector):

    cdef int select_bbox(self, np.float64_t left_edge[3],
                               np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return 1

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        # Return 2 in cases where one or both selectors partially overlap since
        # part of the bounding box could satisfy the condition unless the
        # selectors are identical.
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv1 == rv2:
            if rv1 == 2: 
                # If not identical, part of the bbox will be touched by one
                # selector and not the other.
                # if self.sel1 == self.sel2: return 0  # requires gil
                return 2
            return 0
        if rv1 == 0: return rv2
        if rv2 == 0: return rv1
        return 2  # part of bbox only touched by selector fully covering bbox

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return 1

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv1 == rv2: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        cdef int rv2 = self.sel2.select_point(pos)
        if rv1 == rv2: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv1 == rv2: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("xor",) +
                self.sel2._hash_vals())

cdef class BooleanNEGSelector(BooleanSelector):

    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        # We always return True here, because we don't have a "fully included"
        # check anywhere else.
        return self.sel1.select_bbox(left_edge, right_edge)

    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int rv1 = self.sel1.select_bbox_edge(left_edge, right_edge)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_bbox_edge(left_edge, right_edge)
        if rv2 == 1:
            return 0
        elif rv2 == 0:
            return rv1
        # If sel2 is partial, then sel1 - sel2 will be partial as long
        # as sel1 != sel2
        # if self.sel1 == self.sel2: return 0  # requires gil
        return 2  

    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        return self.sel1.select_grid(left_edge, right_edge, level, o)

    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        cdef int rv1 = self.sel1.select_cell(pos, dds)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_cell(pos, dds)
        if rv2 == 1: return 0
        return 1

    cdef int select_point(self, np.float64_t pos[3]) nogil:
        cdef int rv1 = self.sel1.select_point(pos)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_point(pos)
        if rv2 == 1: return 0
        return 1

    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        cdef int rv1 = self.sel1.select_sphere(pos, radius)
        if rv1 == 0: return 0
        cdef int rv2 = self.sel2.select_sphere(pos, radius)
        if rv2 == 1: return 0
        return 1

    def _hash_vals(self):
        return (self.sel1._hash_vals() +
                ("neg",) +
                self.sel2._hash_vals())

cdef class ChainedBooleanSelector(SelectorObject):
    cdef int n_obj
    cdef np.ndarray selectors
    def __init__(self, dobj):
        # These are data objects, not selectors
        self.n_obj = len(dobj.data_objects)
        self.selectors = np.empty(self.n_obj, dtype="object")
        for i in range(self.n_obj):
            self.selectors[i] = dobj.data_objects[i].selector

cdef class ChainedBooleanANDSelector(ChainedBooleanSelector):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_bbox(
                        left_edge, right_edge) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                              np.float64_t right_edge[3]) nogil:
        cdef int selected = 1
        cdef int ret
        with gil:
            for i in range(self.n_obj):
                ret = (<SelectorObject>self.selectors[i]).select_bbox_edge(
                    left_edge, right_edge)
                if ret == 0:
                    return 0
                elif ret == 2:
                    selected = 2
        return selected

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_grid(
                        left_edge, right_edge, level, o) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_cell(
                        pos, dds) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_point(pos) == 0:
                    return 0
        return 1

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_sphere(
                        pos, radius) == 0:
                    return 0
        return 1

    def _hash_vals(self):
        v = ("chained_and",)
        for s in self.selectors:
            v += s._hash_vals()
        return v

intersection_selector = ChainedBooleanANDSelector

cdef class ChainedBooleanORSelector(ChainedBooleanSelector):
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_bbox(
                        left_edge, right_edge) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_bbox_edge(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3]) nogil:
        cdef int selected = 0
        cdef int ret
        with gil:
            for i in range(self.n_obj):
                ret = (<SelectorObject>self.selectors[i]).select_bbox_edge(
                    left_edge, right_edge)
                if ret == 1:
                    return 1
                elif ret == 2:
                    selected = 2
        return selected

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_grid(self, np.float64_t left_edge[3],
                         np.float64_t right_edge[3], np.int32_t level,
                         Oct *o = NULL) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_grid(
                        left_edge, right_edge, level, o) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_cell(self, np.float64_t pos[3], np.float64_t dds[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_cell(
                        pos, dds) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_point(self, np.float64_t pos[3]) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_point(pos) == 1:
                    return 1
        return 0

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int select_sphere(self, np.float64_t pos[3], np.float64_t radius) nogil:
        with gil:
            for i in range(self.n_obj):
                if (<SelectorObject>self.selectors[i]).select_sphere(
                        pos, radius) == 1:
                    return 1
        return 0

    def _hash_vals(self):
        v = ("chained_or",)
        for s in self.selectors:
            v += s._hash_vals()
        return v

union_selector = ChainedBooleanORSelector
