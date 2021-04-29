# distutils: sources = yt/utilities/lib/tsearch.c
# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Oct container




"""


cimport cython
cimport numpy as np

import numpy as np

from libc.math cimport ceil, floor
from selection_routines cimport AlwaysSelector, SelectorObject

from yt.geometry.oct_visitors cimport (
    NeighbourCellIndexVisitor,
    NeighbourCellVisitor,
    OctPadded,
    StoreIndex,
)

ORDER_MAX = 20
_ORDER_MAX = ORDER_MAX

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

# Here is the strategy for RAMSES containers:
#   * Read each domain individually, creating *all* octs found in that domain
#     file, even if they reside on other CPUs.
#   * Only allocate octs that reside on >= domain
#   * For all octs, insert into tree, which may require traversing existing
#     octs
#   * Note that this does not allow one component of an ObjectPool (an
#     AllocationContainer) to exactly be a chunk, but it is close.  For IO
#     chunking, we can theoretically examine those octs that live inside a
#     given allocator.

cdef class OctreeContainer:

    def __init__(self, oct_domain_dimensions, domain_left_edge,
                 domain_right_edge, partial_coverage = 0,
                 over_refine = 1):
        # This will just initialize the root mesh octs
        self.oref = over_refine
        self.partial_coverage = partial_coverage
        cdef int i, j, k, p
        for i in range(3):
            self.nn[i] = oct_domain_dimensions[i]
        self.num_domains = 0
        self.level_offset = 0
        self.domains = OctObjectPool()
        p = 0
        self.nocts = 0 # Increment when initialized
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid
        self._initialize_root_mesh()
        self.fill_style = "o"

    def _initialize_root_mesh(self):
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = NULL

    @property
    def oct_arrays(self):
        return self.domains.to_arrays()

    @classmethod
    def load_octree(cls, header):
        cdef np.ndarray[np.uint8_t, ndim=1] ref_mask
        ref_mask = header['octree']
        cdef OctreeContainer obj = cls(header['dims'], header['left_edge'],
                header['right_edge'], over_refine = header['over_refine'],
                partial_coverage = header['partial_coverage'])
        # NOTE: We do not allow domain/file indices to be specified.
        cdef SelectorObject selector = AlwaysSelector(None)
        cdef oct_visitors.LoadOctree visitor
        visitor = oct_visitors.LoadOctree(obj, -1)
        cdef int i, j, k, n
        visitor.global_index = -1
        visitor.level = 0
        visitor.oref = 0
        visitor.nz = 1
        assert(ref_mask.shape[0] / float(visitor.nz) ==
            <int>(ref_mask.shape[0]/float(visitor.nz)))
        obj.allocate_domains([ref_mask.shape[0] / visitor.nz])
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (obj.DRE[i] - obj.DLE[i]) / obj.nn[i]
        # Pos is the center of the octs
        cdef OctAllocationContainer *cur = obj.domains.get_cont(0)
        cdef Oct *o
        cdef np.uint64_t nfinest = 0
        visitor.ref_mask = ref_mask
        visitor.octs = cur.my_objs
        visitor.nocts = &cur.n_assigned
        visitor.nfinest = &nfinest
        pos[0] = obj.DLE[0] + dds[0]/2.0
        for i in range(obj.nn[0]):
            pos[1] = obj.DLE[1] + dds[1]/2.0
            for j in range(obj.nn[1]):
                pos[2] = obj.DLE[2] + dds[2]/2.0
                for k in range(obj.nn[2]):
                    if obj.root_mesh[i][j][k] != NULL:
                        raise RuntimeError
                    o = &cur.my_objs[cur.n_assigned]
                    o.domain_ind = o.file_ind = 0
                    o.domain = 1
                    obj.root_mesh[i][j][k] = o
                    cur.n_assigned += 1
                    visitor.pos[0] = i
                    visitor.pos[1] = j
                    visitor.pos[2] = k
                    # Always visit covered
                    selector.recursively_visit_octs(
                        obj.root_mesh[i][j][k],
                        pos, dds, 0, visitor, 1)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]
        obj.nocts = cur.n_assigned
        if obj.nocts * visitor.nz != ref_mask.size:
            raise KeyError(ref_mask.size, obj.nocts, obj.oref,
                obj.partial_coverage, visitor.nz)
        return obj

    def __dealloc__(self):
        if self.root_mesh == NULL: return
        for i in range(self.nn[0]):
            if self.root_mesh[i] == NULL: continue
            for j in range(self.nn[1]):
                if self.root_mesh[i][j] == NULL: continue
                free(self.root_mesh[i][j])
            if self.root_mesh[i] == NULL: continue
            free(self.root_mesh[i])
        free(self.root_mesh)

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        OctVisitor visitor, int vc = -1,
                        np.int64_t *indices = NULL):
        cdef int i, j, k, n
        if vc == -1:
            vc = self.partial_coverage
        visitor.global_index = -1
        visitor.level = 0
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i]) / self.nn[i]
        # Pos is the center of the octs
        pos[0] = self.DLE[0] + dds[0]/2.0
        for i in range(self.nn[0]):
            pos[1] = self.DLE[1] + dds[1]/2.0
            for j in range(self.nn[1]):
                pos[2] = self.DLE[2] + dds[2]/2.0
                for k in range(self.nn[2]):
                    if self.root_mesh[i][j][k] == NULL:
                        raise RuntimeError
                    visitor.pos[0] = i
                    visitor.pos[1] = j
                    visitor.pos[2] = k
                    selector.recursively_visit_octs(
                        self.root_mesh[i][j][k],
                        pos, dds, 0, visitor, vc)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0

    cdef int get_root(self, int ind[3], Oct **o) nogil:
        cdef int i
        for i in range(3):
            if ind[i] < 0 or ind[i] >= self.nn[i]:
                o[0] = NULL
                return 1
        o[0] = self.root_mesh[ind[0]][ind[1]][ind[2]]
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct *get(self, np.float64_t ppos[3], OctInfo *oinfo = NULL,
                  int max_level = 99) nogil:
        #Given a floating point position, retrieve the most
        #refined oct at that time
        cdef int ind32[3]
        cdef np.int64_t ipos[3]
        cdef np.float64_t dds[3]
        cdef np.float64_t cp[3]
        cdef Oct *cur
        cdef Oct *next
        cdef int i
        cur = next = NULL
        cdef np.int64_t ind[3]
        cdef np.int64_t level = -1
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
            ind[i] = <np.int64_t> (floor((ppos[i] - self.DLE[i])/dds[i]))
            cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
            ipos[i] = 0 # We add this to ind later, so it should be zero.
            ind32[i] = ind[i]
        self.get_root(ind32, &next)
        # We want to stop recursing when there's nowhere else to go
        while next != NULL and level < max_level:
            level += 1
            for i in range(3):
                ipos[i] = (ipos[i] << 1) + ind[i]
            cur = next
            for i in range(3):
                dds[i] = dds[i] / 2.0
                if cp[i] > ppos[i]:
                    ind[i] = 0
                    cp[i] -= dds[i] / 2.0
                else:
                    ind[i] = 1
                    cp[i] += dds[i]/2.0
            if cur.children != NULL:
                next = cur.children[cind(ind[0],ind[1],ind[2])]
            else:
                next = NULL
        if oinfo == NULL: return cur
        cdef int ncells = (1 << self.oref)
        cdef np.float64_t factor = 1.0 / (1 << (self.oref-1))
        if self.oref == 0: factor = 2.0
        for i in range(3):
            # We don't normally need to change dds[i] as it has been halved
            # from the oct width, thus making it already the cell width.
            # But, since not everything has the cell width equal to have the
            # width of the oct, we need to apply "factor".
            oinfo.dds[i] = dds[i] * factor # Cell width
            oinfo.ipos[i] = ipos[i]
            oinfo.left_edge[i] = oinfo.ipos[i] * (oinfo.dds[i] * ncells) + self.DLE[i]
        oinfo.level = level
        return cur

    def locate_positions(self, np.float64_t[:,:] positions):
        """
        This routine, meant to be called by other internal routines, returns a
        list of oct IDs and a dictionary of Oct info for all the positions
        supplied.  Positions must be in code_length.
        """
        cdef np.float64_t factor = (1 << self.oref)
        cdef dict all_octs = {}
        cdef OctInfo oi
        cdef Oct* o = NULL
        cdef np.float64_t pos[3]
        cdef np.ndarray[np.uint8_t, ndim=1] recorded
        cdef np.ndarray[np.int64_t, ndim=1] oct_id
        oct_id = np.ones(positions.shape[0], dtype="int64") * -1
        recorded = np.zeros(self.nocts, dtype="uint8")
        cdef np.int64_t i, j, k
        for i in range(positions.shape[0]):
            for j in range(3):
                pos[j] = positions[i,j]
            o = self.get(pos, &oi)
            if o == NULL:
                raise RuntimeError
            if recorded[o.domain_ind] == 0:
                left_edge = np.asarray(<np.float64_t[:3]>oi.left_edge).copy()
                dds = np.asarray(<np.float64_t[:3]>oi.dds).copy()
                right_edge = left_edge + dds*factor
                all_octs[o.domain_ind] = dict(
                    left_edge = left_edge,
                    right_edge = right_edge,
                    level = oi.level
                )
                recorded[o.domain_ind] = 1
            oct_id[i] = o.domain_ind
        return oct_id, all_octs

    def domain_identify(self, SelectorObject selector):
        cdef np.ndarray[np.uint8_t, ndim=1] domain_mask
        domain_mask = np.zeros(self.num_domains, dtype="uint8")
        cdef oct_visitors.IdentifyOcts visitor
        visitor = oct_visitors.IdentifyOcts(self)
        visitor.domain_mask = domain_mask
        self.visit_all_octs(selector, visitor)
        cdef int i
        domain_ids = []
        for i in range(self.num_domains):
            if domain_mask[i] == 1:
                domain_ids.append(i+1)
        return domain_ids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Oct** neighbors(self, OctInfo *oi, np.int64_t *nneighbors, Oct *o,
                         bint periodicity[3]):
        # We are going to do a brute-force search here.
        # This is not the most efficient -- in fact, it's relatively bad.  But
        # we will attempt to improve it in a future iteration, where we will
        # grow a stack of parent Octs.
        # Note that in the first iteration, we will just find the up-to-27
        # neighbors, including the main oct.
        cdef np.int64_t i, j, k, n, level, ii, dlevel
        cdef int ind[3]
        cdef OctList *olist
        cdef OctList *my_list
        my_list = olist = NULL
        cdef Oct *cand
        cdef np.int64_t npos[3]
        cdef np.int64_t ndim[3]
        # Now we get our boundaries for this level, so that we can wrap around
        # if need be.
        # ndim is the oct dimensions of the level, not the cell dimensions.
        for i in range(3):
            ndim[i] = <np.int64_t> ((self.DRE[i] - self.DLE[i]) / oi.dds[i])
            # Here we adjust for oi.dds meaning *cell* width.
            ndim[i] = (ndim[i] >> self.oref)
        my_list = olist = OctList_append(NULL, o)
        for i in range(3):
            npos[0] = (oi.ipos[0] + (1 - i))
            if not periodicity[0] and not \
               (0 <= npos[0] < ndim[0]):
                continue
            elif npos[0] < 0: npos[0] += ndim[0]
            elif npos[0] >= ndim[0]: npos[0] -= ndim[0]
            for j in range(3):
                npos[1] = (oi.ipos[1] + (1 - j))
                if not periodicity[1] and not \
                   (0 <= npos[1] < ndim[1]):
                    continue
                elif npos[1] < 0: npos[1] += ndim[1]
                elif npos[1] >= ndim[1]: npos[1] -= ndim[1]
                for k in range(3):
                    npos[2] = (oi.ipos[2] + (1 - k))
                    if not periodicity[2] and not \
                       (0 <= npos[2] < ndim[2]):
                        continue
                    if npos[2] < 0: npos[2] += ndim[2]
                    if npos[2] >= ndim[2]: npos[2] -= ndim[2]
                    # Now we have our npos, which we just need to find.
                    # Level 0 gets bootstrapped
                    for n in range(3):
                        ind[n] = ((npos[n] >> (oi.level)) & 1)
                    cand = NULL
                    self.get_root(ind, &cand)
                    # We should not get a NULL if we handle periodicity
                    # correctly, but we might.
                    if cand == NULL: continue
                    for level in range(1, oi.level+1):
                        dlevel = oi.level - level
                        if cand.children == NULL: break
                        for n in range(3):
                            ind[n] = (npos[n] >> dlevel) & 1
                        ii = cind(ind[0],ind[1],ind[2])
                        if cand.children[ii] == NULL: break
                        cand = cand.children[ii]
                    if cand.children != NULL:
                        olist = OctList_subneighbor_find(
                            olist, cand, i, j, k)
                    else:
                        olist = OctList_append(olist, cand)
        olist = my_list
        cdef int noct = OctList_count(olist)
        cdef Oct **neighbors
        neighbors = <Oct **> malloc(sizeof(Oct*)*noct)
        for i in range(noct):
            neighbors[i] = olist.o
            olist = olist.next
        OctList_delete(my_list)
        nneighbors[0] = noct
        return neighbors

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mask(self, SelectorObject selector, np.int64_t num_cells = -1,
             int domain_id = -1):
        if num_cells == -1:
            num_cells = selector.count_octs(self, domain_id)
        cdef np.ndarray[np.uint8_t, ndim=4] mask
        cdef oct_visitors.MaskOcts visitor
        visitor = oct_visitors.MaskOcts(self, domain_id)
        cdef int ns = 1 << self.oref
        mask = np.zeros((num_cells, ns, ns, ns), dtype="uint8")
        visitor.mask = mask
        self.visit_all_octs(selector, visitor)
        return mask.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef oct_visitors.ICoordsOcts visitor
        visitor = oct_visitors.ICoordsOcts(self, domain_id)
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="int64")
        visitor.icoords = coords
        self.visit_all_octs(selector, visitor)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        cdef int i
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef oct_visitors.IResOcts visitor
        visitor = oct_visitors.IResOcts(self, domain_id)
        #Return the 'resolution' of each cell; ie the level
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.empty(num_cells, dtype="int64")
        visitor.ires = res
        self.visit_all_octs(selector, visitor)
        if self.level_offset > 0:
            for i in range(num_cells):
                res[i] += self.level_offset
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fwidth(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef oct_visitors.FWidthOcts visitor
        visitor = oct_visitors.FWidthOcts(self, domain_id)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth
        fwidth = np.empty((num_cells, 3), dtype="float64")
        visitor.fwidth = fwidth
        self.visit_all_octs(selector, visitor)
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            fwidth[:,i] *= base_dx
        return fwidth

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fcoords(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef oct_visitors.FCoordsOcts visitor
        visitor = oct_visitors.FCoordsOcts(self, domain_id)
        #Return the floating point unitary position of every cell
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="float64")
        visitor.fcoords = coords
        self.visit_all_octs(selector, visitor)
        cdef int i
        cdef np.float64_t base_dx
        for i in range(3):
            base_dx = (self.DRE[i] - self.DLE[i])/self.nn[i]
            coords[:,i] *= base_dx
            coords[:,i] += self.DLE[i]
        return coords

    def save_octree(self):
        # Get the header
        header = dict(dims = (self.nn[0], self.nn[1], self.nn[2]),
                      left_edge = (self.DLE[0], self.DLE[1], self.DLE[2]),
                      right_edge = (self.DRE[0], self.DRE[1], self.DRE[2]),
                      over_refine = self.oref,
                      partial_coverage = self.partial_coverage)
        cdef SelectorObject selector = AlwaysSelector(None)
        # domain_id = -1 here, because we want *every* oct
        cdef oct_visitors.StoreOctree visitor
        visitor = oct_visitors.StoreOctree(self, -1)
        visitor.oref = 0
        visitor.nz = 1
        cdef np.ndarray[np.uint8_t, ndim=1] ref_mask
        ref_mask = np.zeros(self.nocts * visitor.nz, dtype="uint8") - 1
        visitor.ref_mask = ref_mask
        # Enforce partial_coverage here
        self.visit_all_octs(selector, visitor, 1)
        header['octree'] = ref_mask
        return header

    def selector_fill(self, SelectorObject selector,
                      np.ndarray source,
                      np.ndarray dest = None,
                      np.int64_t offset = 0, int dims = 1,
                      int domain_id = -1):
        # This is actually not correct.  The hard part is that we need to
        # iterate the same way visit_all_octs does, but we need to track the
        # number of octs total visited.
        cdef np.int64_t num_cells = -1
        if dest is None:
            # Note that RAMSES can have partial refinement inside an Oct.  This
            # means we actually do want the number of Octs, not the number of
            # cells.
            num_cells = selector.count_oct_cells(self, domain_id)
            dest = np.zeros((num_cells, dims), dtype=source.dtype,
                            order='C')
        if dims != 1:
            raise RuntimeError
        # Just make sure that we're in the right shape.  Ideally this will not
        # duplicate memory.  Since we're in Cython, we want to avoid modifying
        # the .shape attributes directly.
        dest = dest.reshape((num_cells, 1))
        source = source.reshape((source.shape[0], source.shape[1],
                    source.shape[2], source.shape[3], dims))
        cdef OctVisitor visitor
        cdef oct_visitors.CopyArrayI64 visitor_i64
        cdef oct_visitors.CopyArrayF64 visitor_f64
        if source.dtype != dest.dtype:
            raise RuntimeError
        if source.dtype == np.int64:
            visitor_i64 = oct_visitors.CopyArrayI64(self, domain_id)
            visitor_i64.source = source
            visitor_i64.dest = dest
            visitor = visitor_i64
        elif source.dtype == np.float64:
            visitor_f64 = oct_visitors.CopyArrayF64(self, domain_id)
            visitor_f64.source = source
            visitor_f64.dest = dest
            visitor = visitor_f64
        else:
            raise NotImplementedError
        visitor.index = offset
        # We only need this so we can continue calculating the offset
        visitor.dims = dims
        self.visit_all_octs(selector, visitor)
        if (visitor.global_index + 1) * visitor.nz * visitor.dims > source.size:
            print("GLOBAL INDEX RAN AHEAD.",)
            print (visitor.global_index + 1) * visitor.nz * visitor.dims - source.size
            print(dest.size, source.size, num_cells)
            raise RuntimeError
        if visitor.index > dest.size:
            print("DEST INDEX RAN AHEAD.",)
            print(visitor.index - dest.size)
            print (visitor.global_index + 1) * visitor.nz * visitor.dims, source.size
            print(num_cells)
            raise RuntimeError
        if num_cells >= 0:
            return dest
        return visitor.index - offset

    def domain_ind(self, selector, int domain_id = -1):
        cdef np.ndarray[np.int64_t, ndim=1] ind
        # Here's where we grab the masked items.
        ind = np.full(self.nocts, -1, 'int64')
        cdef oct_visitors.IndexOcts visitor
        visitor = oct_visitors.IndexOcts(self, domain_id)
        visitor.oct_index = ind
        self.visit_all_octs(selector, visitor)
        return ind

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel,
            np.ndarray[np.float64_t, ndim=2] pos,
            int skip_boundary = 1,
            int count_boundary = 0):
        cdef int no, p, i
        cdef int ind[3]
        cdef int nb = 0
        cdef Oct *cur
        cdef np.float64_t pp[3]
        cdef np.float64_t cp[3]
        cdef np.float64_t dds[3]
        no = pos.shape[0] #number of octs
        if curdom > self.num_domains: return 0
        cdef OctAllocationContainer *cont = self.domains.get_cont(curdom - 1)
        cdef int initial = cont.n_assigned
        cdef int in_boundary = 0
        # How do we bootstrap ourselves?
        for p in range(no):
            #for every oct we're trying to add find the
            #floating point unitary position on this level
            in_boundary = 0
            for i in range(3):
                pp[i] = pos[p, i]
                dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
                ind[i] = <np.int64_t> ((pp[i] - self.DLE[i])/dds[i])
                cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
                if ind[i] < 0 or ind[i] >= self.nn[i]:
                    in_boundary = 1
            if skip_boundary == in_boundary == 1:
                nb += count_boundary
                continue
            cur = self.next_root(curdom, ind)
            if cur == NULL: raise RuntimeError
            # Now we find the location we want
            # Note that RAMSES I think 1-findiceses levels, but we don't.
            for _ in range(curlevel):
                # At every level, find the cell this oct
                # lives inside
                for i in range(3):
                    #as we get deeper, oct size halves
                    dds[i] = dds[i] / 2.0
                    if cp[i] > pp[i]:
                        ind[i] = 0
                        cp[i] -= dds[i]/2.0
                    else:
                        ind[i] = 1
                        cp[i] += dds[i]/2.0
                # Check if it has not been allocated
                cur = self.next_child(curdom, ind, cur)
            # Now we should be at the right level
            cur.domain = curdom
            cur.file_ind = p
        return cont.n_assigned - initial + nb

    def allocate_domains(self, domain_counts):
        cdef int count, i
        self.num_domains = len(domain_counts) # 1-indexed
        for i, count in enumerate(domain_counts):
            self.domains.append(count)

    cdef void append_domain(self, np.int64_t domain_count):
        self.num_domains += 1
        self.domains.append(domain_count)

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef Oct *next = self.root_mesh[ind[0]][ind[1]][ind[2]]
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains.get_cont(domain_id - 1)
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_objs[cont.n_assigned]
        cont.n_assigned += 1
        self.root_mesh[ind[0]][ind[1]][ind[2]] = next
        self.nocts += 1
        return next

    cdef Oct* next_child(self, int domain_id, int ind[3], Oct *parent) except? NULL:
        cdef int i
        cdef Oct *next = NULL
        if parent.children != NULL:
            next = parent.children[cind(ind[0],ind[1],ind[2])]
        else:
            # This *8 does NOT need to be made generic.
            parent.children = <Oct **> malloc(sizeof(Oct *) * 8)
            for i in range(8):
                parent.children[i] = NULL
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains.get_cont(domain_id - 1)
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_objs[cont.n_assigned]
        cont.n_assigned += 1
        parent.children[cind(ind[0],ind[1],ind[2])] = next
        self.nocts += 1
        return next

    def file_index_octs(self, SelectorObject selector, int domain_id,
                        num_cells = -1):
        # We create oct arrays of the correct size
        cdef np.int64_t i
        cdef np.ndarray[np.uint8_t, ndim=1] levels
        cdef np.ndarray[np.uint8_t, ndim=1] cell_inds
        cdef np.ndarray[np.int64_t, ndim=1] file_inds
        if num_cells < 0:
            num_cells = selector.count_oct_cells(self, domain_id)
        # Initialize variables with dummy values
        levels = np.full(num_cells, 255, dtype="uint8")
        file_inds = np.full(num_cells, -1, dtype="int64")
        cell_inds = np.full(num_cells, 8, dtype="uint8")
        cdef oct_visitors.FillFileIndicesO visitor_o
        cdef oct_visitors.FillFileIndicesR visitor_r
        if self.fill_style == "r":
            visitor_r = oct_visitors.FillFileIndicesR(self, domain_id)
            visitor_r.levels = levels
            visitor_r.file_inds = file_inds
            visitor_r.cell_inds = cell_inds
            visitor = visitor_r
        elif self.fill_style == "o":
            visitor_o = oct_visitors.FillFileIndicesO(self, domain_id)
            visitor_o.levels = levels
            visitor_o.file_inds = file_inds
            visitor_o.cell_inds = cell_inds
            visitor = visitor_o
        else:
            raise RuntimeError
        self.visit_all_octs(selector, visitor)
        return levels, cell_inds, file_inds

    def morton_index_octs(self, SelectorObject selector, int domain_id,
                          num_cells = -1):
        cdef np.int64_t i
        cdef np.uint8_t[:] levels
        cdef np.uint64_t[:] morton_inds
        if num_cells < 0:
            num_cells = selector.count_oct_cells(self, domain_id)
        levels = np.zeros(num_cells, dtype="uint8")
        morton_inds = np.zeros(num_cells, dtype="uint64")
        for i in range(num_cells):
            levels[i] = 100
            morton_inds[i] = 0
        cdef oct_visitors.MortonIndexOcts visitor
        visitor = oct_visitors.MortonIndexOcts(self, domain_id)
        visitor.level_arr = levels
        visitor.morton_ind = morton_inds
        self.visit_all_octs(selector, visitor)
        return levels, morton_inds

    def domain_count(self, SelectorObject selector):
        # We create oct arrays of the correct size
        cdef np.ndarray[np.int64_t, ndim=1] domain_counts
        domain_counts = np.zeros(self.num_domains, dtype="int64")
        cdef oct_visitors.CountByDomain visitor
        visitor = oct_visitors.CountByDomain(self, -1)
        visitor.domain_counts = domain_counts
        self.visit_all_octs(selector, visitor)
        return domain_counts

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level(self, int level,
                   np.ndarray[np.uint8_t, ndim=1] levels,
                   np.ndarray[np.uint8_t, ndim=1] cell_inds,
                   np.ndarray[np.int64_t, ndim=1] file_inds,
                   dest_fields, source_fields,
                   np.int64_t offset = 0):
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef int i, lvl

        for key in dest_fields:
            dest = dest_fields[key]
            source = source_fields[key]
            for i in range(levels.shape[0]):
                lvl = levels[i]
                if lvl != level: continue
                if file_inds[i] < 0:
                    dest[i + offset] = np.nan
                else:
                    dest[i + offset] = source[file_inds[i], cell_inds[i]]

    def fill_index(self, SelectorObject selector = AlwaysSelector(None)):
        """Get the on-file index of each cell"""
        cdef StoreIndex visitor

        cdef np.int64_t[:, :, :, :] cell_inds

        cell_inds = np.full((self.nocts, 2, 2, 2), -1, dtype=np.int64)

        visitor = StoreIndex(self, -1)
        visitor.cell_inds = cell_inds

        self.visit_all_octs(selector, visitor)

        return np.asarray(cell_inds)

    def fill_octcellindex_neighbours(self, SelectorObject selector, int num_octs=-1, int domain_id=-1, int n_ghost_zones=1):
        """Compute the oct and cell indices of all the cells within all selected octs, extended
        by one cell in all directions (for ghost zones computations).

        Parameters
        ----------
        selector : SelectorObject
            Selector for the octs to compute neighbour of
        num_octs : int, optional
            The number of octs to read in
        domain_id : int, optional
            The domain to perform the selection over

        Returns
        -------
        oct_inds : int64 ndarray (nocts*8, )
            The on-domain index of the octs containing each cell
        cell_inds : uint8 ndarray (nocts*8, )
            The index of the cell in its parent oct

        Note
        ----
        oct_inds/cell_inds
        """
        if num_octs == -1:
            num_octs = selector.count_octs(self, domain_id)

        cdef NeighbourCellIndexVisitor visitor

        cdef np.uint8_t[::1] cell_inds
        cdef np.int64_t[::1] oct_inds

        cell_inds = np.full(num_octs*4**3, 8, dtype=np.uint8)
        oct_inds = np.full(num_octs*4**3, -1, dtype=np.int64)

        visitor = NeighbourCellIndexVisitor(self, -1, n_ghost_zones)
        visitor.cell_inds = cell_inds
        visitor.domain_inds = oct_inds

        self.visit_all_octs(selector, visitor)

        return np.asarray(oct_inds), np.asarray(cell_inds)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def fill_level_with_domain(
                   self, int level,
                   np.uint8_t[:] levels,
                   np.uint8_t[:] cell_inds,
                   np.int64_t[:] file_inds,
                   np.int32_t[:] domains,
                   dict dest_fields,
                   dict source_fields,
                   np.int32_t domain,
                   np.int64_t offset = 0
                   ):
        """Similar to fill_level but accepts a domain argument.

        This is particularly useful for frontends that have buffer zones at CPU boundaries.
        These buffer oct cells have a different domain than the local one and
        are usually not read, but one has to read them e.g. to compute ghost zones.
        """
        cdef np.ndarray[np.float64_t, ndim=2] source
        cdef np.ndarray[np.float64_t, ndim=1] dest
        cdef int i, count, lev
        cdef np.int32_t dom

        for key in dest_fields:
            dest = dest_fields[key]
            source = source_fields[key]
            count = 0
            for i in range(levels.shape[0]):
                lev = levels[i]
                dom = domains[i]
                if lev != level or dom != domain: continue
                count += 1
                if file_inds[i] < 0:
                    dest[i + offset] = np.nan
                else:
                    dest[i + offset] = source[file_inds[i], cell_inds[i]]
        return count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def file_index_octs_with_ghost_zones(
            self, SelectorObject selector, int domain_id,
            int num_cells=1, int n_ghost_zones=1):
        """Similar as file_index_octs, but returns the level, cell index,
        file index and domain of the neighbouring cells as well.

        Arguments
        ---------
        selector : SelectorObject
            The selector object. It is expected to select all cells for a given oct.
        domain_id : int
            The domain to select. Set to -1 to select all domains.
        num_cells : int, optional
            The total number of cells (accounting for a 1-cell thick ghost zone layer).

        Returns
        -------
        levels : uint8, shape (num_cells,)
            The level of each cell of the super oct
        cell_inds : uint8, shape (num_cells, )
            The index of each cell of the super oct within its own oct
        file_inds : int64, shape (num_cells, )
            The on-file position of the cell. See notes below.
        domains : int32, shape (num_cells)
            The domain to which the cells belongs. See notes below.

        Notes
        -----

        The algorithm constructs a "super-oct" around each oct (see sketch below,
        where the original oct cells are marked with an x).

        Note that for sparse octrees (such as RAMSES'), the neighbouring cells
        may belong to another domain (this is stored in `domains`). If the dataset
        provides buffer zones between domains (such as RAMSES), this may be stored
        locally and can be accessed directly.


        +---+---+---+---+
        |   |   |   |   |
        |---+---+---+---|
        |   | x | x |   |
        |---+---+---+---|
        |   | x | x |   |
        |---+---+---+---|
        |   |   |   |   |
        +---+---+---+---+

        """
        cdef np.int64_t i
        cdef int num_octs
        if num_cells < 0:
            num_octs = selector.count_octs(self, domain_id)
            num_cells = num_octs * 4**3
        cdef NeighbourCellVisitor visitor

        cdef np.ndarray[np.uint8_t, ndim=1] levels
        cdef np.ndarray[np.uint8_t, ndim=1] cell_inds
        cdef np.ndarray[np.int64_t, ndim=1] file_inds
        cdef np.ndarray[np.int32_t, ndim=1] domains
        levels = np.full(num_cells, 255, dtype="uint8")
        file_inds = np.full(num_cells, -1, dtype="int64")
        cell_inds = np.full(num_cells, 8, dtype="uint8")
        domains = np.full(num_cells, -1, dtype="int32")

        visitor = NeighbourCellVisitor(self, -1, n_ghost_zones)
        # output: level, file_ind and cell_ind of the neighbouring cells
        visitor.levels = levels
        visitor.file_inds = file_inds
        visitor.cell_inds = cell_inds
        visitor.domains = domains
        # direction to explore and extra parameters of the visitor
        visitor.octree = self
        visitor.last = -1

        # Compute indices
        self.visit_all_octs(selector, visitor)

        return levels, cell_inds, file_inds, domains

    def finalize(self):
        cdef SelectorObject selector = AlwaysSelector(None)
        cdef oct_visitors.AssignDomainInd visitor
        visitor = oct_visitors.AssignDomainInd(self, 1)
        self.visit_all_octs(selector, visitor)
        assert ((visitor.global_index+1)*visitor.nz == visitor.index)

cdef int root_node_compare(const void *a, const void *b) nogil:
    cdef OctKey *ao
    cdef OctKey *bo
    ao = <OctKey *>a
    bo = <OctKey *>b
    if ao.key < bo.key:
        return -1
    elif ao.key == bo.key:
        return 0
    else:
        return 1

cdef class SparseOctreeContainer(OctreeContainer):

    def __init__(self, domain_dimensions, domain_left_edge, domain_right_edge,
                 over_refine = 1):
        cdef int i
        self.partial_coverage = 1
        self.oref = over_refine
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
        self.domains = OctObjectPool()
        self.num_domains = 0
        self.level_offset = 0
        self.nocts = 0 # Increment when initialized
        self.root_mesh = NULL
        self.root_nodes = NULL
        self.tree_root = NULL
        self.num_root = 0
        self.max_root = 0
        # We don't initialize the octs yet
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid
        self.fill_style = "r"

    @classmethod
    def load_octree(cls, header):
        raise NotImplementedError

    def save_octree(self):
        raise NotImplementedError

    cdef int get_root(self, int ind[3], Oct **o) nogil:
        o[0] = NULL
        cdef np.int64_t key = self.ipos_to_key(ind)
        cdef OctKey okey
        cdef OctKey **oresult = NULL
        okey.key = key
        okey.node = NULL
        oresult = <OctKey **> tfind(<void*>&okey,
            &self.tree_root, root_node_compare)
        if oresult != NULL:
            o[0] = oresult[0].node
            return 1
        return 0

    cdef void key_to_ipos(self, np.int64_t key, np.int64_t pos[3]):
        # Note: this is the result of doing
        # for i in range(20):
        #     ukey |= (1 << i)
        cdef np.int64_t ukey = 1048575
        cdef int j
        for j in range(3):
            pos[2 - j] = (<np.int64_t>(key & ukey))
            key = key >> 20

    cdef np.int64_t ipos_to_key(self, int pos[3]) nogil:
        # We (hope) that 20 bits is enough for each index.
        cdef int i
        cdef np.int64_t key = 0
        for i in range(3):
            # Note the casting here.  Bitshifting can cause issues otherwise.
            key |= ((<np.int64_t>pos[i]) << 20 * (2 - i))
        return key

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        OctVisitor visitor,
                        int vc = -1,
                        np.int64_t *indices = NULL):
        cdef int i, j, k, n
        cdef np.int64_t key, ukey
        visitor.global_index = -1
        visitor.level = 0
        if vc == -1:
            vc = self.partial_coverage
        cdef np.float64_t pos[3]
        cdef np.float64_t dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i]) / self.nn[i]
        # Pos is the center of the octs
        cdef Oct *o
        for i in range(self.num_root):
            o = self.root_nodes[i].node
            key = self.root_nodes[i].key
            self.key_to_ipos(key, visitor.pos)
            for j in range(3):
                pos[j] = self.DLE[j] + (visitor.pos[j] + 0.5) * dds[j]
            selector.recursively_visit_octs(
                o, pos, dds, 0, visitor, vc)
            if indices != NULL:
                indices[i] = visitor.index

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0 # We no longer have a domain offset.

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef Oct *next = NULL
        self.get_root(ind, &next)
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains.get_cont(domain_id - 1)
        if cont.n_assigned >= cont.n:
            print("Too many assigned.")
            return NULL
        if self.num_root >= self.max_root:
            print("Too many roots.")
            return NULL
        next = &cont.my_objs[cont.n_assigned]
        cont.n_assigned += 1
        cdef np.int64_t key = 0
        cdef OctKey *ikey = &self.root_nodes[self.num_root]
        key = self.ipos_to_key(ind)
        self.root_nodes[self.num_root].key = key
        self.root_nodes[self.num_root].node = next
        tsearch(<void*>ikey, &self.tree_root, root_node_compare)
        self.num_root += 1
        self.nocts += 1
        return next

    def allocate_domains(self, domain_counts, int root_nodes):
        OctreeContainer.allocate_domains(self, domain_counts)
        self.root_nodes = <OctKey*> malloc(sizeof(OctKey) * root_nodes)
        self.max_root = root_nodes
        for i in range(root_nodes):
            self.root_nodes[i].key = -1
            self.root_nodes[i].node = NULL

    def __dealloc__(self):
        # This gets called BEFORE the superclass deallocation.  But, both get
        # called.
        if self.root_nodes != NULL: free(self.root_nodes)

cdef class ARTOctreeContainer(OctreeContainer):
    def __init__(self, oct_domain_dimensions, domain_left_edge,
                 domain_right_edge, partial_coverage = 0,
                 over_refine = 1):
        OctreeContainer.__init__(self, oct_domain_dimensions,
                domain_left_edge, domain_right_edge, partial_coverage,
                 over_refine)
        self.fill_style = "r"

cdef OctList *OctList_subneighbor_find(OctList *olist, Oct *top,
                                       int i, int j, int k):
    if top.children == NULL: return olist
    # The i, j, k here are the offsets of "top" with respect to
    # the oct for whose neighbors we are searching.
    # Note that this will be recursively called.  We will evaluate either 1, 2,
    # or 4 octs for children and potentially adding them.  In fact, this will
    # be 2**(num_zero) where num_zero is the number of indices that are equal
    # to zero; i.e., the number of dimensions along which we are aligned.
    # For now, we assume we will not be doing this along all three zeros,
    # because that would be pretty tricky.
    if i == j == k == 1: return olist
    cdef np.int64_t n[3]
    cdef np.int64_t ind[3]
    cdef np.int64_t off[3][2]
    cdef np.int64_t ii, ij, ik, ci
    ind[0] = 1 - i
    ind[1] = 1 - j
    ind[2] = 1 - k
    for ii in range(3):
        if ind[ii] == 0:
            n[ii] = 2
            off[ii][0] = 0
            off[ii][1] = 1
        elif ind[ii] == -1:
            n[ii] = 1
            off[ii][0] = 1
        elif ind[ii] == 1:
            n[ii] = 1
            off[ii][0] = 0
    for ii in range(n[0]):
        for ij in range(n[1]):
            for ik in range(n[2]):
                ci = cind(off[0][ii], off[1][ij], off[2][ik])
                cand = top.children[ci]
                if cand.children != NULL:
                    olist = OctList_subneighbor_find(olist,
                        cand, i, j, k)
                else:
                    olist = OctList_append(olist, cand)
    return olist

cdef OctList *OctList_append(OctList *olist, Oct *o):
    cdef OctList *this = olist
    if this == NULL:
        this = <OctList *> malloc(sizeof(OctList))
        this.next = NULL
        this.o = o
        return this
    while this.next != NULL:
        this = this.next
    this.next = <OctList*> malloc(sizeof(OctList))
    this = this.next
    this.o = o
    this.next = NULL
    return this

cdef int OctList_count(OctList *olist):
    cdef OctList *this = olist
    cdef int i = 0 # Count the list
    while this != NULL:
        i += 1
        this = this.next
    return i

cdef void OctList_delete(OctList *olist):
    cdef OctList *next
    cdef OctList *this = olist
    while this != NULL:
        next = this.next
        free(this)
        this = next

cdef class OctObjectPool(ObjectPool):
    # This is an inherited version of the ObjectPool that provides setup and
    # teardown functions for the individually allocated objects.  These allow
    # us to initialize the Octs to default values, and we can also free any
    # allocated memory in them.  Implementing _con_to_array also provides the
    # opportunity to supply views of the octs in Python code.
    def __cinit__(self):
        # Base class will ALSO be called
        self.itemsize = sizeof(Oct)
        assert(sizeof(OctAllocationContainer) == sizeof(AllocationContainer))

    cdef void setup_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                         np.int64_t con_id):
        cdef np.uint64_t i
        cdef Oct* octs = <Oct *> obj
        for n in range(n):
            octs[n].file_ind = octs[n].domain = - 1
            octs[n].domain_ind = n + offset
            octs[n].children = NULL

    cdef void teardown_objs(self, void *obj, np.uint64_t n, np.uint64_t offset,
                           np.int64_t con_id):
        cdef np.uint64_t i, j
        cdef Oct *my_octs = <Oct *> obj
        for i in range(n):
            if my_octs[i].children != NULL:
                free(my_octs[i].children)
        free(obj)

    def _con_to_array(self, int i):
        cdef AllocationContainer *obj = &self.containers[i]
        if obj.n_assigned == 0:
            return None
        cdef OctPadded[:] mm = <OctPadded[:obj.n_assigned]> (
                <OctPadded*> obj.my_objs)
        rv = np.asarray(mm)
        return rv
