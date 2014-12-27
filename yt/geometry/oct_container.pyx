"""
Oct container




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy as np
import numpy as np
from selection_routines cimport SelectorObject
from libc.math cimport floor
cimport selection_routines

ORDER_MAX = 20
_ORDER_MAX = ORDER_MAX

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef OctAllocationContainer *allocate_octs(
        int n_octs, OctAllocationContainer *prev):
    cdef OctAllocationContainer *n_cont
    cdef Oct *oct
    cdef int n, i, j, k
    n_cont = <OctAllocationContainer *> malloc(
        sizeof(OctAllocationContainer))
    if prev == NULL:
        n_cont.offset = 0
    else:
        n_cont.offset = prev.offset + prev.n
    n_cont.my_octs = <Oct *> malloc(sizeof(Oct) * n_octs)
    if n_cont.my_octs == NULL:
        raise MemoryError
    n_cont.n = n_octs
    n_cont.n_assigned = 0
    n_cont.con_id = -1
    for n in range(n_octs):
        oct = &n_cont.my_octs[n]
        oct.file_ind = oct.domain = -1
        oct.domain_ind = n + n_cont.offset
        oct.children = NULL
    if prev != NULL:
        prev.next = n_cont
    n_cont.next = NULL
    return n_cont

cdef void free_octs(
        OctAllocationContainer *first):
    cdef OctAllocationContainer *cur
    while first != NULL:
        cur = first
        for i in range(cur.n):
            if cur.my_octs[i].children != NULL:
                free(cur.my_octs[i].children)
        free(first.my_octs)
        first = cur.next
        free(cur)

# Here is the strategy for RAMSES containers:
#   * Read each domain individually, creating *all* octs found in that domain
#     file, even if they reside on other CPUs.
#   * Only allocate octs that reside on >= domain
#   * For all octs, insert into tree, which may require traversing existing
#     octs
#   * Note that this does not allow OctAllocationContainer to exactly be a
#     chunk, but it is close.  For IO chunking, we can theoretically examine
#     those octs that live inside a given allocator.

cdef class OctreeContainer:

    def __init__(self, oct_domain_dimensions, domain_left_edge,
                 domain_right_edge, partial_coverage = 0,
                 over_refine = 1):
        # This will just initialize the root mesh octs
        self.oref = over_refine
        self.partial_coverage = partial_coverage
        self.cont = NULL
        cdef int i, j, k, p
        for i in range(3):
            self.nn[i] = oct_domain_dimensions[i]
        self.num_domains = 0
        self.level_offset = 0
        self.domains = NULL
        p = 0
        self.nocts = 0 # Increment when initialized
        for i in range(3):
            self.DLE[i] = domain_left_edge[i] #0
            self.DRE[i] = domain_right_edge[i] #num_grid
        self._initialize_root_mesh()
        self.fill_func = oct_visitors.fill_file_indices_oind

    def _initialize_root_mesh(self):
        self.root_mesh = <Oct****> malloc(sizeof(void*) * self.nn[0])
        for i in range(self.nn[0]):
            self.root_mesh[i] = <Oct ***> malloc(sizeof(void*) * self.nn[1])
            for j in range(self.nn[1]):
                self.root_mesh[i][j] = <Oct **> malloc(sizeof(void*) * self.nn[2])
                for k in range(self.nn[2]):
                    self.root_mesh[i][j][k] = NULL

    @classmethod
    def load_octree(cls, header):
        cdef np.ndarray[np.uint8_t, ndim=1] ref_mask
        ref_mask = header['octree']
        cdef OctreeContainer obj = cls(header['dims'], header['left_edge'],
                header['right_edge'], over_refine = header['over_refine'],
                partial_coverage = header['partial_coverage'])
        # NOTE: We do not allow domain/file indices to be specified.
        cdef SelectorObject selector = selection_routines.AlwaysSelector(None)
        cdef OctVisitorData data
        obj.setup_data(&data, -1)
        cdef int i, j, k, n
        data.global_index = -1
        data.level = 0
        data.oref = 0
        data.nz = 1
        assert(ref_mask.shape[0] / float(data.nz) ==
            <int>(ref_mask.shape[0]/float(data.nz)))
        obj.allocate_domains([ref_mask.shape[0] / data.nz])
        cdef np.float64_t pos[3], dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (obj.DRE[i] - obj.DLE[i]) / obj.nn[i]
        # Pos is the center of the octs
        cdef OctAllocationContainer *cur = obj.domains[0]
        cdef Oct *o
        cdef void *p[4]
        cdef np.int64_t nfinest = 0
        p[0] = ref_mask.data
        p[1] = <void *> cur.my_octs
        p[2] = <void *> &cur.n_assigned
        p[3] = <void *> &nfinest
        data.array = p
        pos[0] = obj.DLE[0] + dds[0]/2.0
        for i in range(obj.nn[0]):
            pos[1] = obj.DLE[1] + dds[1]/2.0
            for j in range(obj.nn[1]):
                pos[2] = obj.DLE[2] + dds[2]/2.0
                for k in range(obj.nn[2]):
                    if obj.root_mesh[i][j][k] != NULL:
                        raise RuntimeError
                    o = &cur.my_octs[cur.n_assigned]
                    o.domain_ind = o.file_ind = 0
                    o.domain = 1
                    obj.root_mesh[i][j][k] = o
                    cur.n_assigned += 1
                    data.pos[0] = i
                    data.pos[1] = j
                    data.pos[2] = k
                    # Always visit covered
                    selector.recursively_visit_octs(
                        obj.root_mesh[i][j][k],
                        pos, dds, 0, oct_visitors.load_octree,
                        &data, 1)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]
        obj.nocts = cur.n_assigned
        if obj.nocts * data.nz != ref_mask.size:
            raise KeyError(ref_mask.size, obj.nocts, obj.oref,
                obj.partial_coverage)
        return obj

    cdef void setup_data(self, OctVisitorData *data, int domain_id = -1):
        cdef int i
        data.index = 0
        data.last = -1
        data.global_index = -1
        for i in range(3):
            data.pos[i] = -1
            data.ind[i] = -1
        data.array = NULL
        data.dims = 0
        data.domain = domain_id
        data.level = -1
        data.oref = self.oref
        data.nz = (1 << (data.oref*3))

    def __dealloc__(self):
        free_octs(self.cont)
        if self.root_mesh == NULL: return
        for i in range(self.nn[0]):
            if self.root_mesh[i] == NULL: continue
            for j in range(self.nn[1]):
                if self.root_mesh[i][j] == NULL: continue
                free(self.root_mesh[i][j])
            if self.root_mesh[i] == NULL: continue
            free(self.root_mesh[i])
        free(self.root_mesh)

    def __iter__(self):
        #Get the next oct, will traverse domains
        #Note that oct containers can be sorted 
        #so that consecutive octs are on the same domain
        cdef OctAllocationContainer *cur = self.cont
        cdef Oct *this
        cdef int i
        while cur != NULL:
            for i in range(cur.n_assigned):
                this = &cur.my_octs[i]
                yield (this.file_ind, this.domain_ind, this.domain)
            cur = cur.next

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data,
                        int vc = -1):
        cdef int i, j, k, n
        if vc == -1:
            vc = self.partial_coverage
        data.global_index = -1
        data.level = 0
        cdef np.float64_t pos[3], dds[3]
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
                    data.pos[0] = i
                    data.pos[1] = j
                    data.pos[2] = k
                    selector.recursively_visit_octs(
                        self.root_mesh[i][j][k],
                        pos, dds, 0, func, data, vc)
                    pos[2] += dds[2]
                pos[1] += dds[1]
            pos[0] += dds[0]

    cdef void oct_bounds(self, Oct *o, np.float64_t *corner, np.float64_t *size):
        cdef int i
        #for i in range(3):
        #    size[i] = (self.DRE[i] - self.DLE[i]) / (self.nn[i] << o.level)
        #    corner[i] = o.pos[i] * size[i] + self.DLE[i]

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0

    cdef int get_root(self, int ind[3], Oct **o):
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
                  int max_level = 99):
        #Given a floating point position, retrieve the most
        #refined oct at that time
        cdef int ind32[3]
        cdef np.int64_t ipos[3]
        cdef np.float64_t dds[3], cp[3], pp[3]
        cdef Oct *cur, *next
        cdef int i
        cur = next = NULL
        cdef np.int64_t ind[3], level = -1
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i])/self.nn[i]
            ind[i] = <np.int64_t> (floor((ppos[i] - self.DLE[i])/dds[i]))
            cp[i] = (ind[i] + 0.5) * dds[i] + self.DLE[i]
            ipos[i] = 0 # We add this to ind later, so it should be zero.
            ind32[i] = ind[i]
        self.get_root(ind32, &next)
        # We want to stop recursing when there's nowhere else to go
        while next != NULL and level <= max_level:
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

    def domain_identify(self, SelectorObject selector):
        cdef np.ndarray[np.uint8_t, ndim=1] domain_mask
        domain_mask = np.zeros(self.num_domains, dtype="uint8")
        cdef OctVisitorData data
        self.setup_data(&data)
        data.array = domain_mask.data
        self.visit_all_octs(selector, oct_visitors.identify_octs, &data)
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
        cdef Oct* candidate
        nn = 0
        # We are going to do a brute-force search here.
        # This is not the most efficient -- in fact, it's relatively bad.  But
        # we will attempt to improve it in a future iteration, where we will
        # grow a stack of parent Octs.
        # Note that in the first iteration, we will just find the up-to-27
        # neighbors, including the main oct.
        cdef np.int64_t i, j, k, n, level, ii, nfound = 0, dlevel
        cdef int ind[3]
        cdef OctList *olist, *my_list
        my_list = olist = NULL
        cdef Oct *cand
        cdef np.int64_t npos[3], ndim[3]
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
        cdef np.ndarray[np.uint8_t, ndim=1] coords
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        coords = np.zeros((num_cells*data.nz), dtype="uint8")
        data.array = <void *> coords.data
        self.visit_all_octs(selector, oct_visitors.mask_octs, &data)
        return coords.astype("bool")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def icoords(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        cdef np.ndarray[np.int64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="int64")
        data.array = <void *> coords.data
        self.visit_all_octs(selector, oct_visitors.icoords_octs, &data)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def ires(self, SelectorObject selector, np.int64_t num_cells = -1,
                int domain_id = -1):
        cdef int i
        if num_cells == -1:
            num_cells = selector.count_oct_cells(self, domain_id)
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        #Return the 'resolution' of each cell; ie the level
        cdef np.ndarray[np.int64_t, ndim=1] res
        res = np.empty(num_cells, dtype="int64")
        data.array = <void *> res.data
        self.visit_all_octs(selector, oct_visitors.ires_octs, &data)
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
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth
        fwidth = np.empty((num_cells, 3), dtype="float64")
        data.array = <void *> fwidth.data
        self.visit_all_octs(selector, oct_visitors.fwidth_octs, &data)
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
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        #Return the floating point unitary position of every cell
        cdef np.ndarray[np.float64_t, ndim=2] coords
        coords = np.empty((num_cells, 3), dtype="float64")
        data.array = <void *> coords.data
        self.visit_all_octs(selector, oct_visitors.fcoords_octs, &data)
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
        cdef SelectorObject selector = selection_routines.AlwaysSelector(None)
        # domain_id = -1 here, because we want *every* oct
        cdef OctVisitorData data
        self.setup_data(&data, -1)
        data.oref = 0
        data.nz = 1
        cdef np.ndarray[np.uint8_t, ndim=1] ref_mask
        ref_mask = np.zeros(self.nocts * data.nz, dtype="uint8") - 1
        cdef void *p[1]
        p[0] = ref_mask.data
        data.array = p
        # Enforce partial_coverage here
        self.visit_all_octs(selector, oct_visitors.store_octree, &data, 1)
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
            if dims > 1:
                dest = np.zeros((num_cells, dims), dtype=source.dtype,
                    order='C')
            else:
                dest = np.zeros(num_cells, dtype=source.dtype, order='C')
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        data.index = offset
        # We only need this so we can continue calculating the offset
        data.dims = dims
        cdef void *p[2]
        p[0] = source.data
        p[1] = dest.data
        data.array = &p
        cdef oct_visitor_function *func
        if source.dtype != dest.dtype:
            raise RuntimeError
        if source.dtype == np.int64:
            func = oct_visitors.copy_array_i64
        elif source.dtype == np.float64:
            func = oct_visitors.copy_array_f64
        else:
            raise NotImplementedError
        self.visit_all_octs(selector, func, &data)
        if (data.global_index + 1) * data.nz * data.dims > source.size:
            print "GLOBAL INDEX RAN AHEAD.",
            print (data.global_index + 1) * data.nz * data.dims - source.size
            print dest.size, source.size, num_cells
            raise RuntimeError
        if data.index > dest.size:
            print "DEST INDEX RAN AHEAD.",
            print data.index - dest.size
            print (data.global_index + 1) * data.nz * data.dims, source.size
            print num_cells
            raise RuntimeError
        if num_cells >= 0:
            return dest
        return data.index - offset

    def domain_ind(self, selector, int domain_id = -1):
        cdef np.ndarray[np.int64_t, ndim=1] ind
        # Here's where we grab the masked items.
        ind = np.zeros(self.nocts, 'int64') - 1
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        data.array = ind.data
        self.visit_all_octs(selector, oct_visitors.index_octs, &data)
        return ind

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def add(self, int curdom, int curlevel,
            np.ndarray[np.float64_t, ndim=2] pos,
            int skip_boundary = 1,
            int count_boundary = 0):
        cdef int level, no, p, i, j, k, ind[3]
        cdef int nb = 0
        cdef Oct *cur, *next = NULL
        cdef np.float64_t pp[3], cp[3], dds[3]
        no = pos.shape[0] #number of octs
        if curdom > self.num_domains: return 0
        cdef OctAllocationContainer *cont = self.domains[curdom - 1]
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
            for level in range(curlevel):
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
        cdef OctAllocationContainer *cur = self.cont
        assert(cur == NULL)
        self.num_domains = len(domain_counts) # 1-indexed
        self.domains = <OctAllocationContainer **> malloc(
            sizeof(OctAllocationContainer *) * len(domain_counts))
        for i, count in enumerate(domain_counts):
            cur = allocate_octs(count, cur)
            if self.cont == NULL: self.cont = cur
            self.domains[i] = cur

    cdef void append_domain(self, np.int64_t domain_count):
        self.num_domains += 1
        self.domains = <OctAllocationContainer **> realloc(self.domains, 
                sizeof(OctAllocationContainer *) * self.num_domains)
        if self.domains == NULL: raise RuntimeError
        self.domains[self.num_domains - 1] = NULL
        cdef OctAllocationContainer *cur = NULL
        if self.num_domains > 1:
            cur = self.domains[self.num_domains - 2]
        cur = allocate_octs(domain_count, cur)
        if self.cont == NULL: self.cont = cur
        self.domains[self.num_domains - 1] = cur

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef Oct *next = self.root_mesh[ind[0]][ind[1]][ind[2]]
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        self.root_mesh[ind[0]][ind[1]][ind[2]] = next
        self.nocts += 1
        return next

    cdef Oct* next_child(self, int domain_id, int ind[3], Oct *parent):
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
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n: raise RuntimeError
        next = &cont.my_octs[cont.n_assigned]
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
        levels = np.zeros(num_cells, dtype="uint8")
        file_inds = np.zeros(num_cells, dtype="int64")
        cell_inds = np.zeros(num_cells, dtype="uint8")
        for i in range(num_cells):
            levels[i] = 100
            file_inds[i] = -1
            cell_inds[i] = 9
        cdef OctVisitorData data
        self.setup_data(&data, domain_id)
        cdef void *p[3]
        p[0] = levels.data
        p[1] = file_inds.data
        p[2] = cell_inds.data
        data.array = p
        self.visit_all_octs(selector, self.fill_func, &data)
        return levels, cell_inds, file_inds

    def domain_count(self, SelectorObject selector):
        # We create oct arrays of the correct size
        cdef np.int64_t i, num_octs
        cdef np.ndarray[np.int64_t, ndim=1] domain_counts
        domain_counts = np.zeros(self.num_domains, dtype="int64")
        cdef OctVisitorData data
        self.setup_data(&data, -1)
        data.array = <void*> domain_counts.data
        self.visit_all_octs(selector, oct_visitors.count_by_domain, &data)
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
        cdef int n
        cdef int i, di
        cdef np.int64_t local_pos, local_filled = 0
        cdef np.float64_t val
        for key in dest_fields:
            dest = dest_fields[key]
            source = source_fields[key]
            for i in range(levels.shape[0]):
                if levels[i] != level: continue
                dest[i + offset] = source[file_inds[i], cell_inds[i]]
                local_filled += 1

    def finalize(self):
        cdef SelectorObject selector = selection_routines.AlwaysSelector(None)
        cdef OctVisitorData data
        self.setup_data(&data, 1)
        self.visit_all_octs(selector, oct_visitors.assign_domain_ind, &data)
        assert ((data.global_index+1)*data.nz == data.index)

cdef int root_node_compare(void *a, void *b) nogil:
    cdef OctKey *ao, *bo
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
        cdef int i, j, k, p
        self.partial_coverage = 1
        self.oref = over_refine
        for i in range(3):
            self.nn[i] = domain_dimensions[i]
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
        self.fill_func = oct_visitors.fill_file_indices_rind

    @classmethod
    def load_octree(self, header):
        raise NotImplementedError

    def save_octree(self):
        raise NotImplementedError

    cdef int get_root(self, int ind[3], Oct **o):
        o[0] = NULL
        cdef int i
        cdef np.int64_t key = self.ipos_to_key(ind)
        cdef OctKey okey, **oresult = NULL
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

    cdef np.int64_t ipos_to_key(self, int pos[3]):
        # We (hope) that 20 bits is enough for each index.
        cdef int i
        cdef np.int64_t key = 0
        for i in range(3):
            # Note the casting here.  Bitshifting can cause issues otherwise.
            key |= ((<np.int64_t>pos[i]) << 20 * (2 - i))
        return key

    @cython.cdivision(True)
    cdef void visit_all_octs(self, SelectorObject selector,
                        oct_visitor_function *func,
                        OctVisitorData *data,
                        int vc = -1):
        cdef int i, j, k, n
        cdef np.int64_t key, ukey
        data.global_index = -1
        data.level = 0
        if vc == -1:
            vc = self.partial_coverage
        cdef np.float64_t pos[3], dds[3]
        # This dds is the oct-width
        for i in range(3):
            dds[i] = (self.DRE[i] - self.DLE[i]) / self.nn[i]
        # Pos is the center of the octs
        cdef Oct *o
        for i in range(self.num_root):
            o = self.root_nodes[i].node
            key = self.root_nodes[i].key
            self.key_to_ipos(key, data.pos)
            for j in range(3):
                pos[j] = self.DLE[j] + (data.pos[j] + 0.5) * dds[j]
            selector.recursively_visit_octs(
                o, pos, dds, 0, func, data, vc)

    cdef np.int64_t get_domain_offset(self, int domain_id):
        return 0 # We no longer have a domain offset.

    cdef Oct* next_root(self, int domain_id, int ind[3]):
        cdef int i
        cdef Oct *next = NULL
        self.get_root(ind, &next)
        if next != NULL: return next
        cdef OctAllocationContainer *cont = self.domains[domain_id - 1]
        if cont.n_assigned >= cont.n:
            print "Too many assigned."
            return NULL
        if self.num_root >= self.max_root:
            print "Too many roots."
            return NULL
        next = &cont.my_octs[cont.n_assigned]
        cont.n_assigned += 1
        cdef np.int64_t key = 0
        cdef OctKey *ikey = &self.root_nodes[self.num_root]
        cdef np.int64_t okey = ikey.key
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
        if self.domains != NULL: free(self.domains)

cdef class RAMSESOctreeContainer(SparseOctreeContainer):
    pass

cdef class ARTOctreeContainer(OctreeContainer):
    def __init__(self, oct_domain_dimensions, domain_left_edge,
                 domain_right_edge, partial_coverage = 0,
                 over_refine = 1):
        OctreeContainer.__init__(self, oct_domain_dimensions,
                domain_left_edge, domain_right_edge, partial_coverage,
                 over_refine)
        self.fill_func = oct_visitors.fill_file_indices_rind

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
    cdef np.int64_t n[3], ind[3], off[3][2], ii, ij, ik, ci
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
    cdef OctList *next, *this = olist
    while this != NULL:
        next = this.next
        free(this)
        this = next
