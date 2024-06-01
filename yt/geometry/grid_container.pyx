# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
"""
Matching points on the grid to specific grids



"""


import numpy as np

cimport cython
cimport numpy as np

from yt.utilities.lib.bitarray cimport bitarray


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef GridTreeNode Grid_initialize(np.ndarray[np.float64_t, ndim=1] le,
                                  np.ndarray[np.float64_t, ndim=1] re,
                                  np.ndarray[np.int32_t, ndim=1] dims,
                                  int num_children, int level, int index):

    cdef GridTreeNode node
    cdef int i

    node.index = index
    node.level = level
    for i in range(3):
        node.left_edge[i] = le[i]
        node.right_edge[i] = re[i]
        node.dims[i] = dims[i]
        node.dds[i] = (re[i] - le[i])/dims[i]
        node.start_index[i] = <np.int64_t> rint(le[i] / node.dds[i])
    node.num_children = num_children
    if num_children <= 0:
        node.children = NULL
        return node
    node.children = <GridTreeNode **> malloc(
            sizeof(GridTreeNode *) * num_children)
    for i in range(num_children):
        node.children[i] = NULL

    return node

cdef class GridTree:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, int num_grids,
                  np.ndarray[np.float64_t, ndim=2] left_edge,
                  np.ndarray[np.float64_t, ndim=2] right_edge,
                  np.ndarray[np.int32_t, ndim=2] dimensions,
                  list[list[float]] parent_ind, # Not an array anymore!
                  np.ndarray[np.int64_t, ndim=1] level,
                  np.ndarray[np.int64_t, ndim=1] num_children,
                  refine_by = 2,
                  int min_level = 0):

        cdef int i, j, k
        cdef np.ndarray[np.int64_t, ndim=1] child_ptr
        self.total_size = dimensions.prod(axis=1).sum()
        # This doesn't happen all that often, but it is kind of annoying that
        # we don't have uniform "refine_by-as-a-list" behavior.
        try:
            for i, rf in enumerate(np.squeeze(refine_by)):
                self.refine_by[i] = rf
        except TypeError:
            # Not iterable
            self.refine_by[0] = self.refine_by[1] = self.refine_by[2] = refine_by

        child_ptr = np.zeros(num_grids, dtype='int64')

        self.num_grids = num_grids
        self.num_root_grids = 0
        self.num_leaf_grids = 0

        self.grids = <GridTreeNode *> malloc(
                sizeof(GridTreeNode) * num_grids)

        for i in range(num_grids):
            self.grids[i] = Grid_initialize(left_edge[i,:],
                                            right_edge[i,:],
                                            dimensions[i,:],
                                            num_children[i],
                                            level[i], i)
            if level[i] == min_level:
                self.num_root_grids += 1
            if num_children[i] == 0:
                self.num_leaf_grids += 1

        self.root_grids = <GridTreeNode *> malloc(
                sizeof(GridTreeNode) * self.num_root_grids)
        k = 0
        for i in range(num_grids):
            for j in parent_ind[i]:
                if j >= 0:
                    self.grids[j].children[child_ptr[j]] = &self.grids[i]
                    child_ptr[j] += 1
                else:
                    if level[i] < min_level: continue
                    if k >= self.num_root_grids:
                        raise RuntimeError(k, self.num_root_grids)
                    self.root_grids[k] = self.grids[i]
                    k = k + 1

    def __iter__(self):
        yield self

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def return_tree_info(self):
        cdef int i, j
        levels = []
        indices = []
        nchild = []
        children = []
        for i in range(self.num_grids):
            childs = []
            levels.append(self.grids[i].level)
            indices.append(self.grids[i].index)
            nchild.append(self.grids[i].num_children)
            for j in range(self.grids[i].num_children):
                childs.append(self.grids[i].children[j].index)
            children.append(childs)
        return indices, levels, nchild, children

    @property
    def grid_arrays(self):
        cdef GridTreeNodePadded[:] grids
        grids = <GridTreeNodePadded[:self.num_grids]> \
            (<GridTreeNodePadded*> self.grids)
        grids_basic = np.asarray(grids)
        # This next bit is necessary because as of 0.23.4, Cython can't make
        # nested dtypes automatically where you have a property that is
        # something like float[3].  So we unroll all of those, then re-roll
        # them in a new dtype.
        dtn = {}
        dt = grids_basic.dtype
        for name in dt.names:
            d, o = dt.fields[name]
            n = name
            if name.endswith("_x"):
                f = (d.char, 3)
                n = name[:-2]
            elif name.endswith("_y") or name.endswith("_z"):
                continue
            else:
                f = d.char
            dtn[n] = (f, o)
        return grids_basic.view(dtype=np.dtype(dtn))

    cdef void setup_data(self, GridVisitorData *data) noexcept nogil:
        # Being handed a new GVD object, we initialize it to sane defaults.
        data.index = 0
        data.global_index = 0
        data.n_tuples = 0
        data.child_tuples = NULL
        data.array = NULL
        cdef int i
        for i in range(3):
            data.ref_factor[i] = self.refine_by[i]

    cdef void visit_grids(self, GridVisitorData *data,
                          grid_visitor_function *func,
                          SelectorObject selector):
        # This iterates over all root grids, given a selector+data, and then
        # visits each one and its children.
        cdef int i
        # An issue that arises every so often is what to do
        # about children that have multiple parents.  We
        # address that here by ensuring that each grid is only
        # visited once.  This is a safer solution than
        # computing the edges of the parent and passing those
        # down (recursively) because that requires a lot of
        # additional math *and* arguments between functions.
        cdef bitarray visited_mask = bitarray(self.num_grids)
        cdef GridTreeNode *grid
        for i in range(self.num_root_grids):
            grid = &self.root_grids[i]
            self.recursively_visit_grid(data, func, selector, grid, visited_mask)
        grid_visitors.free_tuples(data)

    cdef void recursively_visit_grid(self, GridVisitorData *data,
                                     grid_visitor_function *func,
                                     SelectorObject selector,
                                     GridTreeNode *grid,
                                     bitarray visited_mask) noexcept nogil:
        # Visit this grid and all of its child grids, with a given grid visitor
        # function.  We early terminate if we are not selected by the selector.
        cdef int i
        data.grid = grid
        if visited_mask._query_value(grid.index):
            return
        if selector.select_bbox(grid.left_edge, grid.right_edge) == 0:
            # Note that this does not increment the global_index.
            return
        grid_visitors.setup_tuples(data)
        selector.visit_grid_cells(data, func)
        for i in range(grid.num_children):
            self.recursively_visit_grid(data, func, selector, grid.children[i], visited_mask)
        visited_mask._set_value(grid.index, 1)

    def count(self, SelectorObject selector):
        # Use the counting grid visitor
        cdef GridVisitorData data
        self.setup_data(&data)
        cdef bitarray mask = bitarray(self.total_size)
        data.array = <void*>mask.buf
        self.visit_grids(&data, grid_visitors.mask_cells, selector)
        size = 0
        self.setup_data(&data)
        data.array = <void*>(&size)
        self.visit_grids(&data,  grid_visitors.count_cells, selector)
        return size

    def select_icoords(self, SelectorObject selector, np.uint64_t size = -1):
        # Fill icoords with a selector
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.int64_t, ndim=2] icoords
        icoords = np.empty((size, 3), dtype="int64")
        data.array = icoords.data
        self.visit_grids(&data, grid_visitors.icoords_cells, selector)
        return icoords

    def select_ires(self, SelectorObject selector, np.uint64_t size = -1):
        # Fill ires with a selector
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.int64_t, ndim=1] ires
        ires = np.empty(size, dtype="int64")
        data.array = ires.data
        self.visit_grids(&data, grid_visitors.ires_cells, selector)
        return ires

    def select_fcoords(self, SelectorObject selector, np.uint64_t size = -1):
        # Fill fcoords with a selector
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.float64_t, ndim=2] fcoords
        fcoords = np.empty((size, 3), dtype="float64")
        data.array = fcoords.data
        self.visit_grids(&data, grid_visitors.fcoords_cells, selector)
        return fcoords

    def select_fwidth(self, SelectorObject selector, np.uint64_t size = -1):
        # Fill fwidth with a selector
        cdef GridVisitorData data
        self.setup_data(&data)
        if size == -1:
            size = 0
            data.array = <void*>(&size)
            self.visit_grids(&data,  grid_visitors.count_cells, selector)
        cdef np.ndarray[np.float64_t, ndim=2] fwidth
        fwidth = np.empty((size, 3), dtype="float64")
        data.array = fwidth.data
        self.visit_grids(&data, grid_visitors.fwidth_cells, selector)
        return fwidth

cdef class MatchPointsToGrids:

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __cinit__(self, GridTree tree,
                  int num_points,
                  np.ndarray[np.float64_t, ndim=1] x,
                  np.ndarray[np.float64_t, ndim=1] y,
                  np.ndarray[np.float64_t, ndim=1] z):

        cdef int i

        self.num_points = num_points
        self.xp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.yp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.zp = <np.float64_t *> malloc(
                sizeof(np.float64_t) * num_points)
        self.point_grids = <np.int64_t *> malloc(
                sizeof(np.int64_t) * num_points)
        for i in range(num_points):
            self.xp[i] = x[i]
            self.yp[i] = y[i]
            self.zp[i] = z[i]
            self.point_grids[i] = -1
        self.tree = tree

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_points_in_tree(self):
        cdef np.ndarray[np.int64_t, ndim=1] pt_grids
        cdef int i, j
        cdef np.uint8_t in_grid
        pt_grids = np.zeros(self.num_points, dtype='int64')
        for i in range(self.num_points):
            in_grid = 0
            for j in range(self.tree.num_root_grids):
                if not in_grid:
                    in_grid = self.check_position(i, self.xp[i], self.yp[i], self.zp[i],
                                                  &self.tree.root_grids[j])
        for i in range(self.num_points):
            pt_grids[i] = self.point_grids[i]
        return pt_grids

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef np.uint8_t check_position(self,
                                   np.int64_t pt_index,
                                   np.float64_t x,
                                   np.float64_t y,
                                   np.float64_t z,
                                   GridTreeNode * grid):
        cdef int i
        cdef np.uint8_t in_grid
        in_grid = self.is_in_grid(x, y, z, grid)
        if in_grid:
            if grid.num_children > 0:
                in_grid = 0
                for i in range(grid.num_children):
                    if not in_grid:
                        in_grid = self.check_position(pt_index, x, y, z, grid.children[i])
                if not in_grid:
                    self.point_grids[pt_index] = grid.index
                    in_grid = 1
            else:
                self.point_grids[pt_index] = grid.index
                in_grid = 1
        return in_grid

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef np.uint8_t is_in_grid(self,
             np.float64_t x,
             np.float64_t y,
             np.float64_t z,
             GridTreeNode * grid):
        if x >= grid.right_edge[0]: return 0
        if y >= grid.right_edge[1]: return 0
        if z >= grid.right_edge[2]: return 0
        if x < grid.left_edge[0]: return 0
        if y < grid.left_edge[1]: return 0
        if z < grid.left_edge[2]: return 0
        return 1
