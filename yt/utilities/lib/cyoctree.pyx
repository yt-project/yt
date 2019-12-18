# distutils: libraries = STD_LIBS
"""
CyOctree building, loading and refining routines



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2019, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport libc.math as math
cimport numpy as np
from yt.geometry.particle_deposit cimport \
    kernel_func, get_kernel_func
from libc.stdlib cimport malloc, free
np.import_array()

cdef struct Octree:
    np.float64_t * node_positions
    np.uint8_t * refined
    np.uint8_t * depth
    np.int64_t * pstart
    np.int64_t * pend
    np.int64_t * children

    np.float64_t * pposx
    np.float64_t * pposy
    np.float64_t * pposz
    np.int64_t * pidx

    np.int64_t n_ref
    np.int64_t num_particles

    np.float64_t * size

    np.uint8_t max_depth
    np.int64_t num_nodes
    np.int64_t max_num_nodes


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_build_node(Octree * tree, long int node_idx):
    cdef np.int64_t splits[9]
    cdef np.int64_t i, j, k, n, start, end
    cdef np.float64_t lx, ly, lz, sz

    if ((tree.pend[node_idx] - tree.pstart[node_idx] > tree.n_ref) and
            (tree.depth[node_idx] < tree.max_depth)):

        # If we are running out of space in our tree, then we *try* to relloacate
        # a tree of double the size
        if tree.num_nodes > tree.max_num_nodes - 16:
            octree_reallocate(tree, tree.max_num_nodes * 2)

        tree.refined[node_idx] = 1

        # Order and split the particles into the children
        # this is hardcoded to 8 children, but can be easily changed
        split_helper(tree, node_idx, splits)

        # figure out the size of the current oct
        sx = tree.size[0] / (2**tree.depth[node_idx])
        sy = tree.size[1] / (2**tree.depth[node_idx])
        sz = tree.size[2] / (2**tree.depth[node_idx])
        lx = tree.node_positions[(3*node_idx)] - sx/2
        ly = tree.node_positions[(3*node_idx)+1] - sy/2
        lz = tree.node_positions[(3*node_idx)+2] - sz/2

        # Loop through and generate the children
        n = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    start = splits[n]
                    end = splits[n + 1]
                    child = tree.num_nodes

                    # store the child location
                    tree.children[8*node_idx + n] = child

                    tree.node_positions[(child*3)] = lx + sx*i
                    tree.node_positions[(child*3)+1] = ly + sy*j
                    tree.node_positions[(child*3)+2] = lz + sz*k

                    tree.refined[child] = 0
                    tree.depth[child] = tree.depth[node_idx] + 1

                    tree.pstart[child] = start
                    tree.pend[child] = end
                    tree.num_nodes += 1

                    # Recursively refine child
                    octree_build_node(tree, child)

                    n += 1

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_allocate(Octree * octree, long int num_nodes):
    octree.node_positions = <np.float64_t *> malloc(num_nodes * 3 * sizeof(np.float64_t))

    octree.size = <np.float64_t *> malloc(3 *sizeof(np.float64_t))
    octree.children = <np.int64_t *> malloc(8 * num_nodes * sizeof(np.int64_t))

    octree.pstart = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    octree.pend = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))

    octree.refined = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    octree.depth = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))


    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_reallocate(Octree * octree, long int num_nodes):
    cdef np.float64_t * old_arr
    cdef np.int64_t * old_arr_int
    cdef np.uint8_t * old_arr_uint
    cdef np.int64_t i

    old_arr = octree.node_positions
    octree.node_positions = <np.float64_t *> malloc(num_nodes * 3 * sizeof(np.float64_t))
    for i in range(3*octree.num_nodes):
        octree.node_positions[i] = old_arr[i]
    free(old_arr)

    old_arr_int = octree.children
    octree.children = <np.int64_t *> malloc(num_nodes * 8 * sizeof(np.int64_t))
    for i in range(8*octree.num_nodes):
        octree.children[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_int = octree.pstart
    octree.pstart = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    for i in range(octree.num_nodes):
        octree.pstart[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_int = octree.pend
    octree.pend = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    for i in range(octree.num_nodes):
        octree.pend[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_uint = octree.refined
    octree.refined = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    for i in range(octree.num_nodes):
        octree.refined[i] = old_arr_uint[i]
    free(old_arr_uint)

    old_arr_uint = octree.depth
    octree.depth = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    for i in range(octree.num_nodes):
        octree.depth[i] = old_arr_uint[i]
    free(old_arr_uint)

    octree.max_num_nodes = num_nodes

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_deallocate(Octree * octree):
    free(octree.node_positions)
    free(octree.size)
    free(octree.children)

    free(octree.pstart)
    free(octree.pend)

    free(octree.refined)
    free(octree.depth)

    free(octree.pidx)

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class CyOctree:
    cdef Octree * c_octree
    cdef np.float64_t[::1, :] input_positions

    cdef np.int64_t n_ref
    cdef np.float64_t[:] left_edge
    cdef np.float64_t[:] right_edge
    cdef np.uint8_t max_depth

    cdef kernel_func kernel

    def __init__(self, np.float64_t[:, :] input_pos, left_edge = None,
                 right_edge = None, np.int64_t n_ref=32, np.uint8_t max_depth=200):

        self.n_ref = n_ref
        self.max_depth = max_depth
        self.input_positions = np.asfortranarray(input_pos, dtype="float64")

        self._allocate_octree()
        self._make_root(left_edge, right_edge)

        octree_build_node(self.c_octree, 0)

        octree_reallocate(self.c_octree, self.c_octree.num_nodes)

    @property
    def num_nodes(self):
        return self.c_octree.num_nodes

    @property
    def node_positions(self):
        cdef np.npy_intp shape[2]
        shape[0] = <np.npy_intp> self.c_octree.num_nodes
        shape[1] = 3
        arr = np.PyArray_SimpleNewFromData(2, &shape[0], np.NPY_FLOAT64, <void *>self.c_octree.node_positions)
        return arr

    @property
    def node_refined(self):
        cdef np.npy_intp shape
        shape = <np.npy_intp> self.c_octree.num_nodes
        arr = np.PyArray_SimpleNewFromData(1, &shape, np.NPY_UINT8, <void *>self.c_octree.refined)
        return arr.astype("bool")

    @property
    def node_depth(self):
        cdef np.npy_intp shape
        shape = <np.npy_intp> self.c_octree.num_nodes
        arr = np.PyArray_SimpleNewFromData(1, &shape, np.NPY_UINT8, <void *>self.c_octree.depth)
        return arr

    @property
    def node_sizes(self):
        cdef np.int64_t i
        sizes = np.zeros((self.c_octree.num_nodes, 3), dtype="float64")
        sizes[:, 0] = self.c_octree.size[0]
        sizes[:, 1] = self.c_octree.size[1]
        sizes[:, 2] = self.c_octree.size[2]

        for i in range(self.c_octree.num_nodes):
            sizes[i, :] /= <np.float64_t> ((2**(<np.float64_t>self.c_octree.depth[i]))/2)

        return sizes

    def _make_root(self, left_edge, right_edge):
        cdef int i = 0

        self.c_octree.num_particles = self.input_positions.shape[0]

        self.c_octree.pidx = <np.int64_t *> malloc(self.c_octree.num_particles * sizeof(np.int64_t))
        for i in range(0, self.c_octree.num_particles):
            self.c_octree.pidx[i] = i

        if left_edge is None:
            left_edge = np.zeros(3, dtype="float64")
            right_edge = np.zeros(3, dtype="float64")

            left_edge[0] = self.c_octree.pposx[0]
            left_edge[1] = self.c_octree.pposy[0]
            left_edge[2] = self.c_octree.pposz[0]
            right_edge[0] = self.c_octree.pposx[0]
            right_edge[1] = self.c_octree.pposy[0]
            right_edge[2] = self.c_octree.pposz[0]

            for i in range(self.c_octree.num_particles):
                if self.c_octree.pposx[i] < left_edge[0]:
                    left_edge[0] = self.c_octree.pposx[i]
                if self.c_octree.pposy[i] < left_edge[1]:
                    left_edge[1] = self.c_octree.pposy[i]
                if self.c_octree.pposz[i] < left_edge[2]:
                    left_edge[2] = self.c_octree.pposz[i]

                if self.c_octree.pposx[i] > right_edge[0]:
                    right_edge[0] = self.c_octree.pposx[i]
                if self.c_octree.pposy[i] > right_edge[1]:
                    right_edge[1] = self.c_octree.pposy[i]
                if self.c_octree.pposz[i] > right_edge[2]:
                    right_edge[2] = self.c_octree.pposz[i]

            self.c_octree.pstart[0] = 0
            self.c_octree.pend[0] = self.input_positions.shape[0]
        else:
            left_edge = left_edge.astype("float64")
            right_edge = right_edge.astype("float64")

            split = select(self.c_octree, left_edge, right_edge, 0, self.input_positions.shape[0])
            self.c_octree.pstart[0] = 0
            self.c_octree.pend[0] = split

        size = (right_edge - left_edge) / 2.0
        center = (right_edge + left_edge) / 2.0
        self.left_edge = left_edge
        self.right_edge = right_edge

        self.c_octree.node_positions[0] = center[0]
        self.c_octree.node_positions[1] = center[1]
        self.c_octree.node_positions[2] = center[2]

        self.c_octree.size[0] = size[0]
        self.c_octree.size[1] = size[1]
        self.c_octree.size[2] = size[2]

        self.c_octree.refined[0] = 0
        self.c_octree.depth[0] = 0

    def _allocate_octree(self):
        self.c_octree = <Octree*> malloc(sizeof(Octree))
        self.c_octree.n_ref = self.n_ref
        self.c_octree.num_nodes = 1

        self.c_octree.max_num_nodes = 1_000_000
        self.c_octree.max_depth = self.max_depth

        self.c_octree.pposx = &self.input_positions[0, 0]
        self.c_octree.pposy = &self.input_positions[0, 1]
        self.c_octree.pposz = &self.input_positions[0, 2]

        octree_allocate(self.c_octree, self.c_octree.max_num_nodes)

    cdef void smooth_onto_cells(self, np.float64_t[:] buff,
                                np.float64_t[:] buff_den, np.float64_t posx,
                                np.float64_t posy, np.float64_t posz,
                                np.float64_t hsml, np.float64_t prefactor,
                                np.float64_t prefactor_norm, long int num_node,
                                int use_normalization=0):

        cdef Octree * tree = self.c_octree
        cdef double q_ij, diff_x, diff_y, diff_z, diff, sx, sy, sz, size
        cdef int i
        cdef long int child_node

        if tree.refined[num_node] == 0:
            diff_x = tree.node_positions[3*num_node] - posx
            diff_y = tree.node_positions[3*num_node+1] - posy
            diff_z = tree.node_positions[3*num_node+2] - posz

            q_ij = math.sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z)
            q_ij /= hsml

            buff[num_node] += (prefactor * self.kernel(q_ij))

            if use_normalization:
                buff_den[num_node] += (prefactor_norm * self.kernel(q_ij))

        else:
            for i in range(8):
                child_node = tree.children[8*num_node + i]

                diff_x = tree.node_positions[3*child_node] - posx
                diff_y = tree.node_positions[3*child_node+1] - posy
                diff_z = tree.node_positions[3*child_node+2] - posz
                diff = math.sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z)

                sx = tree.size[0] / (2**tree.depth[child_node])
                sy = tree.size[1] / (2**tree.depth[child_node])
                sz = tree.size[2] / (2**tree.depth[child_node])
                size = math.sqrt(sx*sx + sy*sy + sz*sz)

                if diff - size - hsml < 0:
                    self.smooth_onto_cells(buff, buff_den, posx, posy, posz,
                                    hsml, prefactor, prefactor_norm,
                                    child_node, use_normalization=use_normalization)


    def interpolate_sph_cells(self, np.float64_t[:] buff,
                              np.float64_t[:] buff_den, np.float64_t[:] posx,
                              np.float64_t[:] posy, np.float64_t[:] posz,
                              np.float64_t[:] pmass, np.float64_t[:] pdens,
                              np.float64_t[:] hsml, np.float64_t[:] field,
                              kernel_name="cubic", int use_normalization=0):

        self.kernel = get_kernel_func(kernel_name)

        cdef int i, j
        cdef double prefactor, prefactor_norm
        for i in range(posx.shape[0]):
            prefactor = pmass[i] / pdens[i] / hsml[i]**3
            prefactor_norm = prefactor
            prefactor *= field[i]

            self.smooth_onto_cells(buff, buff_den, posx[i], posy[i], posz[i],
                                    hsml[i], prefactor, prefactor_norm,
                                    0, use_normalization=use_normalization)

    def mpl_project2d(self, ax, np.int64_t ind1=0, np.int64_t ind2=1, kwargs={}):
        def_kwargs = {'width': 0.1, 'edgecolor': 'k', 'facecolor': 'k'}

        for key, item in kwargs.items():
            def_kwargs[key] = item


        import matplotlib.patches as patches

        cdef np.int64_t i = 0
        rects_plotted = {}
        for i in range(self.c_octree.num_nodes):
            # Create a Rectangle patch
            if i == 0:
                size = (
                    2*self.c_octree.size[ind1]/2**(<np.int64_t>self.c_octree.depth[i]),
                    2*self.c_octree.size[ind2]/2**(<np.int64_t>self.c_octree.depth[i]))
                pos = (
                    self.c_octree.node_positions[(3*i)+ind1]-size[0]/2,
                    self.c_octree.node_positions[(3*i)+ind2]-size[1]/2)

                rect = patches.Rectangle(pos, size[0], size[1], linewidth=def_kwargs["width"],
                                         edgecolor=def_kwargs["edgecolor"], facecolor="none")
                ax.add_patch(rect)

            if self.c_octree.refined[i] == 1:
                size = (
                    self.c_octree.size[ind1]/2**(<np.int64_t>self.c_octree.depth[i]),
                    self.c_octree.size[ind2]/2**(<np.int64_t>self.c_octree.depth[i]))
                pos = (
                    self.c_octree.node_positions[(3*i)+ind1],
                    self.c_octree.node_positions[(3*i)+ind2])

                if pos not in rects_plotted:
                    rect = patches.Arrow(pos[0]-size[0], pos[1], 2*size[0], 0.0, **def_kwargs)
                    ax.add_patch(rect)
                    rect = patches.Arrow(pos[0], pos[1]-size[1], 0, 2*size[1], **def_kwargs)
                    ax.add_patch(rect)

                    rects_plotted[pos] = self.c_octree.depth[i]


    def __del__(self):
        octree_deallocate(self.c_octree)
        free(self.c_octree)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate_x(Octree * octree, double value, np.int64_t start, np.int64_t end) nogil except -1:
    cdef np.int64_t index
    cdef np.int64_t split = start

    cdef np.float64_t * pos = octree.pposx
    cdef np.int64_t * pidx = octree.pidx

    for index in range(start, end):
        if pos[pidx[index]] < value:
            pidx[split], pidx[index] = pidx[index], pidx[split]
            split+=1

    return split


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate_y(Octree * octree, double value, np.int64_t start, np.int64_t end) nogil except -1:
    cdef np.int64_t index
    cdef np.int64_t split = start

    cdef np.float64_t * pos = octree.pposy
    cdef np.int64_t * pidx = octree.pidx

    for index in range(start, end):
        if pos[pidx[index]] < value:
            pidx[split], pidx[index] = pidx[index], pidx[split]
            split+=1

    return split


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate_z(Octree * octree, double value, np.int64_t start, np.int64_t end) nogil except -1:
    cdef np.int64_t index
    cdef np.int64_t split = start

    cdef np.float64_t * pos = octree.pposz
    cdef np.int64_t * pidx = octree.pidx

    for index in range(start, end):
        if pos[pidx[index]] < value:
            pidx[split], pidx[index] = pidx[index], pidx[split]
            split+=1

    return split


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int split_helper(Octree * tree, np.int64_t node_idx, np.int64_t * splits):
    splits[0] = tree.pstart[node_idx]
    splits[8] = tree.pend[node_idx]

    splits[4] = seperate_x(tree, tree.node_positions[(3*node_idx)], splits[0], splits[8])

    splits[2] = seperate_y(tree, tree.node_positions[(3*node_idx)+1], splits[0], splits[4])
    splits[6] = seperate_y(tree, tree.node_positions[(3*node_idx)+1], splits[4], splits[8])

    splits[1] = seperate_z(tree, tree.node_positions[(3*node_idx)+2], splits[0], splits[2])
    splits[3] = seperate_z(tree, tree.node_positions[(3*node_idx)+2], splits[2], splits[4])
    splits[5] = seperate_z(tree, tree.node_positions[(3*node_idx)+2], splits[4], splits[6])
    splits[7] = seperate_z(tree, tree.node_positions[(3*node_idx)+2], splits[6], splits[8])


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int select(
        Octree * octree, np.float64_t[::1] left_edge, np.float64_t[::1] right_edge,
        np.int64_t start, np.int64_t end) nogil except -1:
    cdef np.int64_t index
    cdef np.int64_t split = start

    cdef np.float64_t * posx = octree.pposx
    cdef np.float64_t * posy = octree.pposy
    cdef np.float64_t * posz = octree.pposz
    cdef np.int64_t * pidx = octree.pidx

    for index in range(start, end):
        if posx[pidx[index]] < right_edge[0] and posx[pidx[index]] > left_edge[0]:
            if posy[pidx[index]] < right_edge[1] and posy[pidx[index]] > left_edge[1]:
                if posz[pidx[index]] < right_edge[2] and posz[pidx[index]] > left_edge[2]:
                    if split < index:
                        pidx[split], pidx[index] = pidx[index], pidx[split]
                    split+=1

    return split
