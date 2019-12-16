"""
CyOctree building, loading and refining routines



"""


cimport numpy as np
import numpy as np
cimport cython
import struct

from libcpp.vector cimport vector
from libcpp cimport bool
cimport libc.math as math
from libc.stdlib cimport malloc, free

from yt.geometry.particle_deposit cimport \
    kernel_func, get_kernel_func

################################################################################
#                       OCTREE IMPLEMENTATION DETAILS                          #
################################################################################
# The tree is made of of nodes, which are C structs, containing the left edge,
# the right end and other details to traverse the tree (i.e child and parent
# indexes in the nodes container).
#
# The nodes are stored in an STL vector - as such - it makes sense that the
# parent and child addresses are long's which just describe the index in the STL
# vector. i.e to access a parent,
# &nodes[node.parent] will return the address of the parent node
# The children are contiguous in memory, so the first child can be accessed with
# &nodes[node.children] and the second with,
# &nodes[node.children + 1] etc
# In general we avoid memory addresses in favour of indexes so the reallocation
# of the STL doesn't invalidate those.
#
# The tree is built with a non-recursive algorithm. We start by building the
# root node. We enter the root node, and use the split_helper function to split
# into 2^(3*density_factor) children. We recursively work out which particles
# are in each child by sorting the particle positions array so values less than
# the split value are on one side, and values greater on the other. We then only
# need to store the first and last particle in a node and from that we know
# every particle within the node, in the octree sorted positions array.
# NOTE: Despite being called an octree, the tree is actually of splitting into
# different numbers of children
#
# Once the first lot of children have been generated, they will be stored in the
# STL vector. The STL container will now contain the following,
# nodes = [root, child 1, child 2, child 3, ...]
# where only the root has been processed.
# We then loop through this vector, calling the process_node method and storing
# new children as we generate them and then eventually processing those until no
# new children are generated.
#
# A node will split into children if the node contains more particles than
# n_ref. Below this value, the node will not split, this is called a leaf.
#
# To maintain backwards compatibility we also have the cell structure. This
# means that each leaf is split into 2^(3*over_refine_factor) cells. When an SPH
# field is interpolated onto the octree, we interpolate the value at the centre
# of each cell. The cell locations, and particles they contain, are *NOT*
# stored in memory instead these are calculated on the fly when an interpolation
# or cell position request is made. This has appeared to be a good trade between
# memory and performance. Storing the cells, and all the necessary information
# in memory would increase the memory usage by num_leaves *
# 2^(3*over_refine_factor).

#TODO: Add invalidation and setters
#TODO: Add more deposition functions
#TODO: Add parallel features

cdef struct Node:
    double left_edge[3]
    double right_edge[3]

    np.int64_t start                # First particle we store in pos array
    np.int64_t end                  # Last particle we store

    np.int64_t parent               # Index of parent in nodes container
    np.int64_t children             # Index of 1st child, children are
                                    # contiguous
    bool leaf
    np.int64_t node_id              # Not sure if this is even useful
    np.int64_t leaf_id              # This is used in depositions (maybe)
    unsigned char depth

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class CyOctree:
    cdef vector[Node] nodes         # This is an STL container to store the octs
    cdef double _left_edge[3]       # Boundary conditions for the octree
    cdef double _right_edge[3]
    cdef np.int64_t[:] _idx         # Ordering of particles used by the tree

    cdef np.int64_t _data_version   # Used to decide when to re-build a tree

    cdef int _n_ref                 # Max number of particles per leaf
    cdef np.int64_t _num_particles

    # Cell structure
    cdef int _max_depth
    cdef int _over_refine_factor    # this allows the tree to be built with more
                                    # than 8 cells per leaf
    cdef int _num_cells             # 2**(3*_over_refine_factor)
    cdef int _num_cells_per_dim     # 2**(_over_refine_factor)

    # Children structure
    cdef int _density_factor        # this allows the tree to be built with more
                                    # than 8 children per node
    cdef int _num_children          # 2**(3*_density_factor)
    cdef int _num_children_per_dim  # 2**(_density_factor)

    # This is use for interpolation and is global for the Octree smoothing
    # operations
    cdef kernel_func kernel

    def __init__(self, double[:, ::1] &input_pos = None, left_edge = None,
                 right_edge = None, int n_ref=32, int over_refine_factor=1,
                 int density_factor=1, np.int64_t data_version=0,
                 int max_depth=20):

        # If this is the case, we are very likely just initialising an instance
        # and then going to load an existing Octree from disk, so we don't
        # really need to do anything
        if input_pos is None:
            return

        # These don't have setters as these would invalidate the tree which is
        # a feature we don't have
        # TODO: Add invalidation feature
        self._data_version = data_version
        self._n_ref = n_ref
        self._num_particles = input_pos.shape[0]
        self._max_depth = max_depth

        # Setting the properties which determines how children divide
        self.over_refine_factor = over_refine_factor
        self.density_factor = density_factor

        # Set up the initial idx of the particles, this keeps track of which
        # particle is which in the tree ordered array
        self._idx = np.arange(0, input_pos.shape[0], dtype=np.int64)

        # Set up the bounds and root node
        self.setup_bounds(input_pos, left_edge, right_edge)
        self.setup_root(input_pos)

        # Reserve some space for the nodes to be stored, this is a conversative
        # amount. If we exceed this the STL container will reallocate. This
        # will *not* invalidate any pointers
        # This decreases the conversative amount and keeps retrying, unless we
        # stil fail even with a small reserve, then we error
        cdef int exp_num_nodes = ((2**(3 * self.density_factor) *
                                  self._num_particles) // n_ref + 8)
        cdef int failed = 1
        while exp_num_nodes > 8 and failed == 1:
            try:
                reserve(&self.nodes, exp_num_nodes)
                failed = 0
            except MemoryError:
                exp_num_nodes = exp_num_nodes // 2
                failed = 1

        if failed == 1:
            raise MemoryError("Failed to allocate memory for octree!")

        # Now build the tree
        self.build_tree(&input_pos[0, 0])

        # Give up any excess reserved space
        # NOTE: this doubles the memory usage
        cdef vector[Node] temp
        cdef int i
        temp.reserve(self.nodes.size())
        for i in range(self.nodes.size()):
            temp.push_back(self.nodes[i])
        self.nodes.swap(temp)
        temp.clear()

    cdef int setup_bounds(self, double[:, ::1] &input_pos, left_edge=None,
                      right_edge=None):
        if left_edge is not None:
            for i in range(3):
                self._left_edge[i] = left_edge[i]
        else:
            for i in range(3):
                self._left_edge[i] = np.amin(input_pos[:,i])

        if right_edge is not None:
            for i in range(3):
                self._right_edge[i] = right_edge[i]
        else:
            for i in range(3):
                self._right_edge[i] = np.amax(input_pos[:,i])
        return 0

    cdef int setup_root(self, double[:, ::1] &input_pos):
        cdef Node root
        root.left_edge = self._left_edge
        root.right_edge = self._right_edge
        # Not strictly true and could lead to an issue later
        root.parent = 0

        # Always true in yt context
        root.start = 0
        root.end = input_pos.shape[0]*3

        root.children = 0
        root.leaf = 1
        root.depth = 0
        root.leaf_id = 0
        root.node_id = 0

        # Store the root in the nodes array
        self.nodes.push_back(root)
        return 0

    cdef int reset(self):
        # Clear the big memory users before we load an octree from disk
        self.nodes.clear()
        self._idx = np.zeros(1, dtype=np.int64)-1
        return 0

    cdef int build_tree(self, double * input_pos):
        # Generate an array to store the which particles are in each child oct
        cdef np.int64_t * split_arr
        split_arr = <np.int64_t*> malloc((self._num_children + 1) *
                                                            sizeof(np.int64_t))

        # Loop through the nodes in serial and process them, i.e, sort the
        # particles and create the children, which will increase the node.size
        # then we iterate through those children
        cdef int num_nodes_processed = 0
        while num_nodes_processed < self.nodes.size():
            self.process_node(&self.nodes[num_nodes_processed], input_pos,
                              split_arr)
            num_nodes_processed += 1

        free(split_arr)

        return 0

    cdef int process_node(self, Node * node, double * input_pos,
                          np.int64_t * split_arr) nogil:
        if(node.end - node.start <= 3*self._n_ref or
           node.depth > self._max_depth):
            return 0

        # Node is no longer a leaf
        node.leaf = 0

        # This sorts the children in the node and then stores the position of
        # the first and last particle in each node
        split_helper(node.start, node.end, self._num_children,
                     self._density_factor, split_arr, input_pos, &self._idx[0],
                     &node.left_edge[0], &node.right_edge[0])

        # Generate the node structures for the children, and store the position
        # of the first and last particles they contain
        self.generate_children(node, split_arr)

        return 0

    cdef inline void generate_children(self, Node * node,
                                       np.int64_t * split_arr) nogil:
        cdef int i, j, z, k, split_id
        cdef double dx, dy, dz
        cdef Node temp

        node.children = self.nodes.size()

        # Set the properties which are the same for all children
        temp.parent = node.node_id
        temp.leaf = 1
        temp.depth = node.depth + 1
        temp.leaf_id = 0
        temp.children = 0

        # Set up the values to be used to set the child boundaries
        dx = (node.right_edge[0] - node.left_edge[0]) / self._num_children_per_dim
        dy = (node.right_edge[1] - node.left_edge[1]) / self._num_children_per_dim
        dz = (node.right_edge[2] - node.left_edge[2]) / self._num_children_per_dim

        # Loop through and append the children setting the node dependent values
        z = node.children
        split_id = 0
        for i in range(self._num_children_per_dim):
            for j in range(self._num_children_per_dim):
                for k in range(self._num_children_per_dim):
                    temp.left_edge[0] = node.left_edge[0] + i*dx
                    temp.left_edge[1] = node.left_edge[1] + j*dy
                    temp.left_edge[2] = node.left_edge[2] + k*dz
                    temp.right_edge[0] = node.left_edge[0] + (i+1)*dx
                    temp.right_edge[1] = node.left_edge[1] + (j+1)*dy
                    temp.right_edge[2] = node.left_edge[2] + (k+1)*dz
                    temp.node_id = z

                    temp.start = split_arr[split_id]
                    temp.end = split_arr[split_id + 1]

                    self.nodes.push_back(temp)
                    z+=1
                    split_id+=1

    @property
    def size_bytes(self):
        return sizeof(Node) * self.nodes.size()

    @property
    def max_depth(self):
        return self._max_depth

    @property
    def data_version(self):
        return self._data_version

    @property
    def num_particles(self):
        return self._num_particles

    @property
    def n_ref(self):
        return self._n_ref

    @property
    def num_octs(self):
        return self.nodes.size()

    @property
    def over_refine_factor(self):
        return self._over_refine_factor

    @property
    def idx(self):
        return np.asarray(self._idx)

    @over_refine_factor.setter
    def over_refine_factor(self, value):
        self._over_refine_factor = value
        self._num_cells = 2**(3 * value)
        self._num_cells_per_dim = 2**value

    @property
    def density_factor(self):
        return self._density_factor

    @density_factor.setter
    def density_factor(self, value):
        self._density_factor = value
        self._num_children = 2**(3 * value)
        self._num_children_per_dim  = 2**value

    @property
    def cell_positions(self):
        cdef int i, j, z, k, l, num_leaves

        # Find all the leaves
        num_leaves = 0
        for i in range(self.nodes.size()):
            if self.nodes[i].leaf == 1:
                num_leaves += 1

        cdef np.float64_t[:, :] pos = np.zeros((num_leaves*self._num_cells, 3),
                                               dtype='float64')
        cdef double leftx, lefty, leftz, rightx, righty, rightz
        z = 0
        for i in range(self.nodes.size()):
            if self.nodes[i].leaf == 0:
                continue

            leftx = self.nodes[i].left_edge[0]
            lefty = self.nodes[i].left_edge[1]
            leftz = self.nodes[i].left_edge[2]
            rightx = self.nodes[i].right_edge[0]
            righty = self.nodes[i].right_edge[1]
            rightz = self.nodes[i].right_edge[2]

            self.nodes[i].leaf_id = z

            # we have to generate cell locations
            dx = (rightx - leftx) / self._num_cells_per_dim
            dy = (righty - lefty) / self._num_cells_per_dim
            dz = (rightz - leftz) / self._num_cells_per_dim

            for j in range(self._num_cells_per_dim):
                for k in range(self._num_cells_per_dim):
                    for l in range(self._num_cells_per_dim):
                        pos[z, 0] = leftx + (j + 0.5) * dx
                        pos[z, 1] = lefty + (k + 0.5) * dy
                        pos[z, 2] = leftz + (l + 0.5) * dz
                        z+=1

        return np.asarray(pos)

    # TODO: move these to the location of the rest of the deposition operations
    cdef void smooth_onto_cells(self, np.float64_t[:] buff,
                                np.float64_t[:] buff_den, np.float64_t posx,
                                np.float64_t posy, np.float64_t posz,
                                np.float64_t hsml, np.float64_t prefactor,
                                np.float64_t prefactor_norm, Node * node,
                                int use_normalization=0):

        cdef Node * child
        cdef int i, j, k, l
        cdef double q_ij, diff_x, diff_y, diff_z, dx, dy, dz, voxel_hsml2
        cdef double leftx, lefty, leftz, rightx, righty, rightz

        # If not a leaf, then check if particle is in the children and go check
        # through those. This is recursive - currently
        if node.leaf == 0:
            child = &self.nodes[node.children]
            for i in range(self._num_children):
                leftx = child[i].left_edge[0]
                lefty = child[i].left_edge[1]
                leftz = child[i].left_edge[2]
                rightx = child[i].right_edge[0]
                righty = child[i].right_edge[1]
                rightz = child[i].right_edge[2]

                if leftx - posx < hsml and posx - rightx < hsml:
                    if lefty - posy < hsml and posy - righty < hsml:
                        if leftz - posz < hsml and posz - rightz < hsml:
                            self.smooth_onto_cells(buff, buff_den, posx, posy,
                                                   posz, hsml, prefactor,
                                                   prefactor_norm, &child[i],
                                                   use_normalization)
        else:
            leftx = node.left_edge[0]
            lefty = node.left_edge[1]
            leftz = node.left_edge[2]
            rightx = node.right_edge[0]
            righty = node.right_edge[1]
            rightz = node.right_edge[2]

            # We have to generate cell locations
            dx = (rightx - leftx) / self._num_cells_per_dim
            dy = (righty - lefty) / self._num_cells_per_dim
            dz = (rightz - leftz) / self._num_cells_per_dim

            # Loop through each cell and calculate the contribution, l is the
            # number of the cell we are in
            l = 0
            for i in range(self._num_cells_per_dim):
                for j in range(self._num_cells_per_dim):
                    for k in range(self._num_cells_per_dim):
                        diff_x = (leftx + (i + 0.5) * dx - posx)
                        diff_x *= diff_x
                        diff_y = (lefty + (j + 0.5) * dy - posy)
                        diff_y *= diff_y
                        diff_z = (leftz + (k + 0.5) * dz - posz)
                        diff_z *= diff_z

                        voxel_hsml2 = hsml*hsml
                        q_ij = math.sqrt((diff_x + diff_y + diff_z) /
                                                                    voxel_hsml2)

                        if use_normalization:
                            buff_den[node.leaf_id + l] += (prefactor_norm *
                                                           self.kernel(q_ij))
                        buff[node.leaf_id + l] += prefactor * self.kernel(q_ij)
                        l += 1

    def interpolate_sph_cells(self, np.float64_t[:] buff,
                              np.float64_t[:] buff_den, np.float64_t[:] posx,
                              np.float64_t[:] posy, np.float64_t[:] posz,
                              np.float64_t[:] pmass, np.float64_t[:] pdens,
                              np.float64_t[:] hsml, np.float64_t[:] field,
                              kernel_name="cubic", int use_normalization=0):

        self.kernel = get_kernel_func(kernel_name)

        cdef int i
        cdef double prefactor, prefactor_norm

        for i in range(posx.shape[0]):
            prefactor = pmass[i] / pdens[i] / hsml[i]**3
            prefactor_norm = prefactor
            prefactor *= field[i]

            self.smooth_onto_cells(buff, buff_den, posx[i], posy[i], posz[i],
                                   hsml[i], prefactor, prefactor_norm,
                                   &self.nodes[0],
                                   use_normalization=use_normalization)

    def __richcmp__(self, CyOctree other, op):
        if op == 2:
            return self._is_equal(other)
        elif op ==3:
            return not self._is_equal(other)
        else:
            raise NotImplementedError(("Use == or !=, other comparisons have " +
                                       "not been implemented."))

    def _is_equal(self, CyOctree other):
        cdef bool same = True

        for i in range(3):
            if self._left_edge[i] != other._left_edge[i]:
                same = False
            if self._right_edge[i] != other._right_edge[i]:
                same = False

        if self._n_ref != other._n_ref:
            same = False

        if self._over_refine_factor != other._over_refine_factor:
            same = False

        if self._density_factor != other._density_factor:
            same = False

        if self._data_version != other._data_version:
            same = False

        return same

    # TODO: this code is much slower than I would like, this is likely due to
    # the use of struct -> plan to replace this. A c++ approach is probably
    # faster and more intuitive
    def save(self, fname = None):
        if fname is None:
            raise ValueError("A filename must be specified to save the octree!")

        with open(fname,'wb') as f:
            f.write(struct.pack('2Q3iq', self._num_particles, self.num_octs,
                                      self._n_ref, self.over_refine_factor,
                                      self.density_factor, self._data_version))
            f.write(struct.pack('{}Q'.format(self.num_particles),
                                *self.idx))

            for i in range(self.num_octs):
                f.write(struct.pack('6d4Q?2QB',
                                    self.nodes[i].left_edge[0],
                                    self.nodes[i].left_edge[1],
                                    self.nodes[i].left_edge[2],
                                    self.nodes[i].right_edge[0],
                                    self.nodes[i].right_edge[1],
                                    self.nodes[i].right_edge[2],
                                    self.nodes[i].start,
                                    self.nodes[i].end,
                                    self.nodes[i].parent,
                                    self.nodes[i].children,
                                    self.nodes[i].leaf,
                                    self.nodes[i].node_id,
                                    self.nodes[i].leaf_id,
                                    self.nodes[i].depth))

    def load(self, fname = None):
        if fname is None:
            raise ValueError("A filename must be specified to load the octtree!")
        # clear any current tree we have loaded
        self.reset()

        cdef Node temp
        cdef np.int64_t num_octs
        with open(fname,'rb') as f:
            (self._num_particles, num_octs, self._n_ref,
             self.over_refine_factor, self.density_factor,
             self._data_version) = \
                struct.unpack('2Q3iq', f.read(40))
            self._idx = \
                np.asarray(struct.unpack('{}Q'.format(self.num_particles),
                           f.read(8*self.num_particles)), dtype=np.int64)

            reserve(&self.nodes, num_octs+1)
            for i in range(num_octs):
                (temp.left_edge[0], temp.left_edge[1], temp.left_edge[2],
                 temp.right_edge[0], temp.right_edge[1], temp.right_edge[2],
                 temp.start, temp.end, temp.parent, temp.children, temp.leaf,
                 temp.node_id, temp.leaf_id, temp.depth) = \
                struct.unpack('6d4Q?2QB', f.read(105))
                self.nodes.push_back(temp)

        for i in range(3):
            self._left_edge[i] = self.nodes[0].left_edge[i]
            self._right_edge[i] = self.nodes[0].right_edge[i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int split_helper(int start, int end, int num_children, int max_splits,
                       np.int64_t * split_arr, np.float64_t * pos,
                       np.int64_t * idx, np.float64_t * left,
                       np.float64_t * right) nogil except -1:
    '''
    This takes in the split array and sets up the first and last particle, it
    then calls the spit function.
    '''
    cdef double lower = left[0]
    cdef double upper = right[0]

    split_arr[0] = start
    split_arr[num_children] = end

    split(0, num_children, max_splits, 0, 0, split_arr, idx, pos, left, right,
          lower, upper)

    return 0

cdef int reserve(vector[Node] * vec, int amount) except +MemoryError:
    '''
    This attempts to reserve memory and propagates any errors back.
    '''
    vec.reserve(amount)
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate(np.float64_t * array, np.int64_t * idx, int offset,  double value,
                  np.int64_t start, np.int64_t end) nogil except -1:
    '''
    This takes in an axis, a value and the particles position array.

    It splits the array so all values in the positions array with positions in
    the axis dimension less than value are on the left, and those with a value
    greater are on the right.
    '''
    cdef np.int64_t index
    cdef np.int64_t idx_index = start // 3, idx_split = idx_index
    cdef np.int64_t split = start

    for index in range(start, end, 3):
        idx_index += 1
        if array[index + offset] < value:
            idx[idx_split], idx[idx_index] = idx[idx_index], idx[idx_split]
            array[split], array[index] = array[index], array[split]
            array[split+1], array[index+1] = array[index+1], array[split+1]
            array[split+2], array[index+2] = array[index+2], array[split+2]
            split+=3
            idx_split+=1

    return split

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int split(int start, int end, int max_splits, int splits, int axis,
           np.int64_t * split_arr, np.int64_t * idx, np.float64_t * pos,
           np.float64_t * left_edge, np.float64_t * right_edge,
           np.float64_t lower, np.float64_t upper) nogil except -1:
    '''
    This splits the arrays recursively such that we split the particles in the
    x, y and z axis such that we get 2**(density_factor) children. We store
    which particles are in each children in the split_arr
    '''
    cdef int mid
    cdef double temp_value

    if splits == max_splits and axis < 2:
        splits = 0
        axis += 1
        lower = left_edge[axis]
        upper = right_edge[axis]

    if splits < max_splits:
        splits += 1
        mid = (start + end) // 2
        temp_value = (upper + lower) / 2

        split_arr[mid] = seperate(pos, idx, axis, temp_value, split_arr[start],
                                  split_arr[end])

        split(start, mid, max_splits, splits, axis, split_arr, idx, pos,
              left_edge, right_edge, lower, temp_value)
        split(mid, end, max_splits, splits, axis, split_arr, idx, pos,
              left_edge, right_edge, temp_value, upper)

    return 0
