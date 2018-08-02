cimport numpy as np
import numpy as np
cimport cython

# for writing to file
import struct

# cpp includes
from libcpp.vector cimport vector
from libcpp cimport bool

# c includes
cimport libc.math as math

from yt.geometry.particle_deposit cimport \
    kernel_func, get_kernel_func

cdef struct Node:
    double left_edge[3]
    double right_edge[3]        # may be more efficient to store the width instead
    unsigned long start         # which particles we store
    unsigned long end
    unsigned long parent        # position of parent in Octree.nodes
    unsigned long children      # position of 0th child, children are contiguous
    bool leaf
    unsigned long node_id       # these probably should be longs
    unsigned long leaf_id       # not actually sure how useful this is
    unsigned int depth

@cython.boundscheck(True)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class PyOctree:
    cdef int _state             # 0 if tree is not built, 1 for a built tree
    cdef double _left_edge[3]   # boundary conditions for the oct
    cdef double _right_edge[3]
    cdef long[:] _idx           # ordering of particles used by the tree
    cdef int _n_ref
    cdef int _num_octs
    cdef int _num_particles
    cdef int _max_depth

    cdef vector[Node] nodes     # This is an STL container to store the octs

    # this is use for interpolation and is global for the octree smoothing
    # operations
    cdef kernel_func kernel

    def __init__(self, double[:, ::1] &input_pos = None, left_edge = None,
                 right_edge = None, int n_ref=32):
        # if this is the case, we are very likely just initialising an instance
        # and then going to load an existing Octree from memory, so we don't
        # really need to do anything
        if input_pos is None:
            return

        # max number of particles per oct
        self.n_ref = n_ref
        self.num_particles = input_pos.shape[0]

        # set up the initial idx of the particles
        self.idx = np.arange(0, input_pos.shape[0], dtype=np.int64)

        # set up the bounds
        self.setup_bounds(input_pos, left_edge, right_edge)
        self.setup_root(input_pos)

        # replace this with something more conservative -> probably allow STL to
        # manage its own reallocation
        self.nodes.reserve(10000000)

        # now build the tree, we do this with no gil as we will attempt to
        # parallelize in the future
        with nogil:
            self.process_node(&self.nodes[0], &input_pos[0, 0])

        # setup the final parameters
        self.num_octs = self.nodes.size()

    cdef setup_bounds(self, double[:, ::1] &input_pos, left_edge=None,
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

    cdef setup_root(self, double[:, ::1] &input_pos):
        cdef Node root

        root.left_edge = self._left_edge
        root.right_edge = self._right_edge
        root.parent = 0

        root.start = 0
        root.end = input_pos.shape[0]*3

        root.children = 0
        root.leaf = 1
        root.depth = 0
        root.leaf_id = 0
        root.node_id = 0

        # store the root in an array and make a convenient pointer
        self.nodes.push_back(root)

    cdef reset(self):
        # clear the big memory users
        self.nodes.clear()
        self._idx = np.zeros(1, dtype=int)-1

        # reset some general parameters
        for i in range(3):
            self._left_edge[i] = 0.0
            self._right_edge[i] = 0.0

        self.n_ref = 0
        self.num_octs = 0
        self.num_particles = 0

    # TODO: this this will be a lot more "re-allocation" safe to pass the position
    # of the node in the vector rather than a memory address which can be
    # invalidated upon a reallocation
    cdef int process_node(self, Node * node, double * input_pos) nogil:
        node.leaf = 0

        cdef int i
        cdef int splits[9]
        cdef double temp_value

        # make the children placeholders
        self.generate_children(node)

        # TODO: refactor this section to use an over-refine-factor
        # set up the split integers
        splits[0] = node.start
        splits[8] = node.end

        # split into two partitions based on x value
        temp_value = (node.left_edge[0] + node.right_edge[0]) / 2
        splits[4] = seperate(input_pos, &self._idx[0], 0,
                             temp_value,
                             splits[0], splits[8])

        # split into four partitions using the y value
        temp_value = (node.left_edge[1] + node.right_edge[1]) / 2
        for i in range(0, 2, 1):
            splits[2 + i*4] = seperate(input_pos, &self._idx[0], 1,
                                       temp_value,
                                       splits[4*i], splits[4*i + 4])

        # split into eight partitions using the z value
        temp_value = (node.left_edge[2] + node.right_edge[2]) / 2
        for i in range(0, 4, 1):
            splits[2*i + 1] = seperate(input_pos, &self._idx[0], 2,
                                       temp_value,
                                       splits[2*i], splits[2*i + 2])

        cdef Node * child = &(self.nodes[node.children])
        for i in range(0, 8, 1):
            child[i].start = splits[i]
            child[i].end =  splits[i+1]

        # stop here if no children
        if(node.end - node.start <= 3*self._n_ref):
            return 0

        for i in range(0, 8, 1):
            self.process_node(&child[i], input_pos)

        return 0

    cdef inline void generate_children(self, Node * node) nogil:
        cdef int i, j, z, k
        cdef double dx, dy, dz
        cdef Node temp

        node.children = self.nodes.size()

        # set the properties which are the same for all children
        temp.parent = node.node_id
        temp.leaf = 1
        temp.depth = node.depth + 1
        temp.leaf_id = 0

        dx = (node.right_edge[0] - node.left_edge[0]) / 2
        dy = (node.right_edge[1] - node.left_edge[1]) / 2
        dz = (node.right_edge[2] - node.left_edge[2]) / 2

        z = node.children
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    temp.left_edge[0] = node.left_edge[0] + i*dx
                    temp.left_edge[1] = node.left_edge[1] + j*dy
                    temp.left_edge[2] = node.left_edge[2] + k*dz
                    temp.right_edge[0] = node.left_edge[0] + (i+1)*dx
                    temp.right_edge[1] = node.left_edge[1] + (j+1)*dy
                    temp.right_edge[2] = node.left_edge[2] + (k+1)*dz

                    temp.node_id = z

                    self.nodes.push_back(temp)
                    z+=1

    @property
    def size_bytes(self):
        return sizeof(Node) * self.nodes.size()

    # TODO: fix this implementation
    @property
    def max_depth(self):
        return self.max_depth

    # TODO: replace this with a descriptor
    @property
    def idx(self):
        return np.asarray(self._idx)

    @idx.setter
    def idx(self, array):
        self._idx = array

    @property
    def num_octs(self):
        return self._num_octs

    @num_octs.setter
    def num_octs(self, val):
        self._num_octs = val

    @property
    def n_ref(self):
        return self._n_ref

    @n_ref.setter
    def n_ref(self, val):
        self._n_ref = val

    @property
    def num_particles(self):
        return self._num_particles

    @num_particles.setter
    def num_particles(self, val):
        self._num_particles = val

    # TODO: re-write this to pass through in a much more efficient manner,
    # currently dreadfully slow
    @property
    def leaf_positions(self):
        cdef int i, j, z = 0
        positions = []
        for i in range(self.nodes.size()):
            if self.nodes[i].leaf == 0:
                continue
            else:
                for j in range(3):
                    positions.append((self.nodes[i].left_edge[j] +
                                      self.nodes[i].right_edge[j]) / 2)

                self.nodes[i].leaf_id = z
                z+=1

        positions = np.asarray(positions)
        return positions.reshape((-1,3))

    # TODO: this code is much slower than I would like, this is likely due to
    # the use of struct -> plan to replace this. A c++ approach is probably
    # faster and more intuitive
    def save(self, fname = None, tree_hash = 0):
        if fname is None:
            raise ValueError("A filename must be specified to save the octree!")

        # TODO: we need to save the tree hash as well
        with open(fname,'wb') as f:
            f.write(struct.pack('3d', self.num_particles, self.num_octs,
                                      self.n_ref))
            f.write(struct.pack('{}d'.format(self.num_particles),
                                *self.idx))
            for i in range(self.num_octs):
                f.write(struct.pack('10d?2di',
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
        with open(fname,'rb') as f:
            (self.num_particles, self.num_octs, self.n_ref) = \
                struct.unpack('3d', f.read(24))
            self.idx = \
                np.asarray(struct.unpack('{}d'.format(self.num_particles),
                           f.read(8*self.num_particles)), dtype=np.int64)
            reserve(&self.nodes, self.num_octs+1)
            for i in range(self.num_octs):
                (temp.left_edge[0], temp.left_edge[1], temp.left_edge[2],
                 temp.right_edge[0], temp.right_edge[1], temp.right_edge[2],
                 temp.start, temp.end, temp.parent, temp.children, temp.leaf,
                 temp.node_id, temp.leaf_id, temp.depth) = \
                struct.unpack('10d?2di', f.read(108))

                self.nodes.push_back(temp)

    cdef void smooth_onto_leaves(self, np.float64_t[:] buff, np.float64_t posx,
                                 np.float64_t posy, np.float64_t posz,
                                 np.float64_t hsml, np.float64_t prefactor,
                                 Node * node):

        cdef Node * child
        cdef double q_ij, diff_X, diff_y, diff_z

        if node.leaf == 0:
            child = &self.nodes[node.children]
            for i in range(8):
                if child[i].left_edge[0] - posx < hsml and \
                        posx - child[i].right_edge[0] < hsml:
                    if child[i].left_edge[1] - posy < hsml and \
                            posy - child[i].right_edge[1] < hsml:
                        if child[i].left_edge[2] - posz < hsml and \
                                posz - child[i].right_edge[2] < hsml:
                            self.smooth_onto_leaves(buff, posx, posy, posz,
                                                    hsml, prefactor, &child[i])

        else:
            diff_x = ((node.left_edge[0] + node.right_edge[0]) / 2 - posx)
            diff_x *= diff_x
            diff_y = ((node.left_edge[1] + node.right_edge[1]) / 2 - posy)
            diff_y *= diff_y
            diff_z = ((node.left_edge[2] + node.right_edge[2]) / 2 - posz)
            diff_z *= diff_z

            q_ij = math.sqrt((diff_x + diff_y + diff_z) / (hsml*hsml))

            buff[node.leaf_id] += prefactor * self.kernel(q_ij)

    def interpolate_sph_octs(self, np.float64_t[:] buff, np.float64_t[:] posx,
                             np.float64_t[:] posy, np.float64_t[:] posz,
                             np.float64_t[:] pmass, np.float64_t[:] pdens,
                             np.float64_t[:] hsml, np.float64_t[:] field,
                             kernel_name="cubic"):

        self.kernel = get_kernel_func(kernel_name)

        cdef int i
        cdef double prefactor

        for i in range(posx.shape[0]):
            prefactor = pmass[i] / pdens[i] / hsml[i]**3
            prefactor *= field[i]

            self.smooth_onto_leaves(buff, posx[i], posy[i], posz[i], hsml[i],
                                    prefactor, &self.nodes[0])

# This is a utility function which is able to catch errors when we attempt to
# reserve too much memory
cdef int reserve(vector[Node] * vec, int amount) except +MemoryError:
    vec.reserve(amount)
    return 0

# This is a utility function to separate an array into a region smaller than the
# splitting value and a region larger than the splitting value, this is used to
# determine which particles occupy which octs
@cython.boundscheck(True)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate(double * array, long * idx, char offset,  double value,
                  long start, long end) nogil except -1:
    cdef long index
    cdef long idx_index = start/3, idx_split = idx_index
    cdef long split = start

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
