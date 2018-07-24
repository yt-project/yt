cimport numpy as np
import numpy as np
import struct
cimport cython
from libcpp.vector cimport vector

# this is the underlying c structure of a node
cdef struct Node:
    double left_edge[3]
    double right_edge[3]
    int start
    int end
    int parent
    int children
    int leaf
    int node_id

# this is the underlying c structure for the octree, it is very basic in terms
# of the values that it stores
cdef struct Octree:
    vector[Node] nodes
    double left_edge[3]
    double right_edge[3]
    int n_ref
    Node * root
    int * idx

cdef class PyOctree:
    cdef Octree c_tree
    cdef int[:] _idx
    cdef int _n_ref
    cdef int _num_octs
    cdef int _num_particles

    def __init__(self, double[:, ::1] &input_pos = None, left_edge = None,
                 right_edge = None, int n_ref=32):
        # if this is the case, we are very likely just initialising an instance
        # and then going to load an existing kdtree from memory
        if input_pos is None:
            return

        self.setup_bounds(input_pos, left_edge, right_edge)
        self.setup_ctree(input_pos, n_ref)
        self.setup_root(input_pos)

        # now build the tree, we do this with no gil as we will attempt to
        # parallelize in the future
        with nogil:
            process_node(self.c_tree, self.c_tree.root, &(input_pos[0, 0]))

        # setup the final parameters
        self.n_ref = n_ref
        self.num_octs = self.c_tree.nodes.size()
        self.num_particles = input_pos.shape[0]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef setup_ctree(self, double[:, ::1] &input_pos, int n_ref):
        self.n_ref = n_ref
        self.c_tree.n_ref = n_ref

        self.idx = np.arange(0, input_pos.shape[0], dtype=np.int32)
        self.c_tree.idx = &self._idx[0]

        reserve(&self.c_tree.nodes, 10000000)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef setup_bounds(self, double[:, ::1] &input_pos, left_edge=None,
                     right_edge=None):
        if left_edge is not None:
            for i in range(3):
                self.c_tree.left_edge[i] = left_edge[i]
        else:
            for i in range(3):
                self.c_tree.left_edge[i] = np.amin(input_pos[:,i])

        if right_edge is not None:
            for i in range(3):
                self.c_tree.right_edge[i] = right_edge[i]
        else:
            for i in range(3):
                self.c_tree.right_edge[i] = np.amax(input_pos[:,i])

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef setup_root(self, double[:, ::1] &input_pos):
        cdef Node root
        root.left_edge = self.c_tree.left_edge
        root.right_edge = self.c_tree.right_edge
        root.parent = -1
        root.start = 0
        root.end = input_pos.shape[0]*3
        root.node_id = 0

        # store the root in an array and make a convenient pointer
        self.c_tree.nodes.push_back(root)
        self.c_tree.root = &(self.c_tree.nodes[0])

    def reset(self):
        # reset the c tree
        self.c_tree.nodes.clear()
        for i in range(3):
            self.c_tree.left_edge[i] = 0.0
            self.c_tree.right_edge[i] = 0.0
        self.c_tree.n_ref = 0

        # reset python properties
        self.n_ref = 0
        self.num_octs = 0
        self.idx = np.array([0], dtype=np.int32)
        self.num_particles = 0

    @property
    def idx(self):
        return np.asarray(self._idx)

    @idx.setter
    def idx(self, array):
        self._idx = array
        self.c_tree.idx = &self._idx[0]

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

    # TODO: this code is much slower than I would like, this is likely due to
    # the use of struct -> plan to replace this
    def save(self, fname = None, tree_hash = 0):
        if fname is None:
            raise ValueError("A filename must be specified to save the kdtree!")

        # TODO: we need to save the tree hash as well
        with open(fname,'wb') as f:
            f.write(struct.pack('3i', self.num_particles, self.num_octs,
                                      self.n_ref))
            f.write(struct.pack('{}i'.format(self.num_particles),
                                *self.idx))
            for i in range(self.num_octs):
                f.write(struct.pack('6d5i',
                                    self.c_tree.nodes[i].left_edge[0],
                                    self.c_tree.nodes[i].left_edge[1],
                                    self.c_tree.nodes[i].left_edge[2],
                                    self.c_tree.nodes[i].right_edge[0],
                                    self.c_tree.nodes[i].right_edge[1],
                                    self.c_tree.nodes[i].right_edge[2],
                                    self.c_tree.nodes[i].start,
                                    self.c_tree.nodes[i].end,
                                    self.c_tree.nodes[i].parent,
                                    self.c_tree.nodes[i].children,
                                    self.c_tree.nodes[i].leaf))

    def load(self, fname = None):
        if fname is None:
            raise ValueError("A filename must be specified to load the kdtree!")

        self.reset()

        cdef Node temp
        with open(fname,'rb') as f:
            (self.num_particles, self.num_octs, self.n_ref) = \
                struct.unpack('3i', f.read(12))
            self.idx = \
                np.asarray(struct.unpack('{}i'.format(self.num_particles),
                           f.read(4*self.num_particles)), dtype=np.int32)
            reserve(&self.c_tree.nodes, self.num_octs+1)
            for i in range(self.num_octs):
                (temp.left_edge[0], temp.left_edge[1], temp.left_edge[2],
                 temp.right_edge[0], temp.right_edge[1], temp.right_edge[2],
                 temp.start, temp.end, temp.parent, temp.children,
                 temp.leaf) = struct.unpack('6d5i', f.read(68))

                self.c_tree.nodes.push_back(temp)

# this is a utility function which is able to catch errors when we attempt to
# reserve too much memory
cdef int reserve(vector[Node] * vec, int amount) except +MemoryError:
    vec.reserve(amount)
    return 0

# this makes the children and stores them in the tree.node
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int generate_children(Octree &tree, Node * node) nogil:
    cdef int i, j, z, k
    cdef double dx, dy, dz
    cdef Node temp

    node.children = tree.nodes.size()

    temp.parent = node.node_id

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
                tree.nodes.push_back(temp)
                z+=1
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int process_node(Octree &tree, Node * node, double * input_pos) nogil:
    # skip if not enough children
    if(node.end - node.start <= 3*tree.n_ref):
        node.leaf = 1
        return 0

    cdef int i
    cdef int splits[9]
    cdef double temp_value

    # make the children placeholders
    generate_children(tree, node)

    # TODO: refactor this section to use an over-refine-factor
    # set up the split integers
    splits[0] = node.start
    splits[8] = node.end

    # split into two partitions based on x value
    temp_value = (node.left_edge[0] + node.right_edge[0]) / 2
    splits[4] = seperate(input_pos, &tree.idx[0], 0,
                         temp_value,
                         splits[0], splits[8])

    # split into four partitions using the y value
    temp_value = (node.left_edge[1] + node.right_edge[1]) / 2
    for i in range(0, 2, 1):
        splits[2 + i*4] = seperate(input_pos, &tree.idx[0], 1,
                                   temp_value,
                                   splits[4*i], splits[4*i + 4])

    # split into eight partitions using the z value
    temp_value = (node.left_edge[2] + node.right_edge[2]) / 2
    for i in range(0, 4, 1):
        splits[2*i + 1] = seperate(input_pos, &tree.idx[0], 2,
                                   temp_value,
                                   splits[2*i], splits[2*i + 2])

    # now just need to tell the children which particles they store and repeat
    # the process
    cdef Node * child = &(tree.nodes[node.children])
    for i in range(0, 8, 1):
        child[i].start = splits[i]
        child[i].end =  splits[i+1]
        process_node(tree, &child[i], input_pos)

    return 0

# this is a utility function to separate an array into a region smaller than the
# splitting value and a region larger than the splitting value
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate(double * array, int * idx, int offset,  double &value,
                  int &start, int &end) nogil:
    cdef int index
    cdef int idx_index = start/3, idx_split = idx_index
    cdef int split = start

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
