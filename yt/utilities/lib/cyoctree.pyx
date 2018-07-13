cimport numpy as np
import numpy as np
cimport cython
from libcpp.vector cimport vector

# this is the underlying c structure of a node
cdef struct Node:
    double left_edge[3]
    double right_edge[3]
    int start
    int end
    int parent

# this is the underlying c structure for the octree, it is very basic in terms
# of the values that it stores
cdef struct Octree:
    vector[Node] nodes
    double left_edge[3]
    double right_edge[3]
    int n_ref
    int num_nodes
    Node * root
    vector[int] ids

# this is the python interface to the cython octree
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class PyOctree:
    cdef Octree c_tree

    def __init__(self, double[:, ::1] &input_pos, left_edge = None,
                 right_edge = None, int n_ref=32):
        # the max number of particles per leaf
        self.c_tree.n_ref = n_ref

        # set up the boundaries of the tree
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

        # set up the initial node storage
        self.c_tree.num_nodes = 1

        # attempt to reserve some memory, we don't want to reallocate
        self.c_tree.nodes.reserve(1000000)

        # set up the initial tree id list
        self.c_tree.ids.reserve(input_pos.shape[0]+1)
        for i in range(input_pos.shape[0]):
            self.c_tree.ids.push_back(i)

        # make the root node
        self.setup_root(input_pos.shape[0]*3)

        # now build the tree
        with nogil:
            process_node(self.c_tree, self.c_tree.root, &(input_pos[0, 0]))

        # TODO:
        # shrink the nodes capacity to fit

    def setup_root(self, end):
        cdef Node root
        root.left_edge = self.c_tree.left_edge
        root.right_edge = self.c_tree.right_edge
        root.parent = -1
        root.start = 0
        root.end = end

        # store the root in an array and make a convenient pointer
        self.c_tree.nodes.push_back(root)
        self.c_tree.root = &(self.c_tree.nodes[0])

    @property
    def idx(self):
        return np.asarray(self.c_tree.ids)

    @property
    def num_octs(self):
        return self.c_tree.nodes.size()

# this makes the children and stores them in the tree.node
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int generate_children(Octree &tree, Node * node, int &parent_id) nogil:
    cdef int i, j, z, k
    cdef double dx, dy, dz
    cdef Node * child = &(tree.nodes[parent_id + 1])

    dx = (node.right_edge[0] - node.left_edge[0]) / 2
    dy = (node.right_edge[1] - node.left_edge[1]) / 2
    dz = (node.right_edge[2] - node.left_edge[2]) / 2

    z = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                tree.nodes.push_back(tree.nodes[parent_id])
                child[z].left_edge[0] = node.left_edge[0] + i*dx
                child[z].left_edge[1] = node.left_edge[1] + j*dy
                child[z].left_edge[2] = node.left_edge[2] + k*dz
                child[z].right_edge[0] = node.left_edge[0] + (i+1)*dx
                child[z].right_edge[1] = node.left_edge[1] + (j+1)*dy
                child[z].right_edge[2] = node.left_edge[2] + (k+1)*dz
                child[z].parent  = parent_id
                tree.num_nodes += 1
                z += 1
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int process_node(Octree &tree, Node * node, double * input_pos) nogil:
    # skip if not enough children
    if(node.end - node.start <= 3*tree.n_ref):
        return 0

    cdef int parent_id = tree.num_nodes - 1, i
    cdef int splits[9]
    cdef double temp_value

    # make the children placeholders
    generate_children(tree, node, parent_id)

    # set up the split integers
    splits[0] = node.start
    splits[8] = node.end

    # split into two partitions based on x value
    temp_value = (node.left_edge[0] + node.right_edge[0]) / 2
    splits[4] = seperate(input_pos, &tree.ids[0], 0,
                         temp_value,
                         splits[0], splits[8])

    # split into four partitions using the y value
    temp_value = (node.left_edge[1] + node.right_edge[1]) / 2
    for i in range(0, 2, 1):
        splits[2 + i*4] = seperate(input_pos, &tree.ids[0], 1,
                                   temp_value,
                                   splits[4*i], splits[4*i + 4])

    # split into eight partitions using the z value
    temp_value = (node.left_edge[2] + node.right_edge[2]) / 2
    for i in range(0, 4, 1):
        splits[2*i + 1] = seperate(input_pos, &tree.ids[0], 2,
                                   temp_value,
                                   splits[2*i], splits[2*i + 2])

    # now just need to tell the children which particles they store
    cdef Node * child = &(tree.nodes[parent_id + 1])

    for i in range(0, 8, 1):
        child[i].start = splits[i]
        child[i].end =  splits[i+1]

    for i in range(0, 8, 1):
        process_node(tree, &tree.nodes[parent_id + 1 + i], input_pos)

    return 0

# this is a utility function to separate an array into a region smaller than the
# splitting value and a region larger than the splitting value
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int seperate(double * array, int * ids, int offset,  double &value,
                  int &start, int &end) nogil:
    cdef int index
    cdef int ids_index = start/3, ids_split = ids_index
    cdef int split = start

    for index in range(start, end, 3):
        ids_index += 1
        if array[index + offset] < value:
            ids[ids_split], ids[ids_index] = ids[ids_index], ids[ids_split]
            array[split], array[index] = array[index], array[split]
            array[split+1], array[index+1] = array[index+1], array[split+1]
            array[split+2], array[index+2] = array[index+2], array[split+2]
            split+=3
            ids_split+=1

    return split
