import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free, abs
from fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip
from field_interpolation_tables cimport \
    FieldInterpolationTable, FIT_initialize_table, FIT_eval_transfer,\
    FIT_eval_transfer_with_light
from fixed_interpolator cimport *

from cython.parallel import prange, parallel, threadid

cdef extern from "stdlib.h":
    # NOTE that size_t might not be int
    void *alloca(int)


DEF Nch = 4

cdef struct Split:
    int dim
    np.float64_t pos

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class Node:

    cdef public Node left
    cdef public Node right
    cdef public Node parent
    cdef public int grid
    cdef public int node_id
    cdef np.float64_t * left_edge
    cdef np.float64_t * right_edge
    cdef public data
    cdef Split * split

    def __cinit__(self, 
                  Node parent, 
                  Node left, 
                  Node right, 
                  np.ndarray[np.float64_t, ndim=1] left_edge,
                  np.ndarray[np.float64_t, ndim=1] right_edge,
                  int grid,
                  int node_id):
        self.left = left
        self.right = right
        self.parent = parent
        cdef int i
        self.left_edge = <np.float64_t *> malloc(sizeof(np.float64_t) * 3)
        self.right_edge = <np.float64_t *> malloc(sizeof(np.float64_t) * 3)
        for i in range(3):
            self.left_edge[i] = left_edge[i]
            self.right_edge[i] = right_edge[i]
        self.grid = grid
        self.node_id = node_id

    def print_me(self):
        print 'Node %i' % self.node_id
        print '\t le: %e %e %e' % (self.left_edge[0], self.left_edge[1], 
                                   self.left_edge[2])
        print '\t re: %e %e %e' % (self.right_edge[0], self.right_edge[1], 
                                   self.right_edge[2])
        print '\t grid: %i' % self.grid


def get_left_edge(Node node):
    le = np.empty(3, dtype='float64')
    for i in range(3):
        le[i] = node.left_edge[i]
    return le

def get_right_edge(Node node):
    re = np.empty(3, dtype='float64')
    for i in range(3):
        re[i] = node.right_edge[i]
    return re

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _lchild_id(int node_id):
    return (node_id<<1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _rchild_id(int node_id):
    return (node_id<<1) + 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _parent_id(int node_id):
    return (node_id-1) >> 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int should_i_build(Node node, int rank, int size):
    return 1
    # if (node.node_id < size) or (node.node_id >= 2*size):
    #     return 1 
    # elif node.node_id - size == rank:
    #     return 1 
    # else:
    #     return 0 

def kd_traverse(Node trunk, viewpoint=None):
    if viewpoint is None:
        for node in depth_traverse(trunk):
            if kd_is_leaf(node) and node.grid != -1:
                yield node
    else:
        for node in viewpoint_traverse(trunk, viewpoint):
            if kd_is_leaf(node) and node.grid != -1:
                yield node


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def add_grid(Node node, 
                   np.ndarray[np.float64_t, ndim=1] gle, 
                   np.ndarray[np.float64_t, ndim=1] gre, 
                   int gid, 
                   int rank,
                   int size):

    if not should_i_build(node, rank, size):
        return

    if kd_is_leaf(node):
        insert_grid(node, gle, gre, gid, rank, size)
    else:
        less_id = gle[node.split.dim] < node.split.pos
        if less_id:
            add_grid(node.left, gle, gre,
                     gid, rank, size)

        greater_id = gre[node.split.dim] > node.split.pos
        if greater_id:
            add_grid(node.right, gle, gre,
                     gid, rank, size)
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def insert_grid(Node node, 
                      np.ndarray[np.float64_t, ndim=1] gle, 
                      np.ndarray[np.float64_t, ndim=1] gre, 
                      int grid_id, 
                      int rank,
                      int size):
    if not should_i_build(node, rank, size):
        return

    # If we should continue to split based on parallelism, do so!
    # if should_i_split(node, rank, size):
    #     geo_split(node, gle, gre, grid_id, rank, size)
    #     return
    cdef int contained = 1
    for i in range(3):
        if gle[i] > node.left_edge[i] or\
           gre[i] < node.right_edge[i]:
            contained *= 0

    if contained == 1:
        node.grid = grid_id 
        assert(node.grid != -1)
        return

    # Split the grid
    cdef int check = split_grid(node, gle, gre, grid_id, rank, size)
    # If check is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if check == -1:
        node.grid = -1 
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def add_grids(Node node, 
                    int ngrids,
                    np.ndarray[np.float64_t, ndim=2] gles, 
                    np.ndarray[np.float64_t, ndim=2] gres, 
                    np.ndarray[np.int64_t, ndim=1] gids, 
                    int rank,
                    int size):
    cdef int i, nless, ngreater
    if not should_i_build(node, rank, size):
        return

    if kd_is_leaf(node):
        insert_grids(node, ngrids, gles, gres, gids, rank, size)
        return

    less_ids = gles[:,node.split.dim] < node.split.pos
    greater_ids = gres[:,node.split.dim] > node.split.pos
    nless = 0
    ngreater = 0
    for i in xrange(ngrids):
        nless += less_ids[i]
        ngreater += greater_ids[i]
        
    if nless > 0:
        add_grids(node.left, nless, gles[less_ids], gres[less_ids],
                  gids[less_ids], rank, size)

    if ngreater > 0:
        add_grids(node.right, ngreater, gles[greater_ids], gres[greater_ids],
                  gids[greater_ids], rank, size)
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int should_i_split(Node node, int rank, int size):
    if node.node_id < size:
        return 1
    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void insert_grids(Node node, 
                       int ngrids,
                       np.ndarray[np.float64_t, ndim=2] gles, 
                       np.ndarray[np.float64_t, ndim=2] gres, 
                       np.ndarray[np.int64_t, ndim=1] gids, 
                       int rank,
                       int size):
    
    if not should_i_build(node, rank, size) or ngrids == 0:
        return
    cdef int contained = 1
    cdef int check

    if ngrids == 1:
        # If we should continue to split based on parallelism, do so!
        #if should_i_split(node, rank, size):
        #    geo_split(node, gles, gres, grid_ids, rank, size)
        #    return

        for i in range(3):
            contained *= gles[0, i] <= node.left_edge[i]
            contained *= gres[0, i] >= node.right_edge[i]
    
        if contained == 1:
            # print 'Node fully contained, setting to grid: %i' % gids[0]
            node.grid = gids[0]
            assert(node.grid != -1)
            return

    # Split the grids
    check = split_grids(node, ngrids, gles, gres, gids, rank, size)
    # If check is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if check == -1:
        node.grid = -1
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def split_grid(Node node, 
               np.ndarray[np.float64_t, ndim=1] gle, 
               np.ndarray[np.float64_t, ndim=1] gre, 
               int gid, 
               int rank,
               int size):
    # Find a Split
    data = np.empty((1, 2, 3), dtype='float64')
    for i in range(3):
        data[0, 0, i] = gle[i]
        data[0, 1, i] = gre[i]
        # print 'Single Data: ', gle[i], gre[i]

    le = np.empty(3)
    re = np.empty(3)
    for i in range(3):
        le[i] = node.left_edge[i]
        re[i] = node.right_edge[i]

    best_dim, split_pos, less_id, greater_id = \
        kdtree_get_choices(1, data, le, re)

    # If best_dim is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if best_dim == -1:
        print 'Failed to split grid.'
        return -1

        
    split = <Split *> malloc(sizeof(Split))
    split.dim = best_dim
    split.pos = split_pos

    #del data

    # Create a Split
    divide(node, split)

    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    if less_id == 1:
        insert_grid(node.left, gle, gre,
                     gid, rank, size)

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    if greater_id == 1:
        insert_grid(node.right, gle, gre,
                     gid, rank, size)

    return 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef kdtree_get_choices(int n_grids, 
                        np.ndarray[np.float64_t, ndim=3] data,
                        np.ndarray[np.float64_t, ndim=1] l_corner,
                        np.ndarray[np.float64_t, ndim=1] r_corner):
    cdef int i, j, k, dim, n_unique, best_dim, n_best, addit, my_split
    cdef np.float64_t **uniquedims, *uniques, split
    uniquedims = <np.float64_t **> alloca(3 * sizeof(np.float64_t*))
    for i in range(3):
        uniquedims[i] = <np.float64_t *> \
                alloca(2*n_grids * sizeof(np.float64_t))
    my_max = 0
    my_split = 0
    best_dim = -1
    for dim in range(3):
        n_unique = 0
        uniques = uniquedims[dim]
        for i in range(n_grids):
            # Check for disqualification
            for j in range(2):
                # print "Checking against", i,j,dim,data[i,j,dim]
                if not (l_corner[dim] < data[i, j, dim] and
                        data[i, j, dim] < r_corner[dim]):
                    # print "Skipping ", data[i,j,dim], l_corner[dim], r_corner[dim]
                    continue
                skipit = 0
                # Add our left ...
                for k in range(n_unique):
                    if uniques[k] == data[i, j, dim]:
                        skipit = 1
                        # print "Identified", uniques[k], data[i,j,dim], n_unique
                        break
                if skipit == 0:
                    uniques[n_unique] = data[i, j, dim]
                    n_unique += 1
        if n_unique > my_max:
            best_dim = dim
            my_max = n_unique
            my_split = (n_unique-1)/2
    # I recognize how lame this is.
    cdef np.ndarray[np.float64_t, ndim=1] tarr = np.empty(my_max, dtype='float64')
    for i in range(my_max):
        # print "Setting tarr: ", i, uniquedims[best_dim][i]
        tarr[i] = uniquedims[best_dim][i]
    tarr.sort()
    split = tarr[my_split]
    cdef np.ndarray[np.uint8_t, ndim=1] less_ids = np.empty(n_grids, dtype='uint8')
    cdef np.ndarray[np.uint8_t, ndim=1] greater_ids = np.empty(n_grids, dtype='uint8')
    for i in range(n_grids):
        if data[i, 0, best_dim] < split:
            less_ids[i] = 1
        else:
            less_ids[i] = 0
        if data[i, 1, best_dim] > split:
            greater_ids[i] = 1
        else:
            greater_ids[i] = 0
    # Return out unique values
    return best_dim, split, less_ids.view('bool'), greater_ids.view('bool')

#@cython.boundscheck(True)
#@cython.wraparound(False)
#@cython.cdivision(True)
cdef int split_grids(Node node, 
                       int ngrids,
                       np.ndarray[np.float64_t, ndim=2] gles, 
                       np.ndarray[np.float64_t, ndim=2] gres, 
                       np.ndarray[np.int64_t, ndim=1] gids, 
                       int rank,
                       int size):
    # Find a Split
    cdef int i, j, k

    le = get_left_edge(node)
    re = get_right_edge(node)

    data = np.array([(gles[i,:], gres[i,:]) for i in
        xrange(ngrids)], copy=False)
    best_dim, split_pos, less_ids, greater_ids = \
        kdtree_get_choices(ngrids, data, le, re)
 
    # If best_dim is -1, then we have found a place where there are no choices.
    # Exit out and set the node to None.
    if best_dim == -1:
        print 'Failed to split grids.'
        return -1

    split = <Split *> malloc(sizeof(Split))
    split.dim = best_dim
    split.pos = split_pos

    #del data

    # Create a Split
    divide(node, split)

    nless = np.sum(less_ids)
    ngreat = np.sum(greater_ids)
    # Populate Left Node
    #print 'Inserting left node', node.left_edge, node.right_edge
    insert_grids(node.left, nless, gles[less_ids], gres[less_ids],
                 gids[less_ids], rank, size)

    # Populate Right Node
    #print 'Inserting right node', node.left_edge, node.right_edge
    insert_grids(node.right, ngreat, gles[greater_ids], gres[greater_ids],
                 gids[greater_ids], rank, size)

    return 0

# def geo_split_grid(node, gle, gre, grid_id, rank, size):
#     big_dim = np.argmax(gre-gle)
#     new_pos = (gre[big_dim] + gle[big_dim])/2.
#     old_gre = gre.copy()
#     new_gle = gle.copy()
#     new_gle[big_dim] = new_pos
#     gre[big_dim] = new_pos
# 
#     split = Split(big_dim, new_pos)
# 
#     # Create a Split
#     divide(node, split)
# 
#     # Populate Left Node
#     #print 'Inserting left node', node.left_edge, node.right_edge
#     insert_grid(node.left, gle, gre,
#                 grid_id, rank, size)
# 
#     # Populate Right Node
#     #print 'Inserting right node', node.left_edge, node.right_edge
#     insert_grid(node.right, new_gle, old_gre,
#                 grid_id, rank, size)
#     return
# 
# 
# def geo_split(node, gles, gres, grid_ids, rank, size):
#     big_dim = np.argmax(gres[0]-gles[0])
#     new_pos = (gres[0][big_dim] + gles[0][big_dim])/2.
#     old_gre = gres[0].copy()
#     new_gle = gles[0].copy()
#     new_gle[big_dim] = new_pos
#     gres[0][big_dim] = new_pos
#     gles = np.append(gles, np.array([new_gle]), axis=0)
#     gres = np.append(gres, np.array([old_gre]), axis=0)
#     grid_ids = np.append(grid_ids, grid_ids, axis=0)
# 
#     split = Split(big_dim, new_pos)
# 
#     # Create a Split
#     divide(node, split)
# 
#     # Populate Left Node
#     #print 'Inserting left node', node.left_edge, node.right_edge
#     insert_grids(node.left, gles[:1], gres[:1],
#             grid_ids[:1], rank, size)
# 
#     # Populate Right Node
#     #print 'Inserting right node', node.left_edge, node.right_edge
#     insert_grids(node.right, gles[1:], gres[1:],
#             grid_ids[1:], rank, size)
#     return

cdef new_right(Node node, Split * split):
    new_right = Node.right_edge.copy()
    new_right[split.dim] = split.pos
    return new_right

cdef new_left(Node node, Split * split):
    new_left = Node.left_edge.copy()
    new_left[split.dim] = split.pos
    return new_left

cdef void divide(Node node, Split * split):
    # Create a Split
    node.split = split
    
    lle = np.zeros(3, dtype='float64')
    lre = np.zeros(3, dtype='float64')
    rle = np.zeros(3, dtype='float64')
    rre = np.zeros(3, dtype='float64')

    cdef int i
    for i in range(3):
        lle[i] = node.left_edge[i]
        lre[i] = node.right_edge[i]
        rle[i] = node.left_edge[i]
        rre[i] = node.right_edge[i]
    lre[split.dim] = split.pos
    rle[split.dim] = split.pos

    node.left = Node(node, None, None,
            lle.copy(), lre.copy(), node.grid,
                     _lchild_id(node.node_id))

    node.right = Node(node, None, None,
            rle.copy(), rre.copy(), node.grid,
                      _rchild_id(node.node_id))
    
    return
# 
def kd_sum_volume(Node node):
    cdef np.float64_t vol = 1.0
    if (node.left is None) and (node.right is None):
        if node.grid == -1:
            return 0.0
        for i in range(3):
            vol *= node.right_edge[i] - node.left_edge[i]
        return vol 
    else:
        return kd_sum_volume(node.left) + kd_sum_volume(node.right)
# 
# def kd_sum_cells(node):
#     if (node.left is None) and (node.right is None):
#         if node.grid is None:
#             return 0.0
#         return np.prod(node.right_edge - node.left_edge)
#     else:
#         return kd_sum_volume(node.left) + kd_sum_volume(node.right)
# 
# 

def kd_node_check(Node node):
    assert (node.left is None) == (node.right is None)
    if (node.left is None) and (node.right is None):
        if node.grid != -1:
            return np.prod(node.right_edge - node.left_edge)
        else: return 0.0
    else:
        return kd_node_check(node.left)+kd_node_check(node.right)


def kd_is_leaf(Node node):
    cdef int has_l_child = node.left == None
    cdef int has_r_child = node.right == None
    assert has_l_child == has_r_child
    return has_l_child

def step_depth(Node current, Node previous):
    '''
    Takes a single step in the depth-first traversal
    '''
    if kd_is_leaf(current): # At a leaf, move back up
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down, go left first
        previous = current
        if current.left is not None:
            current = current.left
        elif current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left, go right 
        previous = current
        if current.right is not None:
            current = current.right
        else:
            current = current.parent

    elif current.right is previous: # Moving up from right child, move up
        previous = current
        current = current.parent

    return current, previous
 
def depth_traverse(Node trunk, max_node=None):
    '''
    Yields a depth-first traversal of the kd tree always going to
    the left child before the right.
    '''
    current = trunk
    previous = None
    if max_node is None:
        max_node = np.inf
    while current is not None:
        yield current
        current, previous = step_depth(current, previous)
        if current is None: break
        if current.node_id >= max_node:
            current = current.parent
            previous = current.right
# 
# def depth_first_touch(tree, max_node=None):
#     '''
#     Yields a depth-first traversal of the kd tree always going to
#     the left child before the right.
#     '''
#     current = tree.trunk
#     previous = None
#     if max_node is None:
#         max_node = np.inf
#     while current is not None:
#         if previous is None or previous.parent != current:
#             yield current
#         current, previous = step_depth(current, previous)
#         if current is None: break
#         if current.id >= max_node:
#             current = current.parent
#             previous = current.right
# 
# def breadth_traverse(tree):
#     '''
#     Yields a breadth-first traversal of the kd tree always going to
#     the left child before the right.
#     '''
#     current = tree.trunk
#     previous = None
#     while current is not None:
#         yield current
#         current, previous = step_depth(current, previous)
# 
# 
def viewpoint_traverse(tree, viewpoint):
    '''
    Yields a viewpoint dependent traversal of the kd-tree.  Starts
    with nodes furthest away from viewpoint.
    '''

    current = tree.trunk
    previous = None
    while current is not None:
        yield current
        current, previous = step_viewpoint(current, previous, viewpoint)

def step_viewpoint(Node current, 
                   Node previous, 
                   viewpoint):
    '''
    Takes a single step in the viewpoint based traversal.  Always
    goes to the node furthest away from viewpoint first.
    '''
    if kd_is_leaf(current): # At a leaf, move back up
        previous = current
        current = current.parent
    elif current.split.dim is None: # This is a dead node
        previous = current
        current = current.parent

    elif current.parent is previous: # Moving down
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                previous = current.right
        else:
            if current.left is not None:
                current = current.left
            else:
                previous = current.left

    elif current.right is previous: # Moving up from right 
        previous = current
        if viewpoint[current.split.dim] <= current.split.pos:
            if current.left is not None:
                current = current.left
            else:
                current = current.parent
        else:
            current = current.parent

    elif current.left is previous: # Moving up from left child
        previous = current
        if viewpoint[current.split.dim] > current.split.pos:
            if current.right is not None:
                current = current.right
            else:
                current = current.parent
        else:
            current = current.parent

    return current, previous
# 
# 
# def receive_and_reduce(comm, incoming_rank, image, add_to_front):
#     mylog.debug( 'Receiving image from %04i' % incoming_rank)
#     #mylog.debug( '%04i receiving image from %04i'%(self.comm.rank,back.owner))
#     arr2 = comm.recv_array(incoming_rank, incoming_rank).reshape(
#         (image.shape[0], image.shape[1], image.shape[2]))
# 
#     if add_to_front:
#         front = arr2
#         back = image
#     else:
#         front = image
#         back = arr2
# 
#     if image.shape[2] == 3:
#         # Assume Projection Camera, Add
#         np.add(image, front, image)
#         return image
# 
#     ta = 1.0 - front[:,:,3]
#     np.maximum(ta, 0.0, ta)
#     # This now does the following calculation, but in a memory
#     # conservative fashion
#     # image[:,:,i  ] = front[:,:,i] + ta*back[:,:,i]
#     image = back.copy()
#     for i in range(4):
#         np.multiply(image[:,:,i], ta, image[:,:,i])
#     np.add(image, front, image)
#     return image
# 
# def send_to_parent(comm, outgoing_rank, image):
#     mylog.debug( 'Sending image to %04i' % outgoing_rank)
#     comm.send_array(image, outgoing_rank, tag=comm.rank)
# 
# def scatter_image(comm, root, image):
#     mylog.debug( 'Scattering from %04i' % root)
#     image = comm.mpi_bcast(image, root=root)
#     return image
# 
# def find_node(node, pos):
#     """
#     Find the AMRKDTree node enclosing a position
#     """
#     assert(np.all(node.left_edge <= pos))
#     assert(np.all(node.right_edge > pos))
#     while not kd_is_leaf(node):
#         if pos[node.split.dim] < node.split.pos:
#             node = node.left
#         else:
#             node = node.right
#     return node


