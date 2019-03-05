import cython
import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp cimport bool as cbool
from cpython cimport bool as pybool
from cython.operator cimport dereference

from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

cdef class PyNode:
    r"""A container for leaf info.

    Attributes:
        npts (np.uint64_t): Number of points in this node.
        ndim (np.uint32_t): Number of dimensions in domain.
        num_leaves (np.uint32_t): Number of leaves in the tree containing this
            node.
        start_idx (np.uint64_t): Index where indices for this node begin.
        stop_idx (np.uint64_t): One passed the end of indices for this node.
        left_edge (np.ndarray of float64): Minimum bounds of this node in each
            dimension.
        right_edge (np.ndarray of float64): Maximum bounds of this node in each
            dimension.
        periodic_left (np.ndarray of bool): Periodicity of minimum bounds.
        periodic_right (np.ndarray of bool): Periodicity of maximum bounds.
        domain_width (np.ndarray of float64): Width of the total domain in each
            dimension.
        left_neighbors (list of lists): Indices of neighbor leaves at the
            minimum bounds in each dimension.
        right_neighbors (list of lists): Indices of neighbor leaves at the
            maximum bounds in each dimension.

    """

    cdef void _init_node(self, Node* node, uint32_t num_leaves,
                         double *domain_width):
        cdef np.uint32_t i, j
        self._node = node
        self.id = node.leafid
        self.npts = node.children
        self.ndim = node.ndim
        self.num_leaves = num_leaves
        self.start_idx = node.left_idx
        self.stop_idx = (node.left_idx + node.children)
        self._domain_width = domain_width
        self.left_neighbors = [None for i in range(self.ndim)]
        self.right_neighbors = [None for i in range(self.ndim)]
        for i in range(self.ndim):
            self.left_neighbors[i] = [node.left_neighbors[i][j] for j in
                                      range(node.left_neighbors[i].size())]
            self.right_neighbors[i] = [node.right_neighbors[i][j] for j in
                                       range(node.right_neighbors[i].size())]

    def __cinit__(self):
        # Initialize everthing to NULL/0/None to prevent seg fault
        self._node = NULL
        self.id = 0
        self.npts = 0
        self.ndim = 0
        self.num_leaves = 0
        self.start_idx = 0
        self.stop_idx = 0
        self._domain_width = NULL
        self.left_neighbors = None
        self.right_neighbors = None

    def __init__(self):
        pass

    def __repr__(self):
        nchars = 1 + len(str(self.__class__.__name__))
        return ('%s(id=%i, npts=%i, start_idx=%i, stop_idx=%i,\n' +
                ' ' * nchars + 'left_edge=%s,\n' +
                ' ' * nchars + 'right_edge=%s)') % (
            self.__class__.__name__,
            self.id,
            self.npts,
            self.start_idx,
            self.stop_idx,
            self.left_edge,
            self.right_edge,
        )

    @property
    def periodic_left(self):
        cdef cbool[:] view = <cbool[:self.ndim]> self._node.periodic_left
        return np.asarray(view)
    @property
    def periodic_right(self):
        cdef cbool[:] view = <cbool[:self.ndim]> self._node.periodic_right
        return np.asarray(view)
    @property
    def left_edge(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._node.left_edge
        return np.asarray(view)
    @property
    def right_edge(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._node.right_edge
        return np.asarray(view)
    @property
    def domain_width(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._domain_width
        return np.asarray(view)

    @property
    def slice(self):
        """slice: Slice of kdtree indices contained by this node."""
        return slice(self.start_idx, self.stop_idx)

    @property
    def neighbors(self):
        """list of int: Indices of all neighboring leaves including this
        leaf."""
        cdef np.uint32_t i
        cdef object out
        cdef vector[uint32_t] vout = self._node.all_neighbors
        out = [vout[i] for i in range(<np.uint32_t>vout.size())]
        return out

    def assert_equal(self, PyNode solf):
        """Assert that node properties are equal."""
        np.testing.assert_equal(self.npts, solf.npts)
        np.testing.assert_equal(self.ndim, solf.ndim)
        np.testing.assert_equal(self.num_leaves, solf.num_leaves)
        np.testing.assert_equal(self.id, solf.id)
        np.testing.assert_equal(self.start_idx, solf.start_idx)
        np.testing.assert_equal(self.stop_idx, solf.stop_idx)
        np.testing.assert_array_equal(self.left_edge, solf.left_edge)
        np.testing.assert_array_equal(self.right_edge, solf.right_edge)
        np.testing.assert_array_equal(self.periodic_left, solf.periodic_left)
        np.testing.assert_array_equal(self.periodic_right, solf.periodic_right)
        for i in range(self.ndim):
            np.testing.assert_equal(self.left_neighbors[i], solf.left_neighbors[i])
            np.testing.assert_equal(self.right_neighbors[i], solf.right_neighbors[i])
        np.testing.assert_equal(self.neighbors, solf.neighbors)


cdef class PyKDTree:
    r"""Construct a KDTree for a set of points.

    Args:
        pts (np.ndarray of double): (n,m) array of n coordinates in a
            m-dimensional domain.
        left_edge (np.ndarray of double): (m,) domain minimum in each dimension.
        right_edge (np.ndarray of double): (m,) domain maximum in each dimension.
        periodic (bool or np.ndarray of bool, optional): Truth of the domain
            periodicity overall (if bool), or in each dimension (if np.ndarray).
            Defaults to `False`.
        leafsize (int, optional): The maximum number of points that should be in
            a leaf. Defaults to 10000.
        nleaves (int, optional): The number of leaves that should be in the
            resulting tree. If greater than 0, leafsize is adjusted to produce a
            tree with 2**(ceil(log2(nleaves))) leaves. The leafsize keyword
            argument is ignored if nleaves is greater zero. Defaults to 0.
        data_version (int, optional): An optional user-provided integer that
            can be used to uniquely identify the data used to generate the
            KDTree. This is useful if you save the kdtree to disk and restore
            it later and need to verify that the underlying data is the same.
        use_sliding_midpoint (bool, optional): If True, the sliding midpoint
            rule is used to perform splits. Otherwise, the median is used.
            Defaults to False.

    Raises:
        ValueError: If `leafsize < 2`. This currectly segfaults.

    Attributes:
        npts (uint64): Number of points in the tree.
        ndim (uint32): Number of dimensions points occupy.
        data_version (int64): User-provided version number (defaults to 0)
        num_leaves (uint32): Number of leaves in the tree.
        leafsize (uint32): Maximum number of points a leaf can have.
        leaves (list of `cykdtree.PyNode`): Tree leaves.
        idx (np.ndarray of uint64): Indices sorting the points by leaf.
        left_edge (np.ndarray of double): (m,) domain minimum in each dimension.
        right_edge (np.ndarray of double): (m,) domain maximum in each dimension.
        domain_width (np.ndarray of double): (m,) domain width in each dimension.
        periodic (np.ndarray of bool): Truth of domain periodicity in each
            dimension.

    """

    cdef void _init_tree(self, KDTree* tree):
        self._tree = tree
        self.ndim = tree.ndim
        self.data_version = tree.data_version
        self.npts = tree.npts
        self.num_leaves = tree.num_leaves
        self.leafsize = tree.leafsize
        self._make_leaves()
        self._idx = np.empty(self.npts, 'uint64')
        cdef uint64_t i
        for i in range(self.npts):
            self._idx[i] = tree.all_idx[i]

    def __cinit__(self):
        # Initialize everthing to NULL/0/None to prevent seg fault
        self._tree = NULL
        self.npts = 0
        self.ndim = 0
        self.num_leaves = 0
        self.leafsize = 0
        self._left_edge = NULL
        self._right_edge = NULL
        self._periodic = NULL
        self.leaves = None
        self._idx = None

    def __init__(self, np.ndarray[double, ndim=2] pts = None,
                 left_edge = None,
                 right_edge = None,
                 periodic = False,
                 int leafsize = 10000,
                 int nleaves = 0,
                 data_version = None,
                 use_sliding_midpoint = False):
        # Return with nothing set if points not provided
        if pts is None:
            return
        # Set leafsize of number of leaves provided
        if nleaves > 0:
            nleaves = <int>(2**np.ceil(np.log2(<float>nleaves)))
            leafsize = pts.shape[0]/nleaves + 1
        if (leafsize < 2):
            # This is here to prevent segfault. The cpp code needs modified to
            # support leafsize = 1
            raise ValueError("'leafsize' cannot be smaller than 2.")
        if left_edge is None:
            left_edge = np.min(pts, axis=0)
        else:
            left_edge = np.array(left_edge)
        if right_edge is None:
            right_edge = np.max(pts, axis=0)
        else:
            right_edge = np.array(right_edge)
        if data_version is None:
            data_version = 0
        self.data_version = data_version
        cdef uint32_t k,i,j
        self.npts = <uint64_t>pts.shape[0]
        self.ndim = <uint32_t>pts.shape[1]
        assert(left_edge.size == self.ndim)
        assert(right_edge.size == self.ndim)
        self.leafsize = leafsize
        self._left_edge = <double *>malloc(self.ndim*sizeof(double))
        self._right_edge = <double *>malloc(self.ndim*sizeof(double))
        self._periodic = <cbool *>malloc(self.ndim*sizeof(cbool));
        for i in range(self.ndim):
            self._left_edge[i] = left_edge[i]
            self._right_edge[i] = right_edge[i]
        if isinstance(periodic, pybool):
            for i in range(self.ndim):
                self._periodic[i] = <cbool>periodic
        else:
            for i in range(self.ndim):
                self._periodic[i] = <cbool>periodic[i]
        # Create tree and leaves
        self._make_tree(&pts[0,0], <cbool>use_sliding_midpoint)
        self._make_leaves()

    def __dealloc__(self):
        if self._tree != NULL:
            del self._tree
        if self._left_edge != NULL:
            free(self._left_edge)
        if self._right_edge != NULL:
            free(self._right_edge)
        if self._periodic != NULL:
            free(self._periodic)

    def assert_equal(self, PyKDTree solf, pybool strict_idx = True):
        r"""Compare this tree to another tree.

        Args:
            solf (PyKDTree): Another KDTree to compare with this one.
            strict_idx (bool, optional): If True, the index vectors are
                compared for equality element by element. If False,
                corresponding leaves must contain the same indices, but they
                can be in any order. Defaults to True.

        Raises:
            AssertionError: If there are missmatches between any of the two
                trees' parameters.

        """
        np.testing.assert_equal(self.npts, solf.npts)
        np.testing.assert_equal(self.ndim, solf.ndim)
        np.testing.assert_equal(self.data_version, solf.data_version)
        np.testing.assert_equal(self.leafsize, solf.leafsize)
        np.testing.assert_equal(self.num_leaves, solf.num_leaves)
        np.testing.assert_array_equal(self.left_edge, solf.left_edge)
        np.testing.assert_array_equal(self.right_edge, solf.right_edge)
        np.testing.assert_array_equal(self.periodic, solf.periodic)
        # Compare index at the leaf level since we only care that the leaves
        # contain the same points
        if strict_idx:
            np.testing.assert_array_equal(self._idx, solf._idx)
        for i in range(self.num_leaves):
            self.leaves[i].assert_equal(solf.leaves[i])
            if not strict_idx:
                np.testing.assert_array_equal(
                    np.sort(self._idx[self.leaves[i].slice]),
                    np.sort(solf._idx[solf.leaves[i].slice]))

    cdef void _make_tree(self, double *pts, bool use_sliding_midpoint):
        r"""Carry out creation of KDTree at C++ level."""
        cdef uint64_t[:] idx = np.arange(self.npts).astype('uint64')
        self._tree = new KDTree(pts, &idx[0], self.npts, self.ndim, self.leafsize,
                                self._left_edge, self._right_edge, self._periodic,
                                self.data_version, use_sliding_midpoint)
        self._idx = idx

    cdef void _make_leaves(self):
        r"""Create a list of Python leaf objects from C++ leaves."""
        self.num_leaves = <uint32_t>self._tree.leaves.size()
        self.leaves = [None for _ in xrange(self.num_leaves)]
        cdef Node* leafnode
        cdef PyNode leafnode_py
        cdef object leaf_neighbors = None
        for k in xrange(self.num_leaves):
            leafnode = self._tree.leaves[k]
            leafnode_py = PyNode()
            leafnode_py._init_node(leafnode, self.num_leaves,
                                   self._tree.domain_width)
            self.leaves[leafnode.leafid] = leafnode_py

    @property
    def left_edge(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._tree.domain_left_edge
        return np.asarray(view)
    @property
    def right_edge(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._tree.domain_right_edge
        return np.asarray(view)
    @property
    def domain_width(self):
        cdef np.float64_t[:] view = <np.float64_t[:self.ndim]> self._tree.domain_width
        return np.asarray(view)
    @property
    def periodic(self):
        cdef cbool[:] view = <cbool[:self.ndim]> self._tree.periodic
        # return np.asarray(view)
        cdef object out = np.empty(self.ndim, 'bool')
        cdef np.uint32_t i
        for i in range(self.ndim):
            out[i] = view[i]
        return out

    def leaf_idx(self, np.uint32_t leafid):
        r"""Get array of indices for points in a leaf.

        Args:
            leafid (np.uint32_t): Unique index of the leaf in question.

        Returns:
            np.ndarray of np.uint64_t: Indices of points belonging to leaf.

        """
        cdef np.ndarray[np.uint64_t] out = self._idx[self.leaves[leafid].slice]
        return out

    cdef np.ndarray[np.uint32_t, ndim=1] _get_neighbor_ids(self, np.ndarray[double, ndim=1] pos):
        cdef np.uint32_t i
        cdef vector[uint32_t] vout = self._tree.get_neighbor_ids(&pos[0]);
        cdef np.ndarray[np.uint32_t, ndim=1] out = np.empty(vout.size(), 'uint32')
        for i in xrange(vout.size()):
            out[i] = vout[i]
        return out

    @property
    def idx(self):
        return np.asarray(self._idx)

    def get_neighbor_ids(self, np.ndarray[double, ndim=1] pos):
        r"""Return the IDs of leaves containing & neighboring a given position.

        Args:
            pos (np.ndarray of double): Coordinates.

        Returns:
            np.ndarray of uint32: Leaves containing/neighboring `pos`.

        Raises:
            ValueError: If pos is not contained withing the KDTree.

        """
        return self._get_neighbor_ids(pos)

    cdef np.ndarray[np.uint32_t, ndim=1] _get_neighbor_ids_3(self, np.float64_t pos[3]):
        cdef np.uint32_t i
        cdef vector[uint32_t] vout = self._tree.get_neighbor_ids(&pos[0]);
        cdef np.ndarray[np.uint32_t, ndim=1] out = np.empty(vout.size(), 'uint32')
        for i in xrange(vout.size()):
            out[i] = vout[i]
        return out

    cdef PyNode _get(self, np.ndarray[double, ndim=1] pos):
        assert(<uint32_t>len(pos) == self.ndim)
        cdef Node* leafnode = self._tree.search(&pos[0])
        if leafnode == NULL:
            raise ValueError("Position is not within the kdtree root node.")
        cdef PyNode out = self.leaves[leafnode.leafid]
        return out

    def get(self, np.ndarray[double, ndim=1] pos):
        r"""Return the leaf containing a given position.

        Args:
            pos (np.ndarray of double): Coordinates.

        Returns:
            :class:`cykdtree.PyNode`: Leaf containing `pos`.

        Raises:
            ValueError: If pos is not contained withing the KDTree.

        """
        return self._get(pos)

    def consolidate_edges(self):
        r"""Return arrays of the left and right edges for all leaves in the
        tree on each process.

        Returns:
            tuple(np.ndarray of double, np.ndarray of double): The left (first
                array) and right (second array) edges of each leaf (1st array
                dimension), in each dimension (2nd array dimension).

        """
        cdef np.ndarray[np.float64_t, ndim=2] leaves_le
        cdef np.ndarray[np.float64_t, ndim=2] leaves_re
        leaves_le = np.empty((self.num_leaves, self.ndim), 'float64')
        leaves_re = np.empty((self.num_leaves, self.ndim), 'float64')
        self._tree.consolidate_edges(&leaves_le[0,0], &leaves_re[0,0])
        return (leaves_le, leaves_re)

    def save(self, str filename):
        r"""Saves the PyKDTree to disk as raw binary data.

        Note that this file may not necessarily be portable.

        Args:
            filename (string): Name of the file to serialize the kdtree to

        """
        cdef KDTree* my_tree = self._tree
        cdef ofstream* outputter = new ofstream(filename.encode('utf8'), binary)
        try:
            my_tree.serialize(dereference(outputter))
        finally:
            del outputter

    @classmethod
    def from_file(cls, str filename, data_version=None):
        r"""Create a PyKDTree from a binary file created by ``PyKDTree.save()``

        Note that loading a file created on another machine may create
        a corrupted PyKDTree instance.

        Args:
            filename (string): Name of the file to load the kdtree from
            data_version (int): A unique integer corresponding to the data 
                                being loaded. If the loaded data_version does 
                                not match the data_version supplied here then 
                                an OSError is raised. Optional.
        
        Returns:
            :class:`cykdtree.PyKDTree`: A KDTree restored from the file

        """
        cdef ifstream* inputter = new ifstream(filename.encode(), binary)
        cdef PyKDTree ret = cls()
        if data_version is None:
            data_version = 0
        try:
            ret._init_tree(new KDTree(dereference(inputter)))
        finally:
            del inputter
        return ret
