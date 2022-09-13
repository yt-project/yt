# distutils: libraries = STD_LIBS
"""
CyOctree building, loading and refining routines
"""
cimport cython
cimport libc.math as math
cimport numpy as np

import numpy as np

from libc.stdlib cimport free, malloc

from yt.geometry.particle_deposit cimport get_kernel_func, kernel_func

np.import_array()


cdef struct Octree:
    # Array of 3*num_nodes [x1, y1, z1, x2, y2, z2, ...]
    np.float64_t * node_positions
    # 1 or 0 of whether the oct has refined to make children
    np.uint8_t * refined
    # Each oct stores the depth in the tree: the root node is 0
    np.uint8_t * depth
    # pstart and pend tell us which particles in the pidx are stored in each oct
    np.int64_t * pstart
    np.int64_t * pend
    # This tells us the index of each child
    np.int64_t * children

    # Here we store the coordinates and IDs of all the particles in the tree
    np.float64_t * pposx
    np.float64_t * pposy
    np.float64_t * pposz
    np.int64_t * pidx

    # The max number of particles per leaf, if above, we refine
    np.int64_t n_ref

    # The number of particles in our tree, e.g the length of ppos and pidx
    np.int64_t num_particles

    # The total size of the octree (x, y, z)
    np.float64_t * size

    # The current number of nodes in the octree
    np.int64_t num_nodes

    # The maximum depth before we stop refining, and the maximum number of nodes
    # we can fit in our array
    np.uint8_t max_depth
    np.int64_t max_num_nodes


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_build_node(Octree * tree, long int node_idx):
    """
    This is the main function in the building of the octree. This function takes
    in the tree and an index of a oct to process. If the node has too many
    particles (> n_ref) and the depth is less than the maximum tree depth then
    we refine.

    The `refine` creates 8 sub octs, within our oct. We indenfity the particles
    in each of the subocts. We then recursively call this function on each of
    the subocts.

    Parameters
    ----------
    tree : Octree *
        A pointer to the octree
    node_idx : long int
        The index of the current node we are processing

    Returns
    -------
    int
        Success of tree build
    """
    cdef np.int64_t splits[9]
    cdef np.int64_t i, j, k, n, start, end
    cdef np.float64_t lx, ly, lz, sz, inv_size

    # If we are running out of space in our tree, then we *try* to
    # relloacate a tree of double the size
    if (tree.num_nodes + 8) >= tree.max_num_nodes:
        if octree_reallocate(tree, tree.max_num_nodes * 2):
            return 1

    if (
        (tree.pend[node_idx] - tree.pstart[node_idx] > tree.n_ref) and
        (tree.depth[node_idx] < tree.max_depth)
    ):
        tree.refined[node_idx] = 1

        # As we have decided to refine, we need to know which of the particles
        # in this oct will go into each of the 8 child octs
        split_helper(tree, node_idx, splits)

        # Figure out the size of the current oct
        inv_size = 1. / 2.**tree.depth[node_idx]
        sx = tree.size[0] * inv_size
        sy = tree.size[1] * inv_size
        sz = tree.size[2] * inv_size
        lx = tree.node_positions[(3*node_idx)] - sx/2.
        ly = tree.node_positions[(3*node_idx)+1] - sy/2.
        lz = tree.node_positions[(3*node_idx)+2] - sz/2.

        # Loop through and generate the children AND recursively refine...
        n = 0
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    start = splits[n]
                    end = splits[n + 1]
                    child = tree.num_nodes

                    # Store the child location
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
                    if octree_build_node(tree, child):
                        return 1
                    n += 1

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_allocate(Octree * octree, long int num_nodes):
    """
    This is the main allocation function in the octree. We allocate all of the
    arrays we require to store information about every single oct in the tree

    Parameters
    ----------
    octree : Octree *
        A pointer to the octree
    num_nodes : long int
        The maximum number of nodes to allocate for

    Returns
    -------
    int
        Success of allocation
    """
    octree.node_positions = <np.float64_t *> malloc(
        num_nodes * 3 * sizeof(np.float64_t))
    if octree.node_positions == NULL:
        return 1

    octree.size = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
    if octree.size == NULL:
        return 1

    octree.children = <np.int64_t *> malloc(8 * num_nodes * sizeof(np.int64_t))
    if octree.children == NULL:
        return 1

    octree.pstart = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    if octree.pstart == NULL:
        return 1

    octree.pend = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    if octree.pend == NULL:
        return 1

    octree.refined = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    if octree.refined == NULL:
        return 1

    octree.depth = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    if octree.depth == NULL:
        return 1

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int octree_reallocate(Octree * octree, long int num_nodes):
    """
    This function re-allocates all of the arrays malloc'd in `octree_allocate`
    See Notes for when we want to re-allocate.

    Parameters
    ----------
    octree : Octree *
        A pointer to the octree
    num_nodes : long int
        The maximum number of nodes to (re)allocate for

    Returns
    -------
    int
        Success of the reallocation

    Notes
    -----
    Why do we want to re-allocate?

    Well 2 cases,
    1) The octree is still building and we have ran out of space, so we
    have asked for an increased number of nodes and are reallocating each array

    2) We have finished building the octree and we have used less nodes
    than we originally allocated. We are now shrinking the octree and giving
    the spare memory back.
    """
    cdef np.float64_t * old_arr
    cdef np.int64_t * old_arr_int
    cdef np.uint8_t * old_arr_uint
    cdef np.int64_t i

    old_arr = octree.node_positions
    octree.node_positions = <np.float64_t *> malloc(num_nodes * 3 * sizeof(np.float64_t))
    if octree.node_positions == NULL: return 1
    for i in range(3*octree.num_nodes):
        octree.node_positions[i] = old_arr[i]
    free(old_arr)

    old_arr_int = octree.children
    octree.children = <np.int64_t *> malloc(num_nodes * 8 * sizeof(np.int64_t))
    if octree.children == NULL: return 1
    for i in range(8*octree.num_nodes):
        octree.children[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_int = octree.pstart
    octree.pstart = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    if octree.pstart == NULL: return 1
    for i in range(octree.num_nodes):
        octree.pstart[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_int = octree.pend
    octree.pend = <np.int64_t *> malloc(num_nodes * sizeof(np.int64_t))
    if octree.pend == NULL: return 1
    for i in range(octree.num_nodes):
        octree.pend[i] = old_arr_int[i]
    free(old_arr_int)

    old_arr_uint = octree.refined
    octree.refined = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    if octree.refined == NULL: return 1
    for i in range(octree.num_nodes):
        octree.refined[i] = old_arr_uint[i]
    free(old_arr_uint)

    old_arr_uint = octree.depth
    octree.depth = <np.uint8_t *> malloc(num_nodes * sizeof(np.int8_t))
    if octree.depth == NULL: return 1
    for i in range(octree.num_nodes):
        octree.depth[i] = old_arr_uint[i]
    free(old_arr_uint)

    octree.max_num_nodes = num_nodes

    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void octree_deallocate(Octree * octree):
    """
    Just free-ing every array we allocated to ensure we don't leak.

    Parameter
    ---------
    octree : Octree *
        Pointer to the octree
    """
    free(octree.node_positions)
    free(octree.size)
    free(octree.children)

    free(octree.pstart)
    free(octree.pend)

    free(octree.refined)
    free(octree.depth)

    free(octree.pidx)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class CyOctree:
    """
    This a class to store the underlying octree and particle data that can be
    interacted with from both Cython and Python
    """
    cdef Octree * c_octree
    cdef np.float64_t[::1, :] input_positions

    cdef np.int64_t n_ref
    cdef np.float64_t[:] left_edge
    cdef np.float64_t[:] right_edge
    cdef np.uint8_t max_depth

    cdef kernel_func kernel

    def __init__(
        self,
        np.float64_t[:, :] input_pos,
        left_edge=None,
        right_edge=None,
        np.int64_t n_ref=32,
        np.uint8_t max_depth=200
    ):
        """
        Octree initialiser.

        We copy the inputted particle positions and make the root node. We then
        refine the octree until every leaf has either less particles than n_ref
        or is at the maximum depth.

        Finally, we re-allocate all of the memory required by tree to ensure we
        do not use more memory than required.

        Parameters
        ----------
        input_pos : 2D memory view
            Particles positions in the format (num_particles, 3)
        {left,right}_edge : ndarray
            xyz coordinates of the lower left (upper right) corner of the octree.
        n_ref : int, default: 32
            The maximum number of particles per leaf, if more, the oct
            will refine
        max_depth : int, default: 200
            The maximum depth the octree will refine to. If we set
            this too high then we may hit a stack overflow due to the
            recursive nature of the build
        """
        self.n_ref = n_ref
        self.max_depth = max_depth
        self.input_positions = np.asfortranarray(input_pos, dtype=np.float64)

        if self._allocate_octree():
            raise MemoryError("Unable to allocate memory required for octree build.")

        self._make_root(left_edge, right_edge)

        if octree_build_node(self.c_octree, 0):
            raise MemoryError("Unable to allocate memory required for octree build.")

        if octree_reallocate(self.c_octree, self.c_octree.num_nodes):
            raise MemoryError("Unable to allocate memory required for octree build.")

    def __del__(self):
        """
        Make sure we clean up properly!
        """
        octree_deallocate(self.c_octree)
        free(self.c_octree)

    @property
    def bound_particles(self):
        """
        The particle selection may select SPH particles with smoothing lengths
        which is in the octree domains. However, if the particle center is NOT
        in the octree, they are not included.

        So the number of particles passed to the tree *may* not be equal to the
        number of particles which are bound by the tree.
        """
        return self.c_octree.pend[0]

    @property
    def num_nodes(self):
        """
        The total number of nodes after tree construction
        """
        return self.c_octree.num_nodes

    @property
    def node_positions(self):
        """
        The centre of every node within the octree
        """
        cdef np.npy_intp shape[2]
        shape[0] = <np.npy_intp> self.c_octree.num_nodes
        shape[1] = 3
        arr = np.PyArray_SimpleNewFromData(
            2, &shape[0], np.NPY_FLOAT64, <void *>self.c_octree.node_positions)
        return np.copy(arr)

    @property
    def node_refined(self):
        """
        An array of length num_nodes which contains either True / False for
        whether each cell has refined or not. E.g False for a leaf, True for a
        node
        """
        cdef np.npy_intp shape
        shape = <np.npy_intp> self.c_octree.num_nodes
        arr = np.PyArray_SimpleNewFromData(
            1, &shape, np.NPY_UINT8, <void *>self.c_octree.refined)
        return np.copy(arr).astype(np.bool_)

    @property
    def node_depth(self):
        """
        The depth for each node in the tree. The root node is defined as a
        depth of 0.
        """
        cdef np.npy_intp shape = <np.npy_intp> self.c_octree.num_nodes
        arr = np.PyArray_SimpleNewFromData(
            1, &shape, np.NPY_UINT8, <void *>self.c_octree.depth)
        return np.copy(arr)

    @property
    def node_sizes(self):
        """
        The size of each node in the x, y and z directions. We calculate this
        on the fly. As we know the size of the whole tree and the depth of each
        node
        """
        cdef np.int64_t i
        sizes = np.zeros((self.c_octree.num_nodes, 3), dtype=np.float64)
        sizes[:, 0] = self.c_octree.size[0]
        sizes[:, 1] = self.c_octree.size[1]
        sizes[:, 2] = self.c_octree.size[2]

        for i in range(self.c_octree.num_nodes):
            sizes[i, :] /= <np.float64_t> (2.0**self.c_octree.depth[i] / 2.0)

        return sizes

    def _make_root(self, left_edge, right_edge):
        """
        The root node is the hardest node to build as we need to find out which
        particles we contain, and then sieving them into children is easy.

        In the case that the left_edge/right_edge is not defined then we select
        a tree that is sufficiently big to contain every particle.

        However, if they are defined then we need to loop through and find out
        which particles are actually in the tree.

        Parameters
        ----------
        {left,right}_edge : ndarray
            xyz coordinates of the lower left (upper right) corner of the octree.
            If None, the tree will be made large enough to encompass all particles.
        """
        cdef int i = 0

        # How many particles are there?
        self.c_octree.num_particles = self.input_positions.shape[0]

        # We now number all of the the particles in the tree. This allows us
        # to shuffle the pids and say for example Oct11 contains particles
        # 7 to 11
        # This pidx[7:12] would give the input indices of the particles we
        # store. This proxy allows us to re-arrange the particles without
        # re-arranging the users input data.
        self.c_octree.pidx = <np.int64_t *> malloc(
            self.c_octree.num_particles * sizeof(np.int64_t)
        )
        for i in range(0, self.c_octree.num_particles):
            self.c_octree.pidx[i] = i

        if left_edge is None:
            # If the edges are None, then we can just find the loop through
            # and find them out
            left_edge = np.zeros(3, dtype=np.float64)
            right_edge = np.zeros(3, dtype=np.float64)

            left_edge[0] = self.c_octree.pposx[0]
            left_edge[1] = self.c_octree.pposy[0]
            left_edge[2] = self.c_octree.pposz[0]
            right_edge[0] = self.c_octree.pposx[0]
            right_edge[1] = self.c_octree.pposy[0]
            right_edge[2] = self.c_octree.pposz[0]

            for i in range(self.c_octree.num_particles):
                left_edge[0] = min(self.c_octree.pposx[i], left_edge[0])
                left_edge[1] = min(self.c_octree.pposy[i], left_edge[1])
                left_edge[2] = min(self.c_octree.pposz[i], left_edge[2])

                right_edge[0] = max(self.c_octree.pposx[i], right_edge[0])
                right_edge[1] = max(self.c_octree.pposy[i], right_edge[1])
                right_edge[2] = max(self.c_octree.pposz[i], right_edge[2])

            self.c_octree.pstart[0] = 0
            self.c_octree.pend[0] = self.input_positions.shape[0]
        else:
            # Damn! The user did supply a left and right so we need to find
            # which particles are in the range
            left_edge = left_edge.astype(np.float64)
            right_edge = right_edge.astype(np.float64)

            # We loop through the particles and arrange them such that particles
            # in the tree are to the left of the split and the particles not
            # are to the right
            # e.g.
            # pidx = [1, 2, 3, 5 | 0, 4]
            # where split = 4 and particles 0 and 4 are not in the tree
            split = select(
                self.c_octree, left_edge, right_edge, 0, self.input_positions.shape[0])
            self.c_octree.pstart[0] = 0
            self.c_octree.pend[0] = split

        # Set the total size of the tree
        size = (right_edge - left_edge) / 2.0
        center = (right_edge + left_edge) / 2.0
        self.left_edge = left_edge
        self.right_edge = right_edge

        # Now we add the data about the root node!
        self.c_octree.node_positions[0] = center[0]
        self.c_octree.node_positions[1] = center[1]
        self.c_octree.node_positions[2] = center[2]

        self.c_octree.size[0] = size[0]
        self.c_octree.size[1] = size[1]
        self.c_octree.size[2] = size[2]

        # We are not refined yet
        self.c_octree.refined[0] = 0
        self.c_octree.depth[0] = 0

    def _allocate_octree(self):
        """
        This actually allocates the C struct Octree
        """
        self.c_octree = <Octree*> malloc(sizeof(Octree))
        self.c_octree.n_ref = self.n_ref
        self.c_octree.num_nodes = 1

        # This is sort of an arbitrary guess but it doesn't matter because
        # we will increase this value and attempt to reallocate if it is too
        # small
        self.c_octree.max_num_nodes = max(self.input_positions.shape[0] / self.n_ref, 1)
        self.c_octree.max_depth = self.max_depth

        # Fast C pointers to the particle coordinates
        self.c_octree.pposx = &self.input_positions[0, 0]
        self.c_octree.pposy = &self.input_positions[0, 1]
        self.c_octree.pposz = &self.input_positions[0, 2]

        if octree_allocate(self.c_octree, self.c_octree.max_num_nodes): return 1

    cdef void smooth_onto_cells(
        self,
        np.float64_t[:] buff,
        np.float64_t[:] buff_den,
        np.float64_t posx,
        np.float64_t posy,
        np.float64_t posz,
        np.float64_t hsml,
        np.float64_t prefactor,
        np.float64_t prefactor_norm,
        long int num_node,
        int use_normalization=0
    ):
        """
        We smooth a field onto cells within an octree using SPH deposition. To
        achieve this we loop through every oct in the tree and check if it is a
        leaf. A leaf just means that an oct has not refined, and thus has no
        children.

        Parameters
        ----------
        buff : memoryview
            The array which we are depositing the field onto, it has the
            length of the number of leaves.
        buff_den : memoryview
            The array we deposit just mass onto to allow normalization
        pos<> : float64_t
            The x, y, and z coordinates of the particle we are depositing
        hsml : float64_t
            The smoothing length of the particle
        prefactor(_norm) : float64_t
            This is a pre-computed value, based on the particles
                properties used in the deposition
        num_node : lonmg int
            The current node we are checking to see if refined or not
        use_normalization : int, default: 0
            Do we want a normalized sph field? If so, fill the buff_den.
        """

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
            # All direct children of the current node are the same size, thus
            # we can compute their size once, outside of the loop
            sz_factor = 1.0 / 2.0**(tree.depth[num_node] + 1)
            sqrt_sz_factor = math.sqrt(sz_factor)
            sx = tree.size[0]
            sy = tree.size[1]
            sz = tree.size[2]
            child_node_size = sqrt_sz_factor * math.sqrt(sx*sx + sy*sy + sz*sz)

            for i in range(8):
                child_node = tree.children[8*num_node + i]
                diff_x = tree.node_positions[3*child_node] - posx
                diff_y = tree.node_positions[3*child_node+1] - posy
                diff_z = tree.node_positions[3*child_node+2] - posz
                diff = math.sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z)

                # Could the current particle possibly intersect this child node?
                if diff - child_node_size - hsml < 0:
                    self.smooth_onto_cells(buff, buff_den, posx, posy, posz,
                                    hsml, prefactor, prefactor_norm,
                                    child_node, use_normalization=use_normalization)


    def interpolate_sph_cells(self,
        np.float64_t[:] buff,
        np.float64_t[:] buff_den,
        np.float64_t[:] posx,
        np.float64_t[:] posy,
        np.float64_t[:] posz,
        np.float64_t[:] pmass,
        np.float64_t[:] pdens,
        np.float64_t[:] hsml,
        np.float64_t[:] field,
        kernel_name="cubic",
        int use_normalization=0
    ):
        """
        We loop through every particle in the simulation and begin to deposit
        the particle properties onto all of the leaves that it intersects

        Parameters
        ----------
        buff : memoryview
            The array which we are depositing the field onto, it has the
            length of the number of leaves.
        buff_den : memoryview
            The array we deposit just mass onto to allow normalization
        pos<> : memoryview
            The x, y, and z coordinates of all the partciles
        pmass : memoryview
            The mass of the particles
        pdens : memoryview
            The density of the particles
        hsml : memoryview
            The smoothing lengths of the particles
        field : memoryview
            The field we are depositing for each particle
        kernel_name: str, default: "cubic"
            Choice of kernel for SPH deposition
        use_normalization : int, default: 0
            Do we want a normalized sph field? If so, fill the buff_den.
        """
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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t separate(
        np.float64_t * pos,
        np.int64_t * pidx,
        double value,
        np.int64_t start,
        np.int64_t end
    ) nogil:
    """
    This is a simple utility function which takes a selection of particles and
    re-arranges them such that values below `value` are to the left of split and
    values above are to the right.

    Parameters
    ----------
    pos : float64_t *
        Pointer to the coordinates we are splitting along
    pidx : int64_t &
        Pointer to the corresponding particle IDs
    value : double
        The value to split the data along
    start : int64_t
        Index of first particle in the current node
    end : int64_t
        Index of the last particle in the current node
    """
    cdef np.int64_t index
    cdef np.int64_t split = start
    for index in range(start, end):
        if pos[pidx[index]] < value:
            pidx[split], pidx[index] = pidx[index], pidx[split]
            split +=  1

    return split


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void split_helper(Octree * tree, np.int64_t node_idx, np.int64_t * splits):
    """
    A utility function to split a collection of particles along the x, y and z
    direction such that we identify the particles within the 8 children of an
    oct.

    We first split the particles in the x direction, then we split the particles
    in the y and the z directions. This is currently hardcoded for 8 octs but
    can be readily extended to allow over-refinement.

    The splits[0] and splits[1] tell oct one the start and last particle
    The splits[1] and splits[2] tell oct two the start and last particle and so
    on

    Parameters
    ----------
    tree : Octree *
        A pointer to the octree
    node_idx :  int64_t
        The index of the oct we are splitting
    splits : int64_t *
        Pointer to split array which stores the start and end indices. It needs
        to be N+1 long where N is the number of children
    """
    splits[0] = tree.pstart[node_idx]
    splits[8] = tree.pend[node_idx]

    splits[4] = separate(
        tree.pposx, tree.pidx, tree.node_positions[(3*node_idx)], splits[0], splits[8])

    splits[2] = separate(
        tree.pposy, tree.pidx, tree.node_positions[(3*node_idx)+1], splits[0], splits[4])
    splits[6] = separate(
        tree.pposy, tree.pidx, tree.node_positions[(3*node_idx)+1], splits[4], splits[8])

    splits[1] = separate(
        tree.pposz, tree.pidx, tree.node_positions[(3*node_idx)+2], splits[0], splits[2])
    splits[3] = separate(
        tree.pposz, tree.pidx, tree.node_positions[(3*node_idx)+2], splits[2], splits[4])
    splits[5] = separate(
        tree.pposz, tree.pidx, tree.node_positions[(3*node_idx)+2], splits[4], splits[6])
    splits[7] = separate(
        tree.pposz, tree.pidx, tree.node_positions[(3*node_idx)+2], splits[6], splits[8])


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.int64_t select(
        Octree * octree,
        np.float64_t[::1] left_edge,
        np.float64_t[::1] right_edge,
        np.int64_t start,
        np.int64_t end
    ) nogil:
    """
    Re-arrange the input particles such that those outside the bounds of the
    tree occur after the split index and can thus be ignored for the remainder
    of the tree construction

    Parameters
    ----------
    tree : Octree *
        A pointer to the octree
    left_edge : ndarray
        The coords of the lower left corner of the box
    right_edge : ndarray
        The coords of the upper right corner of the box
    start : int64_t
        The first particle in the bounds (0)
    end : int64_t
        The last particle in the bounds
    """
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
                    split += 1

    return split
