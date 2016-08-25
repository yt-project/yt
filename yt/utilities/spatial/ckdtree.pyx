# Copyright Anne M. Archibald 2008
# Released under the scipy license
import numpy as np
cimport numpy as np
cimport libc.stdlib as stdlib
cimport cython

import kdtree

cdef extern from "platform_dep.h":
    # NOTE that size_t might not be int
    void *alloca(int)

cdef np.float64_t infinity = np.inf

__all__ = ['cKDTree']

# priority queue
cdef union heapcontents:
    int intdata
    char* ptrdata

cdef struct heapitem:
    np.float64_t priority
    heapcontents contents

cdef struct heap:
    int n
    heapitem* heap
    int space

cdef inline heapcreate(heap* self,int initial_size):
    self.space = initial_size
    self.heap = <heapitem*>stdlib.malloc(sizeof(heapitem)*self.space)
    self.n=0

cdef inline heapdestroy(heap* self):
    stdlib.free(self.heap)

cdef inline heapresize(heap* self, int new_space):
    if new_space<self.n:
        raise ValueError("Heap containing %d items cannot be resized to %d" % (self.n, new_space))
    self.space = new_space
    self.heap = <heapitem*>stdlib.realloc(<void*>self.heap,new_space*sizeof(heapitem))

cdef inline heappush(heap* self, heapitem item):
    cdef int i
    cdef heapitem t

    self.n += 1
    if self.n>self.space:
        heapresize(self,2*self.space+1)

    i = self.n-1
    self.heap[i] = item
    while i>0 and self.heap[i].priority<self.heap[(i-1)//2].priority:
        t = self.heap[(i-1)//2]
        self.heap[(i-1)//2] = self.heap[i]
        self.heap[i] = t
        i = (i-1)//2

cdef heapitem heappeek(heap* self):
    return self.heap[0]

cdef heapremove(heap* self):
    cdef heapitem t
    cdef int i, j, k, l

    self.heap[0] = self.heap[self.n-1]
    self.n -= 1
    if self.n < self.space//4 and self.space>40: #FIXME: magic number
        heapresize(self,self.space//2+1)

    i=0
    j=1
    k=2
    while ((j<self.n and
                self.heap[i].priority > self.heap[j].priority or
            k<self.n and
                self.heap[i].priority > self.heap[k].priority)):
        if k<self.n and self.heap[j].priority>self.heap[k].priority:
            l = k
        else:
            l = j
        t = self.heap[l]
        self.heap[l] = self.heap[i]
        self.heap[i] = t
        i = l
        j = 2*i+1
        k = 2*i+2

cdef heapitem heappop(heap* self):
    cdef heapitem it
    it = heappeek(self)
    heapremove(self)
    return it





# utility functions
cdef inline np.float64_t dmax(np.float64_t x, np.float64_t y):
    if x>y:
        return x
    else:
        return y
cdef inline np.float64_t dabs(np.float64_t x):
    if x>0:
        return x
    else:
        return -x
cdef inline np.float64_t dmin(np.float64_t x, np.float64_t y):
    if x<y:
        return x
    else:
        return y
cdef inline np.float64_t _distance_p(np.float64_t*x,np.float64_t*y,np.float64_t p,int k,np.float64_t upperbound,
    np.float64_t*period):
    """Compute the distance between x and y

    Computes the Minkowski p-distance to the power p between two points.
    If the distance**p is larger than upperbound, then any number larger
    than upperbound may be returned (the calculation is truncated).

    Periodicity added by S. Skory.
    """
    cdef int i
    cdef np.float64_t r, m
    r = 0
    if p==infinity:
        for i in range(k):
            m = dmin(dabs(x[i] - y[i]), period[i] - dabs(x[i] - y[i]))
            r = dmax(r,m)
            if r>upperbound:
                return r
    elif p==1:
        for i in range(k):
            m = dmin(dabs(x[i] - y[i]), period[i] - dabs(x[i] - y[i]))
            r += m
            if r>upperbound:
                return r
    elif p==2:
        for i in range(k):
            m = dmin(dabs(x[i] - y[i]), period[i] - dabs(x[i] - y[i]))
            r += m*m
            if r>upperbound:
                return r
    else:
        for i in range(k):
            m = dmin(dabs(x[i] - y[i]), period[i] - dabs(x[i] - y[i]))
            r += m**p
            if r>upperbound:
                return r
    return r



# Tree structure
cdef struct innernode:
    int split_dim
    int n_points
    np.float64_t split
    np.float64_t* maxes
    np.float64_t* mins
    innernode* less
    innernode* greater
cdef struct leafnode:
    int split_dim
    int n_points
    int start_idx
    int end_idx
    np.float64_t* maxes
    np.float64_t* mins

# this is the standard trick for variable-size arrays:
# malloc sizeof(nodeinfo)+self.m*sizeof(np.float64_t) bytes.
cdef struct nodeinfo:
    innernode* node
    np.float64_t side_distances[0]

cdef class cKDTree:
    """kd-tree for quick nearest-neighbor lookup

    This class provides an index into a set of k-dimensional points
    which can be used to rapidly look up the nearest neighbors of any
    point.

    The algorithm used is described in Maneewongvatana and Mount 1999.
    The general idea is that the kd-tree is a binary trie, each of whose
    nodes represents an axis-aligned hyperrectangle. Each node specifies
    an axis and splits the set of points based on whether their coordinate
    along that axis is greater than or less than a particular value.

    During construction, the axis and splitting point are chosen by the
    "sliding midpoint" rule, which ensures that the cells do not all
    become long and thin.

    The tree can be queried for the r closest neighbors of any given point
    (optionally returning only those within some maximum distance of the
    point). It can also be queried, with a substantial gain in efficiency,
    for the r approximate closest neighbors.

    For large dimensions (20 is already large) do not expect this to run
    significantly faster than brute force. High-dimensional nearest-neighbor
    queries are a substantial open problem in computer science.

    Parameters
    ----------
    data : array-like, shape (n,m)
        The n data points of dimension m to be indexed. This array is
        not copied unless this is necessary to produce a contiguous
        array of np.float64_ts, and so modifying this data will result in
        bogus results.
    leafsize : positive integer
        The number of points at which the algorithm switches over to
        brute-force.

    """

    cdef innernode* tree
    cdef readonly object data
    cdef np.float64_t* raw_data
    cdef readonly int n, m
    cdef readonly int leafsize
    cdef readonly object maxes
    cdef np.float64_t* raw_maxes
    cdef readonly object mins
    cdef np.float64_t* raw_mins
    cdef object indices
    cdef np.int64_t* raw_indices
    def __init__(cKDTree self, data, int leafsize=10):
        cdef np.ndarray[np.float64_t, ndim=2] inner_data
        cdef np.ndarray[np.float64_t, ndim=1] inner_maxes
        cdef np.ndarray[np.float64_t, ndim=1] inner_mins
        cdef np.ndarray[np.int64_t, ndim=1] inner_indices
        self.data = np.ascontiguousarray(data,dtype="float64")
        self.n, self.m = np.shape(self.data)
        self.leafsize = leafsize
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.ascontiguousarray(np.amax(self.data,axis=0))
        self.mins = np.ascontiguousarray(np.amin(self.data,axis=0))
        self.indices = np.ascontiguousarray(np.arange(self.n,dtype=np.int64))

        inner_data = self.data
        self.raw_data = <np.float64_t*>inner_data.data
        inner_maxes = self.maxes
        self.raw_maxes = <np.float64_t*>inner_maxes.data
        inner_mins = self.mins
        self.raw_mins = <np.float64_t*>inner_mins.data
        inner_indices = self.indices
        self.raw_indices = <np.int64_t*>inner_indices.data

        self.tree = self.__build(0, self.n, self.raw_maxes, self.raw_mins)

    cdef innernode* __build(cKDTree self, int start_idx, int end_idx, np.float64_t* maxes, np.float64_t* mins):
        cdef leafnode* n
        cdef innernode* ni
        cdef int i, j, t, p, q, d
        cdef np.float64_t size, split, minval, maxval
        cdef np.float64_t*mids
        if end_idx-start_idx<=self.leafsize:
            n = <leafnode*>stdlib.malloc(sizeof(leafnode))
            # Skory
            n.maxes = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
            n.mins = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
            for i in range(self.m):
                n.maxes[i] = maxes[i]
                n.mins[i] = mins[i]
            n.split_dim = -1
            n.start_idx = start_idx
            n.end_idx = end_idx
            return <innernode*>n
        else:
            d = 0
            size = 0
            for i in range(self.m):
                if maxes[i]-mins[i] > size:
                    d = i
                    size =  maxes[i]-mins[i]
            maxval = maxes[d]
            minval = mins[d]
            if maxval==minval:
                # all points are identical; warn user?
                n = <leafnode*>stdlib.malloc(sizeof(leafnode))
                n.split_dim = -1
                n.start_idx = start_idx
                n.end_idx = end_idx
                return <innernode*>n

            split = (maxval+minval)/2

            p = start_idx
            q = end_idx-1
            while p<=q:
                if self.raw_data[self.raw_indices[p]*self.m+d]<split:
                    p+=1
                elif self.raw_data[self.raw_indices[q]*self.m+d]>=split:
                    q-=1
                else:
                    t = self.raw_indices[p]
                    self.raw_indices[p] = self.raw_indices[q]
                    self.raw_indices[q] = t
                    p+=1
                    q-=1

            # slide midpoint if necessary
            if p==start_idx:
                # no points less than split
                j = start_idx
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx+1, end_idx):
                    if self.raw_data[self.raw_indices[i]*self.m+d]<split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[start_idx]
                self.raw_indices[start_idx] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = start_idx+1
                q = start_idx
            elif p==end_idx:
                # no points greater than split
                j = end_idx-1
                split = self.raw_data[self.raw_indices[j]*self.m+d]
                for i in range(start_idx, end_idx-1):
                    if self.raw_data[self.raw_indices[i]*self.m+d]>split:
                        j = i
                        split = self.raw_data[self.raw_indices[j]*self.m+d]
                t = self.raw_indices[end_idx-1]
                self.raw_indices[end_idx-1] = self.raw_indices[j]
                self.raw_indices[j] = t
                p = end_idx-1
                q = end_idx-2

            # construct new node representation
            ni = <innernode*>stdlib.malloc(sizeof(innernode))

            mids = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
            for i in range(self.m):
                mids[i] = maxes[i]
            mids[d] = split
            ni.less = self.__build(start_idx,p,mids,mins)

            for i in range(self.m):
                mids[i] = mins[i]
            mids[d] = split
            ni.greater = self.__build(p,end_idx,maxes,mids)

            stdlib.free(mids)

            ni.split_dim = d
            ni.split = split
            # Skory
            ni.maxes = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
            ni.mins = <np.float64_t*>stdlib.malloc(sizeof(np.float64_t)*self.m)
            for i in range(self.m):
                ni.maxes[i] = maxes[i]
                ni.mins[i] = mins[i]

            return ni

    cdef __free_tree(cKDTree self, innernode* node):
        if node.split_dim!=-1:
            self.__free_tree(node.less)
            self.__free_tree(node.greater)
        stdlib.free(node.maxes) # Skory
        stdlib.free(node.mins)
        stdlib.free(node)

    def __dealloc__(cKDTree self):
        if <int>(self.tree) == 0:
            # should happen only if __init__ was never called
            return
        self.__free_tree(self.tree)

    cdef void __query(cKDTree self,
            np.float64_t*result_distances,
            np.int64_t*result_indices,
            np.float64_t*x,
            int k,
            np.float64_t eps,
            np.float64_t p,
            np.float64_t distance_upper_bound,
            np.float64_t*period):
        assert(p == 2)
        assert(eps == 0.0)
        assert(distance_upper_bound == infinity)
        cdef heap q
        cdef heap neighbors

        cdef int i, j, i2, j2
        cdef np.float64_t t, y
        cdef nodeinfo* inf
        cdef nodeinfo* inf2
        cdef np.float64_t d, di
        cdef np.float64_t m_left, m_right, m
        cdef np.float64_t epsfac
        cdef np.float64_t min_distance
        cdef np.float64_t far_min_distance
        cdef heapitem it, it2, neighbor
        cdef leafnode* node
        cdef innernode* inode
        cdef innernode* near
        cdef innernode* far
        cdef np.float64_t* side_distances

        # priority queue for chasing nodes
        # entries are:
        #  minimum distance between the cell and the target
        #  distances between the nearest side of the cell and the target
        #  the head node of the cell
        heapcreate(&q,12)

        # priority queue for the nearest neighbors
        # furthest known neighbor first
        # entries are (-distance**p, i)
        heapcreate(&neighbors,k)

        # set up first nodeinfo
        inf = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(np.float64_t))
        inf.node = self.tree
        for i in range(self.m):
            inf.side_distances[i] = 0
            t = x[i]-self.raw_maxes[i]
            if t>inf.side_distances[i]:
                inf.side_distances[i] = t
            else:
                t = self.raw_mins[i]-x[i]
                if t>inf.side_distances[i]:
                    inf.side_distances[i] = t
            inf.side_distances[i]=inf.side_distances[i]*inf.side_distances[i]

        # compute first distance
        min_distance = 0.
        for i in range(self.m):
            min_distance += inf.side_distances[i]

        # fiddle approximation factor
        epsfac=1

        while True:
            if inf.node.split_dim==-1:
                node = <leafnode*>inf.node

                # brute-force
                for i in range(node.start_idx,node.end_idx):
                    d = 0.0
                    for i2 in range(self.m):
                        y = self.raw_data[self.raw_indices[i]*self.m + i2]
                        di = dmin(dabs(x[i2] - y), period[i2] - dabs(x[i2] - y))
                        d += di*di
                    if d<distance_upper_bound:
                        # replace furthest neighbor
                        if neighbors.n==k:
                            heapremove(&neighbors)
                        neighbor.priority = -d
                        neighbor.contents.intdata = self.raw_indices[i]
                        heappush(&neighbors,neighbor)

                        # adjust upper bound for efficiency
                        if neighbors.n==k:
                            distance_upper_bound = -heappeek(&neighbors).priority
                # done with this node, get another
                stdlib.free(inf)
                if q.n==0:
                    # no more nodes to visit
                    break
                else:
                    it = heappop(&q)
                    inf = <nodeinfo*>it.contents.ptrdata
                    min_distance = it.priority
            else:
                inode = <innernode*>inf.node

                # we don't push cells that are too far onto the queue at all,
                # but since the distance_upper_bound decreases, we might get
                # here even if the cell's too far
                if min_distance>distance_upper_bound*epsfac:
                    # since this is the nearest cell, we're done, bail out
                    stdlib.free(inf)
                    # free all the nodes still on the heap
                    for i in range(q.n):
                        stdlib.free(q.heap[i].contents.ptrdata)
                    break

                # set up children for searching
                if x[inode.split_dim]<inode.split:
                    near = inode.less
                    far = inode.greater
                else:
                    near = inode.greater
                    far = inode.less

                # near child is at the same distance as the current node
                # we're going here next, so no point pushing it on the queue
                # no need to recompute the distance or the side_distances
                inf.node = near

                # far child is further by an amount depending only
                # on the split value; compute its distance and side_distances
                # and push it on the queue if it's near enough
                inf2 = <nodeinfo*>stdlib.malloc(sizeof(nodeinfo)+self.m*sizeof(np.float64_t))
                it2.contents.ptrdata = <char*> inf2
                inf2.node = far

                # Periodicity added by S Skory
                m_left = dmin( dabs(far.mins[inode.split_dim] - x[inode.split_dim]), \
                    period[inode.split_dim] -  dabs(far.mins[inode.split_dim] - x[inode.split_dim]))
                m_right = dmin( dabs(far.maxes[inode.split_dim] - x[inode.split_dim]), \
                    period[inode.split_dim] -  dabs(far.maxes[inode.split_dim] - x[inode.split_dim]))
                m = dmin(m_left,m_right)

                # most side distances unchanged
                for i in range(self.m):
                    inf2.side_distances[i] = inf.side_distances[i]

                # one side distance changes
                # we can adjust the minimum distance without recomputing
                inf2.side_distances[inode.split_dim] = m*m
                #far_min_distance = min_distance - inf.side_distances[inode.split_dim] + inf2.side_distances[inode.split_dim]
                far_min_distance = m*m

                it2.priority = far_min_distance


                # far child might be too far, if so, don't bother pushing it
                if far_min_distance<=distance_upper_bound*epsfac:
                    heappush(&q,it2)
                else:
                    stdlib.free(inf2)
                    # just in case
                    it2.contents.ptrdata = <char*> 0

        # fill output arrays with sorted neighbors
        for i in range(neighbors.n-1,-1,-1):
            neighbor = heappop(&neighbors) # FIXME: neighbors may be realloced
            result_indices[i] = neighbor.contents.intdata
            result_distances[i] = (-neighbor.priority) #**(1./p) S. Skory

        heapdestroy(&q)
        heapdestroy(&neighbors)

    def query(cKDTree self, object x, int k=1, np.float64_t eps=0, np.float64_t p=2,
            np.float64_t distance_upper_bound=infinity, object period=None):
        """query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf,
           period=None)

        Query the kd-tree for nearest neighbors.

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : int
            The number of nearest neighbors to return.
        eps : non-negative float
            Return approximate nearest neighbors; the k-th returned value
            is guaranteed to be no further than (1 + `eps`) times the
            distance to the real k-th nearest neighbor.
        p : float, 1 <= p <= infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance.
            2 is the usual Euclidean distance.
            infinity is the maximum-coordinate-difference distance.
        distance_upper_bound : non-negative float
            Return only neighbors within this distance.  This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns
        -------
        d : ndarray of floats
            The distances to the nearest neighbors.
            If `x` has shape tuple+(self.m,), then `d` has shape tuple+(k,).
            Missing neighbors are indicated with infinite distances.
        i : ndarray of ints
            The locations of the neighbors in self.data.
            If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
            Missing neighbors are indicated with self.n.

        """
        cdef np.ndarray[np.int64_t, ndim=2] ii
        cdef np.ndarray[np.float64_t, ndim=2] dd
        cdef np.ndarray[np.float64_t, ndim=2] xx
        cdef np.ndarray[np.float64_t, ndim=1] cperiod
        cdef int c
        x = np.asarray(x).astype("float64")
        if period is None:
            period = np.array([np.inf]*self.m)
        else:
            period = np.asarray(period).astype("float64")
        cperiod = np.ascontiguousarray(period)
        if np.shape(x)[-1] != self.m:
            raise ValueError("x must consist of vectors of length %d but has shape %s" % (self.m, np.shape(x)))
        if p<1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        if len(x.shape)==1:
            single = True
            x = x[np.newaxis,:]
        else:
            single = False
        retshape = np.shape(x)[:-1]
        n = np.prod(retshape)
        xx = np.reshape(x,(n,self.m))
        xx = np.ascontiguousarray(xx)
        dd = np.empty((n,k),dtype="float64")
        dd.fill(infinity)
        ii = np.empty((n,k),dtype="int64")
        ii.fill(self.n)
        for c in range(n):
            self.__query(
                    (<np.float64_t*>dd.data)+c*k,
                    (<np.int64_t*>ii.data)+c*k,
                    (<np.float64_t*>xx.data)+c*self.m,
                    k,
                    eps,
                    p,
                    distance_upper_bound,
                    <np.float64_t*>cperiod.data)
        if single:
            if k==1:
                return dd[0,0], ii[0,0]
            else:
                return dd[0], ii[0]
        else:
            if k==1:
                return np.reshape(dd[...,0],retshape), np.reshape(ii[...,0],retshape)
            else:
                return np.reshape(dd,retshape+(k,)), np.reshape(ii,retshape+(k,))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def chainHOP_get_dens(cKDTree self, object omass, int num_neighbors=65, \
            int nMerge=6):
        """ query the tree for the nearest neighbors, to get the density
            of particles for chainHOP.

        Parameters:
        ===========

        mass: A array-like list of the masses of the particles, in the same
            order as the data that went into building the kd tree.

        num_neighbors: Optional, the number of neighbors to search for and to
            use in the density calculation. Default is 65, and is probably what
            one should stick with.

        nMerge: The number of nearest neighbor tags to return for each particle.

        Returns:
        ========

        dens: An array of the densities for each particle, in the same order
            as the input data.

        tags: A two-dimensional array of the indexes, nMerge nearest neighbors
            for each particle.

        """

        # We're no np.int64_ter returning all the tags in this step.
        # We do it chunked, in find_chunk_nearest_neighbors.
        #cdef np.ndarray[np.int64_t, ndim=2] tags
        cdef np.ndarray[np.float64_t, ndim=1] dens
        cdef int i, pj, j
        cdef np.float64_t ih2, fNorm, r2, rs

        #tags = np.empty((self.n, nMerge), dtype="int64")
        dens = np.zeros(self.n, dtype="float64")
        cdef np.ndarray[np.float64_t, ndim=2] local_data = self.data

        cdef np.ndarray[np.float64_t, ndim=1] mass = np.array(omass).astype("float64")
        cdef np.float64_t ipi = 1.0/np.pi

        cdef np.float64_t *query = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * self.m)
        cdef np.float64_t *dist_temp = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * num_neighbors)
        cdef np.int64_t *tags_temp = <np.int64_t *> alloca(
                    sizeof(np.int64_t) * num_neighbors)
        cdef np.float64_t period[3]
        for i in range(3): period[i] = 1.0

        for i in range(self.n):
            for j in range(self.m):
                query[j] = local_data[i,j]
            self.__query(dist_temp, tags_temp,
                         query, num_neighbors, 0.0,
                         2, infinity, period)

            #calculate the density for this particle
            ih2 = -1
            for j in range(num_neighbors):
                ih2 = dmax(ih2, dist_temp[j])
            ih2 = 4.0/ih2
            fNorm = 0.5*(ih2**1.5)*ipi
            for j in range(num_neighbors):
                pj = tags_temp[j]
                r2 = dist_temp[j] * ih2
                rs = 2.0 - (r2**0.5)
                if (r2 < 1.0):
                    rs = (1.0 - 0.75*rs*r2)
                else:
                    rs = 0.25*rs*rs*rs
                rs = rs * fNorm
                dens[i] = dens[i] + rs * mass[pj]
                dens[pj] = dens[pj] + rs * mass[i]

            # store nMerge nearest neighbors
            #tags[i,:] = tags_temp[:nMerge]

        #return (dens, tags)
        return dens

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def find_chunk_nearest_neighbors(cKDTree self, int start, int finish, \
        int num_neighbors=65):
        """ query the tree in chunks, between start and finish, recording the
            nearest neighbors.

        Parameters:
        ===========

        start: The starting point in the dataset for this search.

        finish: The ending point in the dataset for this search.

        num_neighbors: Optional, the number of neighbors to search for.
            The default is 65.

        Returns:
        ========

        chunk_tags: A two-dimensional array of the nearest neighbor tags for the
            points in this search.

        """

        cdef np.ndarray[np.int64_t, ndim=2] chunk_tags
        cdef np.ndarray[np.float64_t, ndim=2] local_data = self.data
        cdef int i, j

        chunk_tags = np.empty((finish-start, num_neighbors), dtype="int64")
        cdef np.float64_t *query = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * self.m)
        cdef np.float64_t *dist_temp = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * num_neighbors)
        cdef np.int64_t *tags_temp = <np.int64_t *> alloca(
                    sizeof(np.int64_t) * num_neighbors)
        cdef np.float64_t period[3]
        for i in range(3): period[i] = 1.0

        for i in range(finish-start):
            for j in range(self.m):
                query[j] = local_data[i+start,j]
            self.__query(dist_temp, tags_temp,
                         query, num_neighbors, 0.0,
                         2, infinity, period)
            for j in range(num_neighbors):
                chunk_tags[i,j] = tags_temp[j]

        return chunk_tags

    def chainHOP_preconnect(self, np.ndarray[np.int64_t, ndim=1] chainID,
                                  np.ndarray[np.float64_t, ndim=1] density,
                                  np.ndarray[np.float64_t, ndim=1] densest_in_chain,
                                  np.ndarray bis_inside,
                                  np.ndarray bsearch_again,
                                  np.float64_t peakthresh,
                                  np.float64_t saddlethresh,
                                  int nn, int nMerge,
                                  object chain_map):
        cdef np.ndarray[np.int32_t, ndim=1] is_inside
        cdef np.ndarray[np.int32_t, ndim=1] search_again
        cdef np.ndarray[np.float64_t, ndim=2] pos
        cdef np.int64_t thisNN, thisNN_chainID, same_count
        cdef np.float64_t *query = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * self.m)
        cdef np.float64_t *dist_temp = <np.float64_t *> alloca(
                    sizeof(np.float64_t) * nn)
        cdef np.int64_t *tags_temp = <np.int64_t *> alloca(
                    sizeof(np.int64_t) * nn)
        cdef np.float64_t period[3]
        cdef np.float64_t thisNN_max_dens, boundary_density
        cdef int i, j, npart, chainID_i, part_mas_dens
        is_inside = bis_inside.astype("int32")
        search_again = bsearch_again.astype("int32")
        pos = self.data
        npart = pos.shape[0]
        for i in range(3): period[i] = 1.0
        for i in xrange(npart):
            # Don't consider this particle if it's not part of a chain.
            if chainID[i] < 0: continue
            chainID_i = chainID[i]
            # If this particle is in the padding, don't make a connection.
            if not is_inside[i]: continue
            # Find this particle's chain max_dens.
            part_max_dens = densest_in_chain[chainID_i]
            # We're only connecting >= peakthresh chains now.
            if part_max_dens < peakthresh: continue
            # Loop over nMerge closest nearest neighbors.
            for j in range(self.m):
                query[j] = pos[i,j]
            self.__query(dist_temp, tags_temp,
                         query, nn, 0.0,
                         2, infinity, period)
            same_count = 0
            for j in xrange(int(nMerge+1)):
                thisNN = tags_temp[j+1] # Don't consider ourselves at tags_temp[0]
                thisNN_chainID = chainID[thisNN]
                # If our neighbor is in the same chain, move on.
                # Move on if these chains are already connected:
                if chainID_i == thisNN_chainID or \
                        thisNN_chainID in chain_map[chainID_i]:
                    same_count += 1
                    continue
                # Everything immediately below is for
                # neighboring particles with a chainID.
                if thisNN_chainID >= 0:
                    # Find thisNN's chain's max_dens.
                    thisNN_max_dens = densest_in_chain[thisNN_chainID]
                    # We're only linking peakthresh chains
                    if thisNN_max_dens < peakthresh: continue
                    # Calculate the two groups boundary density.
                    boundary_density = (density[thisNN] + density[i]) / 2.
                    # Don't connect if the boundary is too low.
                    if boundary_density < saddlethresh: continue
                    # Mark these chains as related.
                    chain_map[thisNN_chainID].add(chainID_i)
                    chain_map[chainID_i].add(thisNN_chainID)
            if same_count == nMerge + 1:
                # All our neighbors are in the same chain already, so
                # we don't need to search again.
                search_again[i] = 0
        return search_again
