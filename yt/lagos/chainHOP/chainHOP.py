import math,sys #, cPickle
from bisect import insort
from sets import Set

from yt.lagos import *
from yt.extensions.kdtree import *
from mpi4py import MPI
from Forthon import *

class RunChainHOP(ParallelAnalysisInterface):
    def __init__(self,period, padding, num_neighbors, bounds,
            xpos, ypos, zpos, index, mass, threshold=160.0):
        self.threshold = threshold
        self.saddlethresh = 2.5 * threshold
        self.peakthresh = 3 * threshold
        self.period = period
        self.padding = padding
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.size = len(self.xpos)
        self.index = na.array(index, dtype='int64')
        self.mass = mass
        self.padded_particles = []
        self.nMerge = 4
        self.d = None
        self.chainID = na.ones(self.size,dtype='int64') * -1
        self._chain_hop()
    

    def _init_kd_tree(self):
        """
        Set up the data objects that get passed to the kD-tree code.
        """
        # Yes, we really do need to initialize this many arrays.
        # They're deleted in _chainHOP.
        fKD.nn_tags = empty((self.nMerge + 2,self.size), dtype='l')
        fKD.dens = zeros(self.size, dtype='d')
        fKD.mass = empty(self.size, dtype='d')
        fKD.pos = empty((3,self.size), dtype='d')
        fKD.qv = empty(3, dtype='d')
        fKD.nn = self.num_neighbors
        fKD.nMerge = self.nMerge + 2
        fKD.nparts = self.size
        fKD.sort = True # Slower, but needed in _connect_chains
        fKD.rearrange = True # Faster, but uses more memory
        # This actually copies the data into the fortran space.
        fKD.pos[0, :] = self.xpos
        fKD.pos[1, :] = self.ypos
        fKD.pos[2, :] = self.zpos
        fKD.mass = self.mass
        # Now call the fortran.
        create_tree()

    def _is_inside(self):
        """
        There are three classes of particles.
        1. Particles inside the 'real' region of each subvolume.
        2. Particles ouside, added in the 'padding' for purposes of having 
           correct particle densities in the real region.
        3. Particles that are one padding distance inside the edges of the
           real region. The chainIDs of these particles are communicated
           to the neighboring tasks so chains can be merged into groups.
        """
        # Test to see if the points are in the 'real' region
        (LE, RE) = self.bounds
        points = fKD.pos.transpose()
        self.is_inside = ( (points >= LE).all(axis=1) * \
            (points < RE).all(axis=1) )
        # Below we find out which particles are in the `annulus', one padding
        # distance inside the boundaries. First we find the particles outside
        # this inner boundary.
        LE += self.padding
        RE -= self.padding
        inner = na.invert( (points >= LE).all(axis=1) * \
            (points < RE).all(axis=1) )
        # After inverting the logic above, we want points that are both
        # inside the real region, but within one padding of the boundary,
        # and this will do it.
        self.is_inside_annulus = na.bitwise_and(self.is_inside, inner)
        # Below we make a mapping of real particle index->local ID
        # Unf. this has to be a dict, because any task can have
        # particles of any particle_index, which means that if it were an
        # array every task would probably end up having this array be as long
        # as the full number of particles.
        self.rev_index = {}
        for i in xrange(self.size):
            # Only padded and annulus particles are needed in the dict.
            if not self.is_inside[i] or self.is_inside_annulus[i]:
                self.rev_index[self.index[i]] = i
        del points, inner

    # These next two functions aren't currently being used. They figure out
    # which direction to send particle data from each task. They may come 
    # in handy in the future, if needed.

    def _find_shift_padded(self, i):
        """
        Find the shift vector for a padded particle.
        """
        # This shouldn't happen, but just in case, an internal particle
        # has no shift.
        if self.is_inside[i]: return na.array([0,0,0], dtype='int64')
        (LE, RE) = self.bounds
        xshift, yshift, zshift = 0, 0, 0
        if self.xpos[i] >= RE[0]:
            xshift = 1
            # this accounts for padded regions that haven't been repositioned
            # around the periodic boundaries
            if (self.xpos[i] - RE[0]) > self.padding:
                xshift = -1
        elif self.xpos[i] < LE[0]:
            xshift = -1
            if (self.xpos[i] - LE[0]) >= self.padding:
                xshift = 1
        if self.ypos[i] >= RE[1]:
            yshift = 1
            if (self.ypos[i] - RE[0]) > self.padding:
                yshift = -1
        elif self.ypos[i] < LE[1]:
            yshift = -1
            if (self.ypos[i] - LE[0]) >= self.padding:
                yshift = 1
        if self.zpos[i] >= RE[2]:
            zshift = 1
            if (self.zpos[i] - RE[0]) > self.padding:
                zshift = -1
        elif self.zpos[i] < LE[2]:
            zshift = -1
            if (self.zpos[i] - LE[0]) >= self.padding:
                zshift = 1
        shift = na.array([xshift,yshift,zshift],dtype='int64')
        return shift

    def _find_shift_real(self, i):
        """
        Find shift vectors for a particle in the 'real' region. A particle
        in the real region, close to the boundary, may be a padded particle
        in several neighbors, so we need to return multiple shift vectors.
        """
        return [ na.array([1,0,0], dtype='int64') ]
        # skip padded particles
        if not self.is_inside[i]: return na.array([0]*3, dtype='int64')
        # Adjust the boundaries
        (LE, RE) = self.bounds
        LE += na.array([self.padding]*3)
        RE -= na.array([self.padding]*3)
        xshift, yshift, zshift = 0, 0, 0
        if self.xpos[i] >= RE[0]:
            xshift = 1
        elif self.xpos[i] < LE[0]:
            xshift = -1
        if self.ypos[i] >= RE[1]:
            yshift = 1
        elif self.ypos[i] < LE[1]:
            yshift = -1
        if self.zpos[i] >= RE[2]:
            zshift = 1
        elif self.zpos[i] < LE[2]:
            zshift = -1
        count = abs(xshift) + abs(yshift) + abs(zshift)
        vectors = []
        base = na.array([xshift,yshift,zshift], dtype='int64')
        diag = na.array([[1,0,0], [0,1,0], [0,0,1]], dtype='int64')
        vectors.append(base)
        if count == 2:
            # find the zero element so we can skip it
            zero_e = na.argmax(na.abs(base) < 1)
            for i in xrange(3):
                if i == zero_e: continue
                vectors.append(diag[i] * base[i])
        elif count == 3:
            for i in xrange(3):
                # the singletons
                vectors.append(diag[i] * base[i])
                for j in xrange(3):
                    # the doubles, the triple is already added, it's 'base'
                    if i == j: continue
                    vectors.append( (diag[i] * base[i]) + (diag[j] * base[j]))
        return vectors
        
        
    def _translate_shift(self, i=None, shift=None):
        """
        Given an integer, return the corresponding shift, and given a shift
        return the integer. There may be a better way to do this.
        """
        if self.d is None:
            d = {0:[-1,-1,-1], 1:[-1,-1,0], 2:[-1,-1,1], 3:[-1,0,-1], 4:[-1,0,0],
                5:[-1,0,1], 6:[-1,1,-1], 7:[-1,1,0], 8:[-1,1,1,], 9:[0,-1,-1],
                10:[0,-1,0], 11:[0,-1,1], 12:[0,0,-1], 13:[0,0,0], 14:[0,0,1],
                15:[0,1,-1], 16:[0,1,0], 17:[0,1,1], 18:[1,-1,-1], 19:[1,-1,0],
                20:[1,-1,1], 21:[1,0,-1], 22:[1,0,0], 23:[1,0,1], 24:[1,1,-1],
                25:[1,1,0], 26:[1,1,1]}
            for key in d:
                d[key] = na.array(d[key], dtype='int64')
            self.d = d
        if i==None and shift==None: return None
        if i!=None:
            return self.d[i]
        if shift!=None:
            shift = na.array(shift)
            for key in self.d:
                if (self.d[key] == shift).all():
                    return key

    def _densestNN(self):
        """
        For all particles, find their densest nearest neighbor. It is done in
        chunks to keep the memory usage down.
        The first search of nearest neighbors (done earlier) did not return all 
        num_neighbor neighbors, so we need to do it again, but we're not
        keeping the all of this data, just using it.
        """
        self.densestNN = na.ones(self.size,dtype='l')
        # We find nearest neighbors in chunks.
        chunksize = 10000
        fKD.chunk_tags = empty((self.num_neighbors, chunksize), dtype='l')
        start = 1 # Fortran counting!
        finish = 0
        while finish < self.size:
            finish = min(finish+chunksize,self.size)
            #fKD.chunk_tags = empty((self.num_neighbors,finish-start+1), dtype='l')
            # Call the fortran. start and finish refer to the data locations
            # in fKD.pos, and specify the range of particles to find nearest
            # neighbors
            fKD.start = start
            fKD.finish = finish
            find_chunk_nearest_neighbors()
            chunk_NNtags = (fKD.chunk_tags[:,:finish-start+1] - 1).transpose()
            # Find the densest nearest neighbors by referencing the already
            # calculated density.
            n_dens = na.take(self.density,chunk_NNtags)
            max_loc = na.argmax(n_dens,axis=1)
            for i in xrange(finish - start + 1): # +1 for fortran counting.
                j = start + i - 1 # -1 for fortran counting.
                self.densestNN[j] = chunk_NNtags[i,max_loc[i]]
            start = finish + 1
        del chunk_NNtags, max_loc, n_dens
    
    def _build_chains(self):
        """
        Build the first round of particle chains. If the particle is too low in
        density, move on.
        """
        chainIDmax = 0
        self.densest_in_chain = {} # chainID->part ID, one to one
        for i in xrange(self.size):
            # If it's already in a group, move on, or if this particle is
            # in the padding, move on because chains can only terminate in
            # the padding, not begin, or if this particle is too low in
            # density, move on.
            if self.chainID[i] > -1 or not self.is_inside[i] or \
                    self.density[i] < self.threshold:
                continue
            chainIDnew = self._recurse_links(i, chainIDmax)
            # If the new chainID returned is the same as we entered, the chain
            # has been named chainIDmax, so we need to start a new chain
            # in the next loop.
            if chainIDnew == chainIDmax:
                chainIDmax += 1
        self.padded_particles = na.array(self.padded_particles, dtype='int64')
        return chainIDmax
    
    def _recurse_links(self, pi, chainIDmax):
        """
        Recurse up the chain to a) a self-highest density particle,
        b) a particle that already has a chainID, then turn it back around
        assigning that chainID to where we came from. If c) which
        is a particle in the padding, terminate the chain right then
        and there, because chains only go one particle deep into the padding.
        """
        nn = self.densestNN[pi]
        inside = self.is_inside[pi]
        nn_chainID = self.chainID[nn]
        # Linking to an already chainID-ed particle (don't make links from 
        # padded particles!)
        if nn_chainID > -1 and inside:
            self.chainID[pi] = nn_chainID
            return nn_chainID
        # If pi is a self-most dense particle or inside the padding, end/create
        # a new chain.
        elif nn == pi or not inside:
            self.chainID[pi] = chainIDmax
            self.densest_in_chain[chainIDmax] = self.density[pi]
            # if this is a padded particle, record it for later
            if not inside:
                self.padded_particles.append(pi)
            return chainIDmax
        # Otherwise, recursively link to nearest neighbors.
        else:
            chainIDnew = self._recurse_links(nn, chainIDmax)
            self.chainID[pi] = chainIDnew
            return chainIDnew

    def _globally_assign_chainIDs(self, chain_count):
        """
        Convert local chainIDs into globally unique chainIDs.
        """
        # First find out the number of chains on each processor.
        self.mine, chain_info = self._mpi_info_dict(chain_count)
        self.nchains = sum(chain_info.values())
        # Figure out our offset.
        self.my_first_id = sum([v for k,v in chain_info.items() if k < self.mine])
        # Change particle IDs, -1 always means no chain assignment.
        select = (self.chainID != -1)
        select = select * self.my_first_id
        self.chainID += select
        del select

    def _create_global_densest_in_chain(self):
        """
        With the globally unique chainIDs, update densest_in_chain.
        """
        # Shift the values over.
        temp = self.densest_in_chain.copy()
        self.densest_in_chain = {}
        for dens in temp:
            self.densest_in_chain[dens + self.my_first_id] = temp[dens]
        # Distribute this.
        self.densest_in_chain = self._mpi_joindict(self.densest_in_chain)
        del temp

    def _communicate_uphill_info(self):
        """
        Communicate the links to the correct neighbors from uphill_info.
        """
        # Find out how many particles we're going to receive, and make arrays
        # of the right size and type to store them.
        to_recv_count = 0
        size_max = 0
        for shift_index in xrange(27):
            # Skip ourselves, [0,0,0]
            if shift_index==13: continue
            shift = self._translate_shift(i=shift_index)
            neighbor = self._mpi_find_neighbor_3d(shift)
            # Skip ourselves
            if neighbor == self.mine: continue
            opp_shift = -1 * shift
            opp_neighbor = self._mpi_find_neighbor_3d(opp_shift)
            to_recv_count += self.global_padded_count[opp_neighbor]
            if size_max < self.global_padded_count[opp_neighbor]:
                size_max = self.global_padded_count[opp_neighbor]
        self.recv_real_indices = na.empty(to_recv_count, dtype='int64')
        self.recv_chainIDs = na.empty(to_recv_count, dtype='int64')
        # It appears that MPI doesn't like it when I slice the buffer, so I 
        # still need temporary arrays, but they don't need to be so large, but
        # still as large as the largest single thing I might receive. This size
        # was calculated above.
        temp_indices = na.empty(size_max, dtype='int64')
        temp_chainIDs = na.empty(size_max, dtype='int64')
        # Send/receive padded particles to our 26 neighbors.
        so_far = 0
        for shift_index in xrange(27):
            # Skip ourselves, [0,0,0]
            if shift_index==13: continue
            shift = self._translate_shift(i=shift_index)
            neighbor = self._mpi_find_neighbor_3d(shift)
            # Skip ourselves
            if neighbor == self.mine: continue
            opp_shift = -1 * shift
            opp_neighbor = self._mpi_find_neighbor_3d(opp_shift)
            opp_size = self.global_padded_count[opp_neighbor]
            MPI.COMM_WORLD.Sendrecv(
                [self.uphill_real_indices, MPI.INT], neighbor, 0,
                [temp_indices, MPI.INT], opp_neighbor, 0)
            MPI.COMM_WORLD.Sendrecv(
                [self.uphill_chainIDs, MPI.INT], neighbor, 0,
                [temp_chainIDs, MPI.INT], opp_neighbor, 0)
            # Only save the part of the buffer that we want to the right places
            # in the full listing.
            self.recv_real_indices[so_far:(so_far + opp_size)] = \
                temp_indices[0:opp_size]
            self.recv_chainIDs[so_far:(so_far + opp_size)] = \
                temp_chainIDs[0:opp_size]
            so_far += opp_size
        del temp_indices, temp_chainIDs

    def _recurse_global_chain_links(self, chainID_translate_map_global, chainID):
        """
        Step up the global chain links until we reach the self-densest chain,
        very similarly to the recursion of particles to densest nearest
        neighbors.
        """
        new_chainID = chainID_translate_map_global[chainID]
        if  new_chainID== chainID:
            return int(chainID)
        else:
            return self._recurse_global_chain_links(chainID_translate_map_global, new_chainID)
        

    def _connect_chains_across_tasks(self):
        """
        Using the uphill links of chains, chains are linked across boundaries.
        Chains that link to a remote chain are recorded, and a complete dict
        of chain connections is created, globally. Then chainIDs are
        reassigned recursively, assigning the ID of the most dense chainID
        to every chain that links to it.
        """
        # Remote (lower dens) chain -> local (higher) chain.
        chainID_translate_map_local = {}
        # Build the stuff to send.
        self.uphill_real_indices = self.index[self.padded_particles]
        self.uphill_chainIDs = self.chainID[self.padded_particles]
        # Now we make a global dict of how many particles each task is
        # sending.
        self.global_padded_count = {self.mine:self.uphill_chainIDs.size}
        self.global_padded_count = self._mpi_joindict(self.global_padded_count)
        # Send/receive 'em.
        self._communicate_uphill_info()
        # Fix the IDs to localIDs.
        for i,real_index in enumerate(self.recv_real_indices):
            try:
                localID = self.rev_index[real_index]
                # We don't want to update the chainIDs of my padded particles.
                # Remember we are supposed to be only considering particles
                # in my *real* region, that are padded in my neighbor.
                if not self.is_inside[localID]:
                    # Make it negative so we can skip it below.
                    self.recv_real_indices[i] = -1
                    continue
                self.recv_real_indices[i] = localID
            except KeyError:
                # This is probably a particle we don't even own, so we want
                # to ignore it.
                self.recv_real_indices[i] = -1
                continue
        # Now relate the local chainIDs to the received chainIDs
        for i,localID in enumerate(self.recv_real_indices):
            # If the 'new' chainID is different that what we already have,
            # we need to record it, but we skip particles that were assigned
            # -1 above.
            if localID != -1:
                if self.recv_chainIDs[i] != self.chainID[localID]:
                    chainID_translate_map_local[self.recv_chainIDs[i]] = \
                        self.chainID[localID]
        # In chainID_translate_map_local, chains may
        # 'point' to only one chain, but a chain may have many that point to
        # it. Therefore each key (a chain) in this dict is unique, but the items
        # the keys point to are not necessarily unique. Most chains do not
        # have an entry in the dict, and that is addressed after the 
        # communication step.
        chainID_translate_map_global = \
            self._mpi_joindict(chainID_translate_map_local)
        for i in xrange(self.nchains):
            try:
                target = chainID_translate_map_global[i]
            except KeyError:
                # The the chain is either self most-dense, or a singleton.
                chainID_translate_map_global[i] = int(i)
        # Build a list of chain densities, sorted smallest to largest.
        dens_temp = []
        for key in self.densest_in_chain:
            item = [self.densest_in_chain[key], key]
            insort(dens_temp, item)
        # Loop over chains, smallest to largest density, recursively until
        # we reach a self-assigned chain. Then we assign that final chainID to
        # the *current* one only.
        for chain in dens_temp:
            new_chainID = \
                self._recurse_global_chain_links(chainID_translate_map_global, chain[1])
            chainID_translate_map_global[chain[1]] = new_chainID
            # At the same time, remove chains from densest_in_chain that have
            # been reassigned.
            if chain[1] != new_chainID:
                del self.densest_in_chain[chain[1]]
                # Also fix nchains to keep up.
                self.nchains -= 1
        # Convert local particles to their new chainID
        for i in xrange(self.size):
            old_chainID = self.chainID[i]
            if old_chainID == -1: continue
            new_chainID = chainID_translate_map_global[old_chainID]
            self.chainID[i] = new_chainID
        del chainID_translate_map_local, dens_temp, self.recv_chainIDs
        del self.recv_real_indices, self.uphill_real_indices, self.uphill_chainIDs

    def _communicate_annulus_chainIDs(self):
        """
        Transmit all of our chainID-ed particles that are within self.padding
        of the boundaries to all of our neighbors. Tests show that this is
        faster than trying to figure out which of the neighbors to send the data
        to.
        """
        # Pick the particles in the annulus.
        real_indices = self.index[self.is_inside_annulus]
        chainIDs = self.chainID[self.is_inside_annulus]
        # Eliminate un-assigned particles.
        select = (chainIDs != -1)
        real_indices = real_indices[select]
        chainIDs = chainIDs[select]
        send_count = real_indices.size
        # Here distribute the counts globally. Unfortunately, it's a barrier(), 
        # but there's so many places in this that need to be globally synched
        # that it's not worth the effort right now to make this one spot better.
        global_annulus_count = {self.mine:send_count}
        global_annulus_count = self._mpi_joindict(global_annulus_count)

        size_max = max(global_annulus_count.values())
        recv_real_indices = na.empty(size_max, dtype='int64')
        recv_chainIDs = na.empty(size_max, dtype='int64')

        for shift_index in xrange(27):
            if shift_index==13: continue
            shift = self._translate_shift(i=shift_index)
            neighbor = self._mpi_find_neighbor_3d(shift)
            if neighbor == self.mine: continue
            opp_shift = -1 * shift
            opp_neighbor = self._mpi_find_neighbor_3d(opp_shift)
            opp_size = global_annulus_count[opp_neighbor]
            # Send/receive the data.
            status = MPI.COMM_WORLD.Sendrecv(
                [real_indices, MPI.INT], neighbor, 0,
                [recv_real_indices, MPI.INT], opp_neighbor, 0)
            status = MPI.COMM_WORLD.Sendrecv(
                [chainIDs, MPI.INT], neighbor, 0,
                [recv_chainIDs, MPI.INT], opp_neighbor, 0)
            # Use the data immediately.
            for i,real_index in enumerate(recv_real_indices[0:opp_size]):
                try:
                    localID = self.rev_index[real_index]
                    # We are only updating our particles that are in our
                    # padding, so to be rigorous we will skip particles
                    # that are in our real region.
                    if self.is_inside[localID]:
                        continue
                    self.chainID[localID] = recv_chainIDs[i]
                except KeyError:
                    # We ignore data that's not for us.
                    continue
        del recv_real_indices, recv_chainIDs, real_indices, chainIDs, select

    def _connect_chains(self):
        """
        With the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        self.chain_densest_n = {} # chainID -> {chainIDs->boundary dens}
        self.reverse_map = {} # chainID -> groupID, one to one
        for chainID in self.densest_in_chain:
            self.chain_densest_n[int(chainID)] = {}
            self.reverse_map[int(chainID)] = -1
        for i in xrange(self.size):
            # Don't consider this particle if it's not part of a chain.
            if self.chainID[i] < 0: continue
            # If this particle is in the padding, don't make a connection.
            if not self.is_inside[i]: continue
            # Find this particle's chain max_dens.
            part_max_dens = self.densest_in_chain[self.chainID[i]]
            # Loop over nMerge closest nearest neighbors.
            for j in xrange(self.nMerge+2):
                thisNN = self.NNtags[i,j]
                thisNN_chainID = self.chainID[thisNN]
                # If our neighbor is in the same chain, move on.
                if self.chainID[i] == thisNN_chainID: continue
                # No introspection, nor our connected NN.
                # This is probably the same as above, but it's OK.
                # It can be removed later.
                if thisNN==i or thisNN==self.densestNN[i]: continue
                # Everything immediately below is for
                # neighboring particles with a chainID. 
                if thisNN_chainID >= 0:
                    # Find thisNN's chain's max_dens.
                    thisNN_max_dens = self.densest_in_chain[thisNN_chainID]
                    # Calculate the two groups boundary density.
                    boundary_density = (self.density[thisNN] + self.density[i]) / 2.
                    # Find out who's denser.
                    if thisNN_max_dens >= part_max_dens:
                        higher_chain = thisNN_chainID
                        lower_chain = self.chainID[i]
                    else:
                        higher_chain = self.chainID[i]
                        lower_chain = thisNN_chainID
                    # See if this boundary density is higher than
                    # previously recorded for this pair of chains.
                    # Links only go 'downhill'.
                    try:
                        old = self.chain_densest_n[higher_chain][lower_chain]
                        if old < boundary_density:
                            # make this the new densest boundary between this pair
                            self.chain_densest_n[higher_chain][lower_chain] = \
                                boundary_density
                    except KeyError:
                        # we haven't seen this pairing before, record this as the
                        # new densest boundary between chains
                        self.chain_densest_n[higher_chain][lower_chain] = \
                            boundary_density
                else:
                    continue

    def _make_global_chain_densest_n(self):
        """
        We want to record the maximum boundary density between all chains on
        all tasks.
        """
        # Remember that every task has every active chainID as a key in 
        # chain_densest_n (see the beginning of _connect_chains())
        # Uses barrier(), but that's needed.
        for higher_chain in self.chain_densest_n:
            self.chain_densest_n[higher_chain] = \
                self._mpi_maxdict(self.chain_densest_n[higher_chain])
    
    def _build_groups(self):
        """
        With the collection of possible chain links, build groups.
        """
        g_high = []
        g_low = []
        g_dens = []
        densestbound = {} # chainID -> boundary density
        for chainID in self.densest_in_chain:
            densestbound[chainID] = -1.0
        groupID = 0
        # First assign a group to all chains with max_dens above peakthresh.
        for chainID in self.densest_in_chain:
            if self.densest_in_chain[chainID] >= self.peakthresh:
                self.reverse_map[chainID] = groupID
                groupID += 1
        # Loop over all of the chain linkages.
        for chain_high in self.chain_densest_n:
            for chain_low in self.chain_densest_n[chain_high]:
                dens = self.chain_densest_n[chain_high][chain_low]
                max_dens_high = self.densest_in_chain[chain_high]
                max_dens_low = self.densest_in_chain[chain_low]
                # If neither are peak density groups, mark them for later
                # consideration.
                if max_dens_high < self.peakthresh and \
                    max_dens_low < self.peakthresh:
                        g_high.append(chain_high)
                        g_low.append(chain_low)
                        g_dens.append(dens)
                        continue
                # If both are peak density groups, and have a boundary density
                # that is high enough, make them into a group, otherwise
                # move onto another linkage.
                if max_dens_high >= self.peakthresh and \
                        max_dens_low >= self.peakthresh:
                    if dens < self.saddlethresh:
                        continue
                    else:
                        group_high = self.reverse_map[chain_high]
                        group_low = self.reverse_map[chain_low]
                        # Both are already identified as groups, so we need
                        # to re-assign the less dense group to the denser
                        # groupID.
                        for chID in self.reverse_map:
                            if self.reverse_map[chID] == group_low:
                                self.reverse_map[chID] = group_high
                        continue
                # Else, one is above peakthresh, the other below
                # find out if this is the densest boundary seen so far for
                # the lower chain.
                group_high = self.reverse_map[chain_high]
                if dens > densestbound[chain_low]:
                    densestbound[chain_low] = dens
                    self.reverse_map[chain_low] = group_high
                # Done double loop over links.
        """
        Now the fringe chains are connected to the proper group
        (>peakthresh) with the largest boundary.  But we want to look
        through the boundaries between fringe groups to propagate this
        along.  Connections are only as good as their smallest boundary
        """
        changes = 1
        while changes:
            changes = 0
            for j,dens in enumerate(g_dens):
                chain_high = g_high[j]
                chain_low = g_low[j]
                # If the density of this boundary and the densestbound of
                # the other chain is higher than a chain's densestbound, then
                # replace it.
                if dens > densestbound[chain_low] and \
                        densestbound[chain_high] > densestbound[chain_low]:
                    changes += 1
                    if dens < densestbound[chain_high]:
                        densestbound[chain_low] = dens
                    else:
                        densestbound[chain_low] = densestbound[chain_high]
                    self.reverse_map[chain_low] = self.reverse_map[chain_high]
        # Now we have to find the unique groupIDs, since they may have been
        # merged.
        temp = []
        for chain in self.reverse_map:
            temp.append(self.reverse_map[chain])
        # Uniquify the list.
        temp = list(Set(temp))
        # Remove -1 from the list.
        temp.pop(temp.index(-1))
        # Make a secondary map to make the IDs consecutive.
        secondary_map = {}
        for i,ID in enumerate(temp):
            secondary_map[ID] = i
        # Update reverse_map
        for chain in self.reverse_map:
            # Don't attempt to fix non-assigned chains.
            if self.reverse_map[chain] == -1: continue
            self.reverse_map[chain] = secondary_map[self.reverse_map[chain]]
        group_count = len(temp)
        del secondary_map, temp, g_high, g_low, g_dens, densestbound
        return group_count

    def _translate_groupIDs(self, group_count):
        """
        Using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        for i in xrange(self.size):
            # Don't translate non-affiliated particles.
            if self.chainID[i] == -1: continue
            # We want to remove the group tag from padded particles,
            # so when we return it to HaloFinding, there is no duplication.
            if self.is_inside[i]:
                self.chainID[i] = self.reverse_map[self.chainID[i]]
            else:
                self.chainID[i] = -1
        # Create a densest_in_group, analogous to densest_in_chain.
        self.densest_in_group = {}
        for i in xrange(group_count):
            self.densest_in_group[i] = 0.0
        for chainID in self.densest_in_chain:
            groupID = self.reverse_map[chainID]
            if groupID == -1: continue
            max_dens = self.densest_in_chain[chainID]
            if self.densest_in_group[groupID] < max_dens:
                self.densest_in_group[groupID] = max_dens

    def _chain_hop(self):
        mylog.info('Building kd tree for %d particles...' % \
            self.size)
        self._init_kd_tree()
        # Mark particles in as being in/out of the domain.
        self._is_inside()
        # Loop over the particles to find NN for each.
        mylog.info('Finding nearest neighbors/density...')
        chainHOP_tags_dens()
        self.density = fKD.dens
        self.NNtags = (fKD.nn_tags - 1).transpose()
        # We can free these right now, the rest later.
        del fKD.dens, fKD.nn_tags, fKD.mass, fKD.dens
        del self.mass, self.xpos, self.ypos, self.zpos
        count = len(na.where(self.density >= self.threshold)[0])
        print 'count above thresh', count
        # Now each particle has NNtags, and a local self density.
        # Let's find densest NN
        mylog.info('Finding densest nearest neighbors...')
        self._densestNN()
        # Now we can free these.
        del fKD.pos, fKD.chunk_tags
        free_tree() # Frees the kdtree object.
        # Build the chain of links.
        mylog.info('Building particle chains...')
        chain_count = self._build_chains()
        mylog.info('Gobally assigning chainIDs...')
        self._globally_assign_chainIDs(chain_count)
        mylog.info('Globally finding densest in chains...')
        self._create_global_densest_in_chain()
        mylog.info('Building chain connections across tasks...')
        self._connect_chains_across_tasks()
        mylog.info('Communicating connected chains...')
        self._communicate_annulus_chainIDs()
        mylog.info('Connecting %d chains into groups...' % self.nchains)
        self._connect_chains()
        mylog.info('Communicating group links globally...')
        self._make_global_chain_densest_n()
        mylog.info('Building final groups...')
        group_count = self._build_groups()
        self.group_count = group_count
        mylog.info('Remapping particles to final groups...')
        self._translate_groupIDs(group_count)
        mylog.info('Found %d groups...' % group_count)

