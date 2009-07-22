import math,sys #, cPickle
from bisect import insort
from sets import Set

from yt.lagos import *
from yt.extensions.kdtree import *
from yt.lagos.ParallelTools import _send_array, _recv_array
from mpi4py import MPI
from Forthon import *

class PaddedPart(object):
    def __init__(self, local_index, particle_index, chainID, shift):
        self.local_index = local_index
        self.particle_index = particle_index
        self.chainID = chainID
        self.shift = shift

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
        self.index = index
        self.rev_index = {}
        for i in xrange(self.size):
            self.rev_index[self.index[i]] = i
        self.mass = mass
        self.sends = []
        self.recvs = {}
        self.padded_particles = []
        self.nMerge = 4
        self.chainID = na.ones(self.size,dtype='l') * -1
        self._chain_hop()
        self.padded_particles = na.array(self.padded_particles)
    

    def _init_kd_tree(self):
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
        # This actually copies the data into the fortran space
        fKD.pos[0, :] = self.xpos
        fKD.pos[1, :] = self.ypos
        fKD.pos[2, :] = self.zpos
        fKD.mass = self.mass
        # Now call the fortran
        create_tree()

    def _is_inside(self):
        # Test to see if the points are in the 'real' region
        (LE, RE) = self.bounds
        points = fKD.pos.transpose()
        self.is_inside = ( (points >= LE).all(axis=1) * \
            (points < RE).all(axis=1) )
        # Test to see if the points are one padding distance inside the
        # regular boundaries.
        LE += self.padding
        RE -= self.padding
        inner = na.invert( (points >= LE).all(axis=1) * \
            (points < RE).all(axis=1) )
        # After inverting the logic above, we want points that are both
        # inside the padding, but within one padding of the boundary, and this
        # will do it.
        self.is_inside_real = na.bitwise_and(self.is_inside, inner)

    def _find_shift_padded(self, i):
        # Find the shift vector for a padded particle
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
        d = {0:[-1,-1,-1], 1:[-1,-1,0], 2:[-1,-1,1], 3:[-1,0,-1], 4:[-1,0,0],
            5:[-1,0,1], 6:[-1,1,-1], 7:[-1,1,0], 8:[-1,1,1,], 9:[0,-1,-1],
            10:[0,-1,0], 11:[0,-1,1], 12:[0,0,-1], 13:[0,0,0], 14:[0,0,1],
            15:[0,1,-1], 16:[0,1,0], 17:[0,1,1], 18:[1,-1,-1], 19:[1,-1,0],
            20:[1,-1,1], 21:[1,0,-1], 22:[1,0,0], 23:[1,0,1], 24:[1,1,-1],
            25:[1,1,0], 26:[1,1,1]}
        for key in d:
            d[key] = na.array(d[key], dtype='int64')
        if i==None and shift==None: return None
        if i!=None:
            return d[i]
        if shift!=None:
            shift = na.array(shift)
            for key in d:
                if (d[key] == shift).all():
                    return key

    def _densestNN(self):
        """
        The first search of nearest neighbors did not return all 
        num_neighbor neighbors, so we need to do it again, but we're not
        keeping the all of this data, just using it.
        """
        self.densestNN = na.ones(self.size,dtype='l')
        # We find nearest neighbors in chunks
        chunksize = 10000
        
        start = 1 # fortran counting!
        finish = 0
        while finish < self.size:
            finish = min(finish+chunksize,self.size)
            fKD.chunk_tags = empty((self.num_neighbors,finish-start+1), dtype='l')
            # call the fortran
            fKD.start = start
            fKD.finish = finish
            find_chunk_nearest_neighbors()
            chunk_NNtags = (fKD.chunk_tags - 1).transpose()
            # find the densest nearest neighbors
            n_dens = na.take(self.density,chunk_NNtags)
            max_loc = na.argmax(n_dens,axis=1)
            for i in xrange(finish - start + 1): # +1 for fortran counting
                j = start + i - 1 # -1 for fortran counting
                self.densestNN[j] = chunk_NNtags[i,max_loc[i]]
            start = finish + 1
    
    def _build_chains(self):
        """
        Build the first round of particle chains. If the particle is too low in
        density, move on.
        """
        chainIDmax = 0
        self.densest_in_chain = {} # chainID->part ID, one to one
        for i in xrange(self.size):
            # if it's already in a group, move on, or if this particle is
            # in the padding, move on because chains can only terminate in
            # the padding, not begin, or if this particle is too low in
            # density, move on
            if self.chainID[i] > -1 or not self.is_inside[i] or \
                self.density[i] < self.threshold: continue
            chainIDnew = self._recurse_links(i, chainIDmax)
            # if the new chainID returned is the same as we entered, the chain
            # has been named chainIDmax, so we need to start a new chain
            # in the next loop
            if chainIDnew == chainIDmax:
                chainIDmax += 1
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
        # linking to an already chainID-ed particle (don't make links from 
        # padded particles!)
        if nn_chainID > -1 and inside:
            self.chainID[pi] = nn_chainID
            return nn_chainID
        # if pi is a self-most dense particle or inside the padding, end/create
        # a new chain
        elif nn == pi or not inside:
            self.chainID[pi] = chainIDmax
            self.densest_in_chain[chainIDmax] = self.density[pi]
            # if this is a padded particle, record it for later
            if not inside:
                self.padded_particles.append(pi)
            return chainIDmax
        # otherwise, recursively link to nearest neighbors
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
        # Change particle IDs, -1 still means no chain.
        for i in xrange(self.size):
            if self.chainID[i] != -1:
                self.chainID[i] += self.my_first_id

    def _create_global_densest_in_chain(self):
        """
        With the globally unique chainIDs, update densest_in_chain.
        """
        # shift the values over
        temp = self.densest_in_chain.copy()
        self.densest_in_chain = {}
        for dens in temp:
            self.densest_in_chain[dens + self.my_first_id] = temp[dens]
        # distribute this
        self.densest_in_chain = self._mpi_joindict(self.densest_in_chain)

    def _build_uphill_info(self):
        """
        Use self.padded_particles to build the data that will be transferred
        to neighboring regions.
        """
        self.uphill_real_indices = []
        self.uphill_chainIDs = []
        for part in self.padded_particles:
            # figure out the real particle_index
            real_index = self.index[part]
            chainID = self.chainID[part]
            self.uphill_real_indices.append(real_index)
            self.uphill_chainIDs.append(chainID)
        self.uphill_real_indices = na.array(self.uphill_real_indices, dtype='int64')
        self.uphill_chainIDs = na.array(self.uphill_chainIDs, dtype='int64')

    def _communicate_uphill_info(self):
        """
        Communicate the links to the correct neighbors from uphill_info.
        """
        self.recv_real_indices = []
        self.recv_chainIDs = []
        size_max = max(self.global_padded_count.values())
        temp_indices = na.empty(size_max, dtype='int64')
        temp_chainIDs = na.empty(size_max, dtype='int64')
        # Send/receive padded particles to our 26 neighbors
        for shift_index in xrange(27):
            # Skip ourselves, [0,0,0]
            if shift_index==13: continue
            shift = self._translate_shift(i=shift_index)
            neighbor = self._find_neighbor_3d(shift)
            # Skip ourselves
            if neighbor == self.mine: continue
            opp_shift = -1 * shift
            opp_neighbor = self._find_neighbor_3d(opp_shift)
            opp_size = self.global_padded_count[opp_neighbor]
            status = MPI.COMM_WORLD.Sendrecv(
                [self.uphill_real_indices, MPI.LONG_INT], neighbor, 0,
                [temp_indices, MPI.LONG_INT], opp_neighbor, 0)
            status = MPI.COMM_WORLD.Sendrecv(
                [self.uphill_chainIDs, MPI.LONG_INT], neighbor, 0,
                [temp_chainIDs, MPI.LONG_INT], opp_neighbor, 0)
            self.recv_real_indices.extend(temp_indices[0:opp_size].tolist())
            self.recv_chainIDs.extend(temp_chainIDs[0:opp_size].tolist())
        self.recv_real_indices = na.array(self.recv_real_indices, dtype='int64')
        self.recv_chainIDs = na.array(self.recv_chainIDs, dtype='int64')
        #print MPI.COMM_WORLD.rank,'max',max(self.recv_real_indices),'min',min(self.recv_real_indices)

    def _recurse_global_chain_links(self, chainID_translate_map_global, chainID):
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
        # remote (lower dens) chain -> local (higher) chain
        chainID_translate_map_local = {}
        # build the stuff to send
        self._build_uphill_info()
        # now we make a global dict of how many particles each task is
        # sending
        self.global_padded_count = {self.mine:len(self.uphill_chainIDs)}
        self.global_padded_count = self._mpi_joindict(self.global_padded_count)
        # send/receive 'em
        self._communicate_uphill_info()
        # fix the IDs to localIDs
        for i,real_index in enumerate(self.recv_real_indices):
            try:
                localID = self.rev_index[real_index]
                self.recv_real_indices[i] = localID
            except KeyError:
                pass
        # Now relate the local chainIDs to the received chainIDs
        for i,localID in enumerate(self.recv_real_indices):
            # if the 'new' chainID is different that what we already have,
            # we need to record it
            if self.recv_chainIDs[i] != self.chainID[localID]:
                chainID_translate_map_local[self.recv_chainIDs[i]] = \
                    self.chainID[localID]
        # chainID_translate_map_local is now a 'sparse' dict. Chains may
        # 'point' to only one chain, but a chain may have many that point to
        # it. Therefore each key in this dict is unique, but the items
        # the keys point to are not necessarily unique. Most chains do not
        # have an entry in the sparse dict, and that is addressed.
        chainID_translate_map_global = \
            self._mpi_joindict(chainID_translate_map_local)
        for i in xrange(self.nchains):
            try:
                target = chainID_translate_map_global[i]
            except KeyError:
                chainID_translate_map_global[i] = int(i)
        # Build a list of chain densities, sorted smallest to largest
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
        # I may eventually do this more efficently, but for now...
        # Convert local particles to their new chainID
        for i in xrange(self.size):
            old_chainID = self.chainID[i]
            if old_chainID == -1: continue
            new_chainID = chainID_translate_map_global[old_chainID]
            self.chainID[i] = new_chainID

    def _communicate_annulus_chainIDs(self):
        """
        Transmit all of our chainID-ed particles that are within self.padding
        of the boundaries to all of our neighbors. Tests show that this is
        faster than trying to figure out which of the neighbors to send the data
        to.
        """
        real_indices = []
        chainIDs = []
        for i in xrange(self.size):
            # We only want the particles in the 'annulus'
            if not self.is_inside_real[i]: continue
            chainID = self.chainID[i]
            # non-assigned particles are not interesting
            if chainID == -1: continue
            real_index = self.index[i]
            real_indices.append(real_index)
            chainIDs.append(chainID)
        real_indices = na.array(real_indices, dtype='int64')
        chainIDs = na.array(chainIDs, dtype='int64')
        send_count = len(real_indices)
        global_annulus_count = {self.mine:send_count}
        global_annulus_count = self._mpi_joindict(global_annulus_count)
        #print global_annulus_count

        size_max = max(global_annulus_count.values())
        recv_real_indices = na.empty(size_max, dtype='int64')
        recv_chainIDs = na.empty(size_max, dtype='int64')

        for shift_index in xrange(27):
            if shift_index==13: continue
            # send the data
            shift = self._translate_shift(i=shift_index)
            neighbor = self._find_neighbor_3d(shift)
            # there's no need to send to myself...
            if neighbor == self.mine: continue
            opp_shift = -1 * shift
            opp_neighbor = self._find_neighbor_3d(opp_shift)
            opp_size = global_annulus_count[opp_neighbor]
            status = MPI.COMM_WORLD.Sendrecv(
                [real_indices, MPI.LONG_INT], neighbor, 0,
                [recv_real_indices, MPI.LONG_INT], opp_neighbor, 0)
            status = MPI.COMM_WORLD.Sendrecv(
                [chainIDs, MPI.LONG_INT], neighbor, 0,
                [recv_chainIDs, MPI.LONG_INT], opp_neighbor, 0)
            # use the data
            for i,real_index in enumerate(recv_real_indices[0:opp_size]):
                try:
                    localID = self.rev_index[real_index]
                    self.chainID[localID] = recv_chainIDs[i]
                except KeyError:
                    continue
        del recv_real_indices, recv_chainIDs
        sys.exit()

    def _connect_chains(self,chain_count):
        """
        With the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        self.chain_densest_n = {} # chainID -> {chainIDs->boundary dens}
        for i in xrange(chain_count):
            self.chain_densest_n[i] = {}
        self.reverse_map = {} # chainID -> groupID, one to one
        for i in xrange(chain_count):
            self.reverse_map[i] = -1
        # This prevents key problems, and is deleted below.
        self.reverse_map[-1] = -1
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
                # particles with a chainID. 
                if thisNN_chainID >= 0:
                    # If our neighbor is not part of a chain, move on
                    #if thisNN_chainID < 0: continue
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
        del self.reverse_map[-1]
    
    def _build_groups(self, chain_count):
        """
        With the collection of possible chain links, build groups.
        """
        g_high = []
        g_low = []
        g_dens = []
        densestbound = {} # chainID -> boundary density
        for i in xrange(chain_count):
            densestbound[i] = -1.0
        groupID = 0
        # first assign a group to all chains with max_dens above peakthresh
        for chainID in self.densest_in_chain:
            if self.densest_in_chain[chainID] >= self.peakthresh:
                self.reverse_map[chainID] = groupID
                groupID += 1
        # loop over all of the chain linkages
        for chain_high in self.chain_densest_n:
            for chain_low in self.chain_densest_n[chain_high]:
                dens = self.chain_densest_n[chain_high][chain_low]
                max_dens_high = self.densest_in_chain[chain_high]
                max_dens_low = self.densest_in_chain[chain_low]
                # if neither are peak density groups, mark them for later
                # consideration
                if max_dens_high < self.peakthresh and \
                    max_dens_low < self.peakthresh:
                        g_high.append(chain_high)
                        g_low.append(chain_low)
                        g_dens.append(dens)
                        continue
                # if both are peak density groups, and have a boundary density
                # that is high enough, make them into a group, otherwise
                # move onto another linkage
                if max_dens_high >= self.peakthresh and \
                        max_dens_low >= self.peakthresh:
                    if dens < self.saddlethresh:
                        continue
                    else:
                        group_high = self.reverse_map[chain_high]
                        group_low = self.reverse_map[chain_low]
                        # both are already identified as groups, so we need
                        # to re-assign the less dense group to the denser
                        # groupID
                        for chID in self.reverse_map:
                            if self.reverse_map[chID] == group_low:
                                self.reverse_map[chID] = group_high
                        continue
                # else, one is above peakthresh, the other below
                # find out if this is the densest boundary seen so far for
                # the lower chain
                group_high = self.reverse_map[chain_high]
                if dens > densestbound[chain_low]:
                    densestbound[chain_low] = dens
                    self.reverse_map[chain_low] = group_high
                # done double loop over links
        
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
        temp = list(Set(temp)) # uniquify temp
        temp.sort()
        # Make the groupIDs run 0 to group_count-1
        for t,gID in enumerate(temp[1:]): # -1 is always the first item, skip it
            if gID != t:
                for chain in self.reverse_map:
                    if self.reverse_map[chain] == gID:
                        self.reverse_map[chain] = t
                
        group_count = len(temp) - 1 # don't count the -1 groupID
        return group_count

    
    def _pare_groups_by_max_dens(self):
        """
        For all groups, check to see that the densest particle in the group is
        dense enough.
        ***** Warning! **********
        This is more or less wrong, and doesn't belong here anyway!
        """
        for groupID in self.chain_connections:
            max_dens = max( (self.densest_in_Chain[chainID] for chainID in 
                            self.chain_connections[groupID] ) )
            # Too low, elminate the group by emptying it.
            if max_dens < 3 * self.threshold:
                # at the same time elminate the associated chains
                for chainID in self.chain_connections[groupID]:
                    self.reverse_map[chainID] = -1
                self.chain_connections[groupID] = []
                print 'elminated a group'
    
    def _translate_groupIDs(self, group_count):
        """
        Using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        for i in xrange(self.size):
            # Don't translate non-affiliated particles.
            if self.chainID[i] == -1: continue
            self.chainID[i] = self.reverse_map[self.chainID[i]]
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
        count = len(na.where(self.density >= self.threshold)[0])
        print 'count above thresh', count
        # Now each particle has NNtags, and a local self density.
        # Let's find densest NN
        mylog.info('Finding densest nearest neighbors...')
        self._densestNN()
        # Now we can free these
        del fKD.pos, fKD.chunk_tags
        free_tree() # Frees the kdtree object.
        chain_count = 0
        for i in xrange(self.size):
            if i == self.densestNN[i]: chain_count += 1
        mylog.info('there are %d self-densest particles' % chain_count)
        # Build the chain of links.
        mylog.info('Building particle chains...')
        chain_count = self._build_chains()
        mylog.info('Gobally assigning chainIDs...')
        self._globally_assign_chainIDs(chain_count)
        mylog.info('Global densest in chain...')
        self._create_global_densest_in_chain()
        #self._build_uphill_info()
        mylog.info('Connecting chains across tasks...')
        self._connect_chains_across_tasks()
        mylog.info('Communicating padded particle chainIDs...')
        self._communicate_annulus_chainIDs()
        #mylog.info('Connecting %d chains into groups...' % chain_count)
        #self._connect_chains(chain_count)
        #group_count = self._build_groups(chain_count)
        #mylog.info('Paring %d groups...' % group_count)
        # Don't pare groups here, it has to happen later.
        #self._pare_groups_by_max_dens()
        #self._translate_groupIDs(group_count)
        #mylog.info('Found %d groups...' % group_count)


