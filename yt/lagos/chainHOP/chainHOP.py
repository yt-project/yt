import math,sys #, cPickle
from bisect import insort

from yt.extensions.kdtree import *
from Forthon import *

class RunChainHOP(object):
    def __init__(self,period, padding, num_neighbors, bounds,
            xpos, ypos, zpos, mass, threshold=160.0):
        self.threshold = threshold
        self.period = period
        self.padding = padding
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.mass = mass
        self.padded_particles = []
        self.size = len(self.xpos)
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

    def _connect_chains(self,chain_count):
        """
        With the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        self.chain_densest_n = {} # chainID -> [chainID, boundary dens]
        for i in xrange(chain_count):
            self.chain_densest_n[i] = [-1, 0.0]
        
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
                # If our neighbor is not part of a chain, move on.
                if thisNN_chainID < 0: continue
                # Find thisNN's chain's max_dens.
                thisNN_max_dens = self.densest_in_chain[thisNN_chainID]
                # Calculate the two groups boundary density.
                boundary_density = (self.density[thisNN] + self.density[i]) / 2.
                # If both these chains have higher than peakdens max_dens, they
                # can only be linked if the boundary_density is high enough.
                if thisNN_max_dens >= (3 * self.threshold) and \
                    part_max_dens >= (3 * self.threshold) and \
                    boundary_density < (2.5 * self.threshold):
                        # In this case, no link is made.
                        continue
                # Find out who's denser.
                max_dens = max(thisNN_max_dens, part_max_dens)
                if max_dens == thisNN_max_dens:
                    higher_chain = thisNN_chainID
                    lower_chain = self.chainID[i]
                else:
                    higher_chain = self.chainID[i]
                    lower_chain = thisNN_chainID
                # See if this boundary density is higher than
                # previously recorded.
                # Links only go 'uphill'.
                old = self.chain_densest_n[lower_chain]
                if old[1] < boundary_density:
                    self.chain_densest_n[lower_chain] = \
                        [higher_chain,boundary_density]
    
    
    def _build_groups(self, chain_count):
        """
        Using the links between chains by boundary density, recurse uphill
        to the most dense chain. This is very similar to _build_chains.
        """
        groupIDmax = 0
        self.reverse_map = {} # chainID -> groupID, one to one
        for i in xrange(chain_count):
            self.reverse_map[i] = -1
        # This prevents key problems, and is deleted below.
        self.reverse_map[-1] = -1
        for chainID in self.chain_densest_n:
            # If this chain is already in a group, skip it.
            if self.reverse_map[chainID] != -1: continue
            groupIDnew = self._recurse_chain_links(chainID, groupIDmax)
            # If the new groupID returned is the same as we entered, the group
            # has been named groupIDmax, so we need to start a new group
            # in the next loop.
            if groupIDnew == groupIDmax:
                groupIDmax += 1
        del self.reverse_map[-1]
        return groupIDmax

    def _recurse_chain_links(self, chainID, groupIDmax):
        """
        Recurse up to a a) un-linked chain, b) a chain that already has a
        groupID.
        """
        entry = self.chain_densest_n[chainID]
        n_groupID = self.reverse_map[entry[0]]
        # Linking to an already groupID-ed chain.
        if n_groupID != -1:
            self.reverse_map[chainID] = n_groupID
            #print 'already IDed assigning ',n_groupID,' to chain ',chainID
            return n_groupID
        # If this chain itself has no neighbors, end the linking, and create
        # a new group.
        elif entry[0] == -1:
            self.reverse_map[chainID] = groupIDmax
            #print 'no neigh assigning ',groupIDmax,' to chain ',chainID
            return groupIDmax
        # Otherwise recurse to other chains.
        else:
            groupIDnew = self._recurse_chain_links(entry[0], groupIDmax)
            self.reverse_map[chainID] = groupIDnew
            #print 'after recurse assigning ',groupIDnew,' to chain ',chainID
            return groupIDnew

    
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
        # Connect the chains into groups.
        mylog.info('Connecting %d chains into groups...' % chain_count)
        self._connect_chains(chain_count)
        group_count = self._build_groups(chain_count)
        #mylog.info('Paring %d groups...' % group_count)
        # Don't pare groups here, it has to happen later.
        #self._pare_groups_by_max_dens()
        self._translate_groupIDs(group_count)
        mylog.info('Found %d groups...' % group_count)


