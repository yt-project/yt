import math,sys #, cPickle
from bisect import insort

from yt.extensions.kdtree import *

from Forthon import *

class partNN(object):
    def __init__(self, NNtags, density, densestNN, chainID, is_inside, order_index):
        self.NNtags = NNtags # numpy array of the tags from the kdtree
        self.density = density # this particles local density
        self.densestNN = densestNN # this particles densest nearest neighbor
        self.chainID = chainID # this particles current chain identifier, changes
        self.is_inside = is_inside # True/False is inside the main volume
        self.order_index = order_index # order_index for this particle

class Run_chain_HOP(object):
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
        self.__initialize_NN()
        self._chainHOP()
        self.padded_particles = na.array(self.padded_particles)
    

    def __init_kd_tree(self,size):
        # Yes, we really do need to initialize this many arrays.
        # They're deleted in _chainHOP.
        fKD.nn_tags = empty((self.num_neighbors,size),dtype='l')
        fKD.nn_dist = empty((self.num_neighbors,size),dtype='d')
        fKD.dens = zeros(size,dtype='d')
        fKD.mass = empty(size,dtype='d')
        fKD.pos = empty((3,size),dtype='d')
        fKD.qv = empty(3,dtype='d')
        fKD.nn = self.num_neighbors
        fKD.nparts = size
        fKD.sort = True # slower, but needed in _connect_chains
        fKD.rearrange = True # faster, more memory
        # this actually copies the data into the fortran space
        fKD.pos[0,:] = self.xpos
        fKD.pos[1,:] = self.ypos
        fKD.pos[2,:] = self.zpos
        fKD.mass = self.mass
        # now call the fortran
        create_tree()

    def __initialize_NN(self):
        size = len(self.xpos)
        self.NN = []
        for i in xrange(size):
            NNtemp = partNN(na.empty([self.num_neighbors]),0.0,-1,-1,False,-1)
            self.NN.append(NNtemp)

    def _is_inside(self,point):
        # test to see if this point is inside the 'real' region
        (LE,RE) = self.bounds
        point = na.array(point)
        if (point >= LE).all() and (point < RE).all():
            return True
        else:
            return False
    
    def _densestNN(self, pi, NNtags):
        # find the densest nearest neighbor
        NNdens = 0.
        for tags in NNtags:
            if self.NN[tags].density > NNdens:
                NNdens = self.NN[tags].density
                NNindex = tags
        self.NN[pi].densestNN = NNindex
    
    def _build_chains(self):
        """
        build the first round of particle chains. If the particle is too low in
        density, move on.
        """
        chainIDmax = 0
        self.densest_in_chain = {} # chainID->part ID, one to one
        for part in self.NN:
            # if it's already in a group, move on, or if this particle is in the padding,
            # move on because chains can only terminate in the padding, not
            # begin, or if this particle is too low in density, move on
            if part.chainID > -1 or part.is_inside is False or part.density < self.threshold: continue
            chainIDnew = self._recurse_links(part.order_index, chainIDmax)
            # if the new chainID returned is the same as we entered, the chain
            # has been named chainIDmax, so we need to start a new chain
            # in the next loop
            if chainIDnew == chainIDmax:
                chainIDmax += 1
        return chainIDmax
    
    def _recurse_links(self, pi, chainIDmax):
        """
        recurse up the chain to a) a self-highest density particle,
        b) a particle that already has a chainID, then turn it back around
        assigning that chainID to where we came from. If c) which is a particle
        in the padding, terminate the chain right then and there, because chains
        only go one particle deep into the padding.
        """
        nn = self.NN[pi].densestNN
        inside = self.NN[pi].is_inside
        # linking to an already chainID-ed particle (don't make links from 
        # padded particles!)
        if self.NN[nn].chainID > -1 and inside is True:
            self.NN[pi].chainID = self.NN[nn].chainID
            return self.NN[nn].chainID
        # if pi is a self-most dense particle or inside the padding, end/create
        # a new chain
        elif nn == pi or inside is False:
            self.NN[pi].chainID = chainIDmax
            self.densest_in_chain[chainIDmax] = self.NN[pi].density
            # if this is a padded particle, record it
            if inside is False:
                self.padded_particles.append(pi)
            return chainIDmax
        # recursively link to nearest neighbors
        else:
            chainIDnew = self._recurse_links(nn, chainIDmax)
            self.NN[pi].chainID = chainIDnew
            return chainIDnew

    def _connect_chains(self,chain_count):
        """
        with the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        self.chain_densest_n = {} # chainID -> [chainID, boundary dens]
        for i in range(chain_count):
            self.chain_densest_n[i] = [-1,0.0]
        
        for part in self.NN:
            # don't consider this particle if it's not part of a chain
            if part.chainID < 0: continue
            # if this particle is in the padding, don't make a connection
            if part.is_inside is False: continue
            # find this particle's chain max_dens
            part_max_dens = self.densest_in_chain[part.chainID]
            # loop over nMerge closest nearest neighbors
            for i in range(self.nMerge+2):
                thisNN = part.NNtags[i]
                thisNN_chainID = self.NN[thisNN].chainID
                # if our neighbor is in the same chain, move on
                if part.chainID == thisNN_chainID: continue
                # no introspection, nor our connected NN
                # this is probably the same as above, but it's OK. Can be removed
                # later.
                if thisNN==part.order_index or thisNN==part.densestNN: continue
                # if our neighbor is not part of a chain, move on
                if thisNN_chainID < 0: continue
                # find thisNN's chain's max_dens
                thisNN_max_dens = self.densest_in_chain[thisNN_chainID]
                # calculate the two groups boundary density
                boundary_density = (self.NN[thisNN].density + part.density) / 2.
                # if both these chains have higher than peakdens max_dens, they
                # can only be linked if the boundary_density is high enough.
                if thisNN_max_dens >= (3 * self.threshold) and \
                    part_max_dens >= (3 * self.threshold) and \
                    boundary_density < (2.5 * self.threshold):
                        # in this case, no link is made
                        continue
                # find out who's denser
                max_dens = max(thisNN_max_dens,part_max_dens)
                if max_dens == thisNN_max_dens:
                    higher_chain = thisNN_chainID
                    lower_chain = part.chainID
                else:
                    higher_chain = part.chainID
                    lower_chain = thisNN_chainID
                # see if this boundary density is higher than previously recorded
                # links only go 'uphill'
                old = self.chain_densest_n[lower_chain]
                if old[1] < boundary_density:
                    self.chain_densest_n[lower_chain] = [higher_chain,boundary_density]
    
    
    def _build_groups(self,chain_count):
        """
        Using the links between chains by boundary density, recurse uphill
        to the most dense chain. This is very similar to _build_chains.
        """
        groupIDmax = 0
        self.reverse_map = {} # chainID -> groupID, one to one
        for i in range(chain_count):
            self.reverse_map[i] = -1
        # this prevents key problems, is deleted below
        self.reverse_map[-1] = -1
        for chainID in self.chain_densest_n:
            # if this chain is already in a group, skip it
            if self.reverse_map[chainID] != -1: continue
            groupIDnew = self._recurse_chain_links(chainID, groupIDmax)
            # if the new groupID returned is the same as we entered, the group
            # has been named groupIDmax, so we need to start a new group
            # in the next loop
            if groupIDnew == groupIDmax:
                groupIDmax += 1
        del self.reverse_map[-1]
        return groupIDmax

    def _recurse_chain_links(self, chainID, groupIDmax):
        """
        recurse up to a a) un-linked chain, b) a chain that already has a
        groupID.
        """
        entry = self.chain_densest_n[chainID]
        n_groupID = self.reverse_map[entry[0]]
        # linking to an already groupID-ed chain
        if n_groupID != -1:
            self.reverse_map[chainID] = n_groupID
            #print 'already IDed assigning ',n_groupID,' to chain ',chainID
            return n_groupID
        # if this chain itself has no neighbors, end the linking, and create
        # a new group
        elif entry[0] == -1:
            self.reverse_map[chainID] = groupIDmax
            #print 'no neigh assigning ',groupIDmax,' to chain ',chainID
            return groupIDmax
        # otherwise recurse to other chains
        else:
            groupIDnew = self._recurse_chain_links(entry[0],groupIDmax)
            self.reverse_map[chainID] = groupIDnew
            #print 'after recurse assigning ',groupIDnew,' to chain ',chainID
            return groupIDnew

    
    def _pare_groups_by_max_dens(self):
        """
        for all groups, check to see that the densest particle in the group is
        dense enough.
        """
        for groupID in self.chain_connections:
            max_dens = max( (self.densest_in_Chain[chainID] for chainID 
                            in self.chain_connections[groupID] ) )
            # too low, elminate the group by emptying it
            if max_dens < 3 * self.threshold:
                # at the same time elminate the associated chains
                for chainID in self.chain_connections[groupID]:
                    self.reverse_map[chainID] = -1
                self.chain_connections[groupID] = []
                print 'elminated a group'
    
    def _translate_groupIDs(self,group_count):
        """
        using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        for nn in self.NN:
            # don't translate non-affiliated particles
            if nn.chainID == -1: continue
            nn.chainID = self.reverse_map[nn.chainID]
        # create a densest_in_group, analogous to densest_in_chain
        self.densest_in_group = {}
        for i in range(group_count):
            self.densest_in_group[i] = 0.0
        for chainID in self.densest_in_chain:
            groupID = self.reverse_map[chainID]
            max_dens = self.densest_in_chain[chainID]
            if self.densest_in_group[groupID] < max_dens:
                self.densest_in_group[groupID] = max_dens

    def _chainHOP(self):
        size = len(self.xpos)
        self.nMerge = 4
        mylog.info('Building kd tree for %d particles...' % \
            size)
        self.__init_kd_tree(size)
        # loop over the particles to find NN for each
        mylog.info('Finding nearest neighbors/density...')
        chainHOP_tags_dens()
        mylog.info('Copying results...')
        for i, nn in enumerate(self.NN):
            nn.order_index = i
            nn.NNtags = fKD.nn_tags[:,i] - 1
            nn.density = fKD.dens[i]
        # when done with the tree free the memory, del these first
        del fKD.pos, fKD.dens, fKD.nn_dist, fKD.nn_tags, fKD.mass, fKD.dens
        free_tree() # frees the kdtree object
        count = 0
        for part in self.NN:
            if part.density >= (self.threshold): count += 1
        print 'count above thresh', count
        # now each particle has NNtags/dist, and a local self density
        # let's find densest NN
        mylog.info('Finding densest nearest neighbors...')
        for i in range(size):
            self._densestNN(i, self.NN[i].NNtags)
        # mark particles in self.NN as being in/out of the domain
        for i in range(size):
            self.NN[i].is_inside = self._is_inside([self.xpos[i], self.ypos[i],
                self.zpos[i]])
        chain_count = 0
        for i, nn in enumerate(self.NN):
            if i == nn.densestNN: chain_count += 1
        mylog.info('there are %d self-densest particles' % chain_count)
        # build the chain of links
        mylog.info('Building particle chains...')
        chain_count = self._build_chains()
        # connect the chains into groups
        mylog.info('Connecting %d chains into groups...' % chain_count)
        self._connect_chains(chain_count)
        group_count = self._build_groups(chain_count)
        #mylog.info('Paring %d groups...' % group_count)
        # don't pare groups here, it has to happen later
        #self._pare_groups_by_max_dens()
        self._translate_groupIDs(group_count)
        mylog.info('Converting %d groups...' % group_count)
        # convert self.NN for returning
        gIDs = []
        dens = []
        for part in self.NN:
            gIDs.append(part.chainID)
            dens.append(part.density)
        self.gIDs = na.array(gIDs)
        self.dens = na.array(dens)


