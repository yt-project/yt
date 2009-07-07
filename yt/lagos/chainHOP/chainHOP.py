import math #, cPickle
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
        by linking particles that have high enough average combined density.
        """
        self.chain_connections = {} # groupID -> [chainIDs], one to many
        self.reverse_map = {} # chainID -> groupID, one to one
        self.connected_groups = [] # list of lists of groupIDs that are connected
                                   # each entry is a pair [ , ]
        for i in range(chain_count):
            self.reverse_map[i] = -1
        groupID = 0 # inreased with every newly linked chain pair
        for part in self.NN:
            # don't consider this particle if it's not part of a chain
            if part.chainID < 0: continue
            # if this particle is in the padding, don't make a connection
            if part.is_inside is False: continue
            # find this particles chain max_dens
            part_dens = self.densest_in_chain[part.chainID]
            # loop over nMerge nearest neighbors
            for i in range(self.nMerge+2):
                thisNN = part.NNtags[i]
                # if our neighbor is in the same chain, move on
                if part.chainID == self.NN[thisNN].chainID: continue
                # no introspection, nor our connected NN
                if thisNN==part.order_index or thisNN==part.densestNN: continue
                # if our neighbor is not part of a group, move on
                if self.NN[thisNN].chainID < 0: continue
                # find thisNN's chain's max_dens
                thisNN_dens = self.densest_in_chain[self.NN[thisNN].chainID]
                # if the average density is high enough, and they're both
                # themselves groups (dens>peakdens), link the chains. Chains
                # with both dens > peakdens are *only* linked if their boundary
                # dens > saddle. Otherwise they remain separate for all time.
                if ((self.NN[thisNN].density + part.density) / 2.) >= (self.threshold * 2.5) and
                    part_dens >= (3*self.threshold) and thisNN_dens >= (3*self.threshold):
                    # find out if either particle chainID has already been mapped
                    group1 = self.reverse_map[part.chainID]
                    group2 = self.reverse_map[self.NN[thisNN].chainID]
                    # if they're already connected, move on
                    if group1 == group2 and group1 != -1:
                        #print 'already connected'
                        continue
                    # if one chain has been connected into a group, assign
                    # the other
                    if group1 > -1 and group2 == -1:
                        self.chain_connections[group1].append(self.NN[thisNN].chainID)
                        self.reverse_map[self.NN[thisNN].chainID] = group1
                        #print 'group1 already assigned'
                        continue
                    # visa-versa
                    if group2 > -1 and group1 == -1:
                        self.chain_connections[group2].append(part.chainID)
                        self.reverse_map[part.chainID] = group2
                        #print 'group2 already assigned'
                        continue
                    # link the different groups
                    if group1 > -1 and group2 > -1:
                        # max/min for unique representation
                        # this will eventually be reversed in _connect_chains()
                        connection = [max(group1,group2),min(group1,group2)]
                        #connection = [min(group1,group2),max(group1,group2)]
                        if connection not in self.connected_groups:
                            self.connected_groups.append(connection)
                        #print 'both already assigned'
                        continue
                    # if neither chain has been linked into a group yet, we've
                    # come down this far
                    self.chain_connections[groupID] = [part.chainID,self.NN[thisNN].chainID]
                    self.reverse_map[part.chainID] = groupID
                    self.reverse_map[self.NN[thisNN].chainID] = groupID
                    #print 'neither already assigned'
                    groupID += 1
                # end if avg dens & peak densities are high enough.
                
                # now we need to record *possible* connections, replacing old
                # entries if the new one is higher in boundary threshold
                # This is analagous to the particle chain hopping, where we
                # need to find the densest nearest neighbors, but here we're
                # finding the densest nearest neighbor chains.
                
        # chains that haven't been linked to another chain increase the count
        # by themselves and change label. Also add to chain_connections
        for i in range(chain_count):
            if self.reverse_map[i] == -1:
                self.reverse_map[i] = groupID
                groupID += 1
        return groupID
    
    def _connect_groups(self):
        # sort connected_groups decending by first group, which is always larger
        # this is to prevent pathological situations of multiply linked groups
        # unf., Cython doesn't support sorting by a key so we do it manually
        index_map = []
        for index, pair in enumerate(self.connected_groups):
            datum = [pair[0], pair[1], index]
            insort(index_map,datum)
        connected_groups_temp = self.connected_groups[:]
        for index, triplet in enumerate(index_map):
            connected_groups_temp[index] = self.connected_groups[triplet[2]]
        self.connected_groups = connected_groups_temp[:]
        self.connected_groups.reverse()
        # below is the way we would sort this if it worked in Cython.
        # self.connected_groups.sort(key=lambda x: -1*x[0])
        # now reverse the entries
        for pair in self.connected_groups:
            pair.reverse()
        # now update the chain_connections, merging groups
        chain_connections_temp = self.chain_connections.copy()
        added = []
        for groups in self.connected_groups:
            # always add the larger indexed group of chains into the smaller indexed
            chain_connections_temp[groups[0]].extend(self.chain_connections[groups[1]])
            added.append(groups[1])
        # copy the temp back, skipping over added groups of chains
        fresh_count = 0
        self.chain_connections = {}
        for groupID in chain_connections_temp:
            if groupID in added: continue
            self.chain_connections[fresh_count] = chain_connections_temp[groupID]
            # update reverse_map at the same time
            for chainID in chain_connections_temp[groupID]:
                self.reverse_map[chainID] = fresh_count
            fresh_count += 1
        return fresh_count
    
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
    
    def _translate_groupIDs(self):
        """
        using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        for i in range(len(self.NN)):
            # don't translate non-affiliated particles
            if self.NN[i].chainID == -1: continue
            self.NN[i].chainID = self.reverse_map[self.NN[i].chainID]
        # also translate self.densest_in_chain
        temp = {}
        for groupID, connections in self.chain_connections.items():
            temp[groupID] = max((self.densest_in_chain[chainID] for chainID
                           in connections))
        # copy it back
        self.densest_in_chain = temp

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
        group_count = self._connect_chains(chain_count)
        mylog.info('Connecting %d groups...' % group_count)
        group_count = self._connect_groups()
        #mylog.info('Paring %d groups...' % group_count)
        # don't pare groups here, it has to happen later
        #self._pare_groups_by_max_dens()
        self._translate_groupIDs()
        mylog.info('Converting %d groups...' % group_count)
        # convert self.NN for returning
        gIDs = []
        dens = []
        for part in self.NN:
            gIDs.append(part.chainID)
            dens.append(part.density)
        self.gIDs = na.array(gIDs)
        self.dens = na.array(dens)


