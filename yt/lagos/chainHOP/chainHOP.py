import math #, cPickle
from bisect import insort

#from yt.lagos.kd import *
from yt.extensions.kdtree import *

from Forthon import *

class partNN:
    def __init__(self, NNlist, density, densestNN, chainID, is_inside, order_index):
        self.NNlist = NNlist # a list of [distance, order_index (not particle_index)]
        self.density = density # this particles local density
        self.densestNN = densestNN # this particles densest nearest neighbor
        self.chainID = chainID # this particles current chain identifier, changes
        self.is_inside = is_inside # True/False is inside the main volume
        self.order_index = order_index # order_index for this particle

class Run_chain_HOP:
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
        fKD.tags = empty((self.num_neighbors),dtype='i')
        fKD.dist = empty((self.num_neighbors),dtype='f')
        fKD.pos = empty((3,size),dtype='f')
        fKD.qv = empty(3,dtype='f')
        fKD.nn = self.num_neighbors
        fKD.nparts = size
        fKD.sort = True
        fKD.rearrange = True
        for i in range(size):
            fKD.pos[0][i] = self.xpos[i]
            fKD.pos[1][i] = self.ypos[i]
            fKD.pos[2][i] = self.zpos[i]
        # now call the fortran
        create_tree()

    def __initialize_NN(self):
        size = len(self.xpos)
        self.NN = []
        for i in range(size):
            NNtemp = partNN([],0.0,-1,-1,False,-1)
            self.NN.append(NNtemp)
    
    def _is_inside(self,point):
        # test to see if this point is inside the 'real' region
        (LE,RE) = self.bounds
        point = na.array(point)
        if (point >= LE).all() and (point < RE).all():
            return True
        else:
            return False
    
    def _smDensitySym(self, pi, NNlist):
        # calculate the density for particle pi
        # this is giving different values than HOP because the py-kd tree is
        # giving slightly different distances between particles than the c-kd
        # tree.
        ih2 = 4.0/NNlist[self.num_neighbors-1][0]
        fNorm = 0.5*math.sqrt(ih2)*ih2/math.pi
        sum = self.mass.sum()
        for i,neighbor in enumerate(NNlist):
            pj = neighbor[1]
            r2 = neighbor[0]*ih2
            rs = 2.0 - math.sqrt(r2)
            if (r2 < 1.0): rs = (1.0 - 0.75*rs*r2)
            else: rs = 0.25*rs*rs*rs
            rs *= fNorm
            self.NN[pi].density += rs * self.mass[pj] / sum
            self.NN[pj].density += rs * self.mass[pi] / sum
    
    def _densestNN(self, pi, NNlist):
        # find the densest nearest neighbor
        NNdens = 0.
        for neighbor in NNlist:
            if self.NN[neighbor[1]].density > NNdens:
                NNdens = self.NN[neighbor[1]].density
                NNindex = neighbor[1]
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
            # loop over nearest neighbors
            for i,neighbor in enumerate(part.NNlist):
                # no introspection, nor our connected NN
                if i==part.order_index or neighbor[1]==part.densestNN: continue
                # if our neighbor is not part of a group, move on
                if self.NN[neighbor[1]].chainID < 0: continue
                # if our neighbor is in the padding, move on
                #inside = self.NN[neighbor[1]].is_inside
                #if inside is False: continue # I think this might be OK
                # if the average density is high enough, link the chains
                if ((self.NN[neighbor[1]].density + part.density) / 2.) >= (self.threshold * 2.5):
                    # find out if either particle chainID has already been mapped
                    chain1 = self.reverse_map[part.chainID]
                    chain2 = self.reverse_map[self.NN[neighbor[1]].chainID]
                    # if they're already connected, move on
                    if chain1 == chain2 and chain1 != -1:
                        #print 'already connected'
                        continue
                    # if one chain has been connected into a group, assign
                    # the other
                    if chain1 > -1 and chain2 == -1:
                        self.chain_connections[chain1].append(self.NN[neighbor[1]].chainID)
                        self.reverse_map[self.NN[neighbor[1]].chainID] = chain1
                        #print 'chain1 already assigned'
                        continue
                    # visa-versa
                    if chain2 > -1 and chain1 == -1:
                        self.chain_connections[chain2].append(part.chainID)
                        self.reverse_map[part.chainID] = chain2
                        #print 'chain2 already assigned'
                        continue
                    # link the different groups
                    if chain1 > -1 and chain2 > -1:
                        # max/min for unique representation
                        # this will eventually be reversed in _connect_chains()
                        connection = [max(chain1,chain2),min(chain1,chain2)]
                        #connection = [min(chain1,chain2),max(chain1,chain2)]
                        if connection not in self.connected_groups:
                            self.connected_groups.append(connection)
                        #print 'both already assigned'
                        continue
                    # if neither chain has been linked into a group yet, we've
                    # come down this far
                    self.chain_connections[groupID] = [part.chainID,self.NN[neighbor[1]].chainID]
                    self.reverse_map[part.chainID] = groupID
                    self.reverse_map[self.NN[neighbor[1]].chainID] = groupID
                    #print 'neither already assigned'
                    groupID += 1
        # chains that haven't been linked to another chain should just link back
        # to themselves
        for i in range(chain_count):
            if self.reverse_map[i] == -1:
                self.reverse_map[i] = i
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
            max_dens = 0.0
            for chainID in self.chain_connections[groupID]:
                if max_dens < self.densest_in_chain[chainID]:
                    max_dens = self.densest_in_chain[chainID]
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
        for part in self.NN:
            # don't translate non-affiliated particles
            if part.chainID == -1: continue
            part.chainID = self.reverse_map[part.chainID]
        # also translate self.densest_in_chain
        temp = {}
        for groupID in self.chain_connections:
            max_dens = 0.
            for chainID in self.chain_connections[groupID]:
                if max_dens < self.densest_in_chain[chainID]:
                    max_dens = self.densest_in_chain[chainID]
            temp[groupID] = max_dens
        # copy it back
        self.densest_in_chain = temp

    def _chainHOP(self):
        size = len(self.xpos)
        mylog.info('Building kd tree for %d particles...' % \
            size)
        self.__init_kd_tree(size)
        # loop over the particles to find NN for each
        mylog.info('Finding nearest neighbors/density...')
        for i in range(size):
            self.NN[i].order_index = i
            if i % 1000==0:
                mylog.info('loop %d of %d' % (i,size))
            # call the fortran after setting the query vector
            fKD.qv = [ self.xpos[i], self.ypos[i], self.zpos[i] ]
            find_nn_nearest_neighbors()
            tags = fKD.tags - 1 # -1 for fortran counting
            dist = fKD.dist
            n_points = []
            for n in range(self.num_neighbors):
                n_points.append([dist[n],tags[n]])
            self.NN[i].NNlist = n_points
            # find the density at this particle
            self._smDensitySym(i,self.NN[i].NNlist)
        # when done with the tree free the memory
        free_tree()
        count = 0
        for part in self.NN:
            if part.density >= (self.threshold): count += 1
        print 'count above thresh', count
        # now each particle has NNlist, and a local self density
        # let's find densest NN
        mylog.info('Finding densest nearest neighbors...')
        for i in range(size):
            self._densestNN(i,self.NN[i].NNlist)
        # mark particles in self.NN as being in/out of the domain
        for i in range(size):
            self.NN[i].is_inside = self._is_inside([self.xpos[i], self.ypos[i],
                self.zpos[i]])
        chain_count = 0
        for i in range(size):
            if i == self.NN[i].densestNN:
                chain_count += 1
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



#     def _chainHOP(self):
#         self.dataset = []
#         size = len(self.xpos)
#         for i in range(size):
#             p = Point()
#             p.data = [self.xpos[i],
#                       self.ypos[i],
#                       self.zpos[i]]
#             p.order_index = i
#             self.dataset.append(p)
#         mylog.info('Building kd tree for %d particles...' % \
#             size)
#         kd = buildKdHyperRectTree(self.dataset[:],self.num_neighbors)
#         # loop over the particles to find NN for each
#         mylog.info('Finding nearest neighbors/density...')
#         for i in self.dataset:
#             pi = i.order_index
#             self.NN[pi].order_index = pi
#             if pi % 1000==0:
#                 mylog.info('loop %d of %d' % (pi,size))
#             # make the neighbors object
#             neighbors = Neighbors()
#             neighbors.k = self.num_neighbors # +1 for the central point
#             # the search radius being larger than the padding doesn't help
#             # padded particles very close to the real reigion, so making it only
#             # slightly larger is all that is necessary, and that's set in
#             # HaloFinding.py.
#             neighbors.minDistanceSquared = self.search_radius * self.search_radius
#             neighbors.points = []
#             getKNN(i.data, kd, neighbors, 0., self.period.tolist())
#             n_points = []
#             for n in neighbors.points:
#                 n_points.append([n[0],n[1].order_index])
#                 if pi==0:
#                     print n[0],n[1].order_index
#             self.NN[pi].NNlist = n_points
#             # find the density at this particle
#             self._smDensitySym(pi,self.NN[pi].NNlist)
#         count = 0
#         for part in self.NN:
#             if part.density >= (self.threshold): count += 1
#         print 'count above thresh', count
#         # now each particle has NNlist, and a local self density
#         # let's find densest NN
#         mylog.info('Finding densest nearest neighbors...')
#         for i in self.dataset:
#             pi = i.order_index
#             self._densestNN(pi,self.NN[pi].NNlist)
# #         fp = open('NN.txt','w')
# #         for part in self.NN:
# #             data = self.dataset[part.order_index].data
# #             line = str(data) + ' ' + str(self.dataset[part.densestNN].data) + '\n'
# #             fp.write(line)
# #         fp.close()
#         # mark particles in self.NN as being in/out of the domain
#         for i in self.dataset:
#             self.NN[i.order_index].is_inside = self._is_inside(i.data)
#         chain_count = 0
#         for i in self.dataset:
#             pi = i.order_index
#             if pi == self.NN[pi].densestNN:
#                 chain_count += 1
#         mylog.info('there are %d self-densest particles' % chain_count)
#         # build the chain of links
#         mylog.info('Building particle chains...')
#         chain_count = self._build_chains()
#         # connect the chains into groups
#         mylog.info('Connecting %d chains into groups...' % chain_count)
#         group_count = self._connect_chains(chain_count)
#         mylog.info('Connecting %d groups...' % group_count)
#         group_count = self._connect_groups()
#         #mylog.info('Paring %d groups...' % group_count)
#         # don't pare groups here, it has to happen later
#         #self._pare_groups_by_max_dens()
#         self._translate_groupIDs()
#         mylog.info('Converting %d groups...' % group_count)
#         # convert self.NN for returning
#         gIDs = []
#         dens = []
#         for part in self.NN:
#             gIDs.append(part.chainID)
#             dens.append(part.density)
#         self.gIDs = na.array(gIDs)
#         self.dens = na.array(dens)
# 
