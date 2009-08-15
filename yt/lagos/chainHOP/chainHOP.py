import math,sys
from collections import defaultdict

from yt.lagos import *
from yt.extensions.kdtree import *
from yt.performance_counters import yt_counters, time_function
from mpi4py import MPI
from Forthon import *

class RunChainHOP(ParallelAnalysisInterface):
    def __init__(self,period, padding, num_neighbors, bounds,
            xpos, ypos, zpos, index, mass, threshold=160.0, rearrange=True):
        self.threshold = threshold
        self.rearrange = rearrange
        self.saddlethresh = 2.5 * threshold
        self.peakthresh = 3 * threshold
        self.period = period
        self.padding = padding
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.xpos = xpos
        self.ypos = ypos
        self.zpos = zpos
        self.real_size = len(self.xpos)
        self.index = na.array(index, dtype='int64')
        self.mass = mass
        self.padded_particles = []
        self.nMerge = 4
        self.d = None
        yt_counters("chainHOP")
        self._chain_hop()
        yt_counters("chainHOP")

    def _global_bounds_neighbors(self):
        """
        Build a dict of the boundaries of all the tasks, and figure out which
        tasks are our neighbors.
        """
        self.mine, global_bounds = self._mpi_info_dict(self.bounds)
        # A task is a neighbor if it has a corner (vertex) inside a sphere
        # centered at my center with radius equal to my volume radius 
        # (center to a corner) plus max padding. Due to various geometrical
        # arrangements that can happen, the list of neighbors is made
        # symmetric using communication. This is not the most efficient way,
        # but it's simple and communication is not the performance bottle-
        # neck.
        my_LE, my_RE = self.bounds
        my_center = (my_LE + my_RE) / 2.
        radius = na.sqrt( (my_RE[0] - my_LE[0])**2. + \
            (my_RE[1] - my_LE[1])**2. + (my_RE[2] - my_LE[2])**2) / 2. + \
            self.max_padding
        # Put the vertices into a big list, each row is
        # array[x,y,z, taskID]
        vertices = []
        for taskID in global_bounds:
            thisLE, thisRE = global_bounds[taskID]
            vertices.append(na.array([thisLE[0], thisLE[1], thisLE[2], taskID]))
            vertices.append(na.array([thisLE[0], thisLE[1], thisRE[2], taskID]))
            vertices.append(na.array([thisLE[0], thisRE[1], thisLE[2], taskID]))
            vertices.append(na.array([thisRE[0], thisLE[1], thisLE[2], taskID]))
            vertices.append(na.array([thisLE[0], thisRE[1], thisRE[2], taskID]))
            vertices.append(na.array([thisRE[0], thisLE[1], thisRE[2], taskID]))
            vertices.append(na.array([thisRE[0], thisRE[1], thisLE[2], taskID]))
            vertices.append(na.array([thisRE[0], thisRE[1], thisRE[2], taskID]))
        # Look for vertices within one radius (periodic) of our center, and
        # add them as neighbors
        neighbors = set([])
        for vertex in vertices:
            distx = na.abs(vertex[0] - my_center[0])
            disty = na.abs(vertex[1] - my_center[1])
            distz = na.abs(vertex[2] - my_center[2])
            distx = na.minimum(distx, 1 - distx)**2.
            disty = na.minimum(disty, 1 - disty)**2.
            distz = na.minimum(distz, 1 - distz)**2.
            dist = na.sqrt(distx + disty + distz)
            if dist <= radius:
                neighbors.add(int(vertex[3]))
        # Now we build a global dict of neighbor sets, and if a remote task
        # lists us as their neighbor, we add them as our neighbor.
        self.mine, global_neighbors = self._mpi_info_dict(neighbors)
        self.neighbors = neighbors.copy()
        for taskID in global_neighbors:
            if taskID == self.mine: continue
            if self.mine in global_neighbors[taskID]:
                self.neighbors.add(taskID)
        # We can remove ourselves from the set.
        self.neighbors.discard(self.mine)
        # Clean up.
        del global_neighbors, global_bounds, vertices, neighbors
        
    def _global_padding(self):
        """
        Find the maximum padding of all our neighbors, used to send our
        annulus data.
        """
        self.mine, global_padding = self._mpi_info_dict(self.padding)
        self.max_padding = max(global_padding.itervalues())
        del global_padding

    def _communicate_raw_padding_data(self):
        """
        Send the raw particle data (x,y,zpos, mass and index) from our
        'annulus' to our neighbors. This is how each task builds up their
        padded particles. On the receive end, we discriminate against the
        data, only keeping data we want.
        """
        yt_counters("Communicate raw padding")
        (LE, RE) = self.bounds
        temp_LE = na.empty(3, dtype='float64')
        temp_RE = na.empty(3, dtype='float64')
        # Pick the particles in the annulus that will be sent.
        yt_counters("Picking padding data to send.")
        send_count = len(na.where(self.is_inside_annulus == True)[0])
        points = na.empty((send_count, 3), dtype='float64')
        points[:,0] = self.xpos[self.is_inside_annulus]
        points[:,1] = self.ypos[self.is_inside_annulus]
        points[:,2] = self.zpos[self.is_inside_annulus]
        real_indices = self.index[self.is_inside_annulus].astype('int64')
        mass = self.mass[self.is_inside_annulus].astype('float64')
        # Distribute the send sizes globally.
        global_annulus_count = {self.mine:send_count}
        global_annulus_count = self._mpi_joindict(global_annulus_count)
        yt_counters("Picking padding data to send.")
        # Initialize the arrays to receive data.
        recv_real_indices = {}
        recv_points = {}
        recv_mass = {}
        yt_counters("Initializing recv arrays.")
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            recv_real_indices[opp_neighbor] = na.empty(opp_size, dtype='int64')
            recv_points[opp_neighbor] = na.empty((opp_size,3), dtype='float64')
            recv_mass[opp_neighbor] = na.empty(opp_size, dtype='float64')
        yt_counters("Initializing recv arrays.")
        yt_counters("MPI stuff.")
        # Now we receive the data, but we don't actually use it.
        hooks = []
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            hooks.append(MPI.COMM_WORLD.Irecv([recv_real_indices[opp_neighbor], MPI.LONG], opp_neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Irecv([recv_points[opp_neighbor], MPI.DOUBLE], opp_neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Irecv([recv_mass[opp_neighbor], MPI.DOUBLE], opp_neighbor, 0))
        # Now we send the data.
        for neighbor in self.neighbors:
            hooks.append(MPI.COMM_WORLD.Isend([real_indices, MPI.LONG], neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Isend([points, MPI.DOUBLE], neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Isend([mass, MPI.DOUBLE], neighbor, 0))
        # We need to define our expanded boundaries.
        temp_LE = LE - self.padding
        temp_RE = RE + self.padding
        # Now we use the data, after all the recvs are done. This can probably,
        # and stuff below, be turned into a while loop, that exits once all the
        # data has been received and processed. The processing order doesn't
        # matter.
        MPI.Request.Waitall(hooks)
        yt_counters("MPI stuff.")
        yt_counters("Discrimination.")
        # We can now delete our sent data.
        del points, mass, real_indices
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            # Adjust the values of the positions if needed. I think there's
            # a better way to do this!
            for i,point in enumerate(recv_points[opp_neighbor][:opp_size]):
                if point[0] < temp_LE[0] and point[0] < temp_RE[0]:
                    recv_points[opp_neighbor][i][0] += self.period[0]
                if point[0] > temp_LE[0] and point[0] > temp_RE[0]:
                    recv_points[opp_neighbor][i][0] -= self.period[0]
                if point[1] < temp_LE[1] and point[1] < temp_RE[1]:
                    recv_points[opp_neighbor][i][1] += self.period[1]
                if point[1] > temp_LE[1] and point[1] > temp_RE[1]:
                    recv_points[opp_neighbor][i][1] -= self.period[1]
                if point[2] < temp_LE[2] and point[2] < temp_RE[2]:
                    recv_points[opp_neighbor][i][2] += self.period[2]
                if point[2] > temp_LE[2] and point[2] > temp_RE[2]:
                    recv_points[opp_neighbor][i][2] -= self.period[2]
            is_inside = ( (recv_points[opp_neighbor][:opp_size] >= temp_LE).all(axis=1) * \
                (recv_points[opp_neighbor][:opp_size] < temp_RE).all(axis=1) )
            self.index = na.concatenate((self.index, recv_real_indices[opp_neighbor][:opp_size][is_inside]))
            # Clean up immediately to reduce peak memory usage.
            del recv_real_indices[opp_neighbor]
            self.xpos = na.concatenate((self.xpos, recv_points[opp_neighbor][:opp_size,0][is_inside]))
            self.ypos = na.concatenate((self.ypos, recv_points[opp_neighbor][:opp_size,1][is_inside]))
            self.zpos = na.concatenate((self.zpos, recv_points[opp_neighbor][:opp_size,2][is_inside]))
            del recv_points[opp_neighbor]
            self.mass = na.concatenate((self.mass, recv_mass[opp_neighbor][:opp_size][is_inside]))
            del recv_mass[opp_neighbor]
        yt_counters("Discrimination.")
        self.size = self.index.size
        # Now that we have the full size, initialize the chainID array
        self.chainID = na.ones(self.size,dtype='int64') * -1
        # Clean up explicitly, but these should be empty dicts by now.
        del recv_real_indices, hooks
        del recv_points, recv_mass
        yt_counters("Communicate raw padding")
        

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
        fKD.rearrange = self.rearrange # True is faster, but uses more memory
        # This actually copies the data into the fortran space.
        fKD.pos[0, :] = self.xpos
        fKD.pos[1, :] = self.ypos
        fKD.pos[2, :] = self.zpos
        fKD.mass = self.mass
        # Now call the fortran.
        create_tree()

    def _is_inside(self, round):
        """
        There are three classes of particles.
        1. Particles inside the 'real' region of each subvolume.
        2. Particles ouside, added in the 'padding' for purposes of having 
           correct particle densities in the real region.
        3. Particles that are one padding distance inside the edges of the
           real region. The chainIDs of these particles are communicated
           to the neighboring tasks so chains can be merged into groups.
        The input *round* is either 'first' or 'second.' First is before the
        padded particles have been communicated, and second after.
        """
        # Test to see if the points are in the 'real' region
        (LE, RE) = self.bounds
        if round == 'first':
            points = na.empty((self.real_size, 3), dtype='float64')
        elif round == 'second':
            points = na.empty((self.size, 3), dtype='float64')
        points[:,0] = self.xpos
        points[:,1] = self.ypos
        points[:,2] = self.zpos
        self.is_inside = ( (points >= LE).all(axis=1) * \
            (points < RE).all(axis=1) )
        # Below we find out which particles are in the `annulus', one padding
        # distance inside the boundaries. First we find the particles outside
        # this inner boundary.
        if round == 'first':
            temp_LE = LE + self.max_padding
            temp_RE = RE - self.max_padding
        if round == 'second':
            temp_LE = LE + self.padding
            temp_RE = RE - self.padding
        inner = na.invert( (points >= temp_LE).all(axis=1) * \
            (points < temp_RE).all(axis=1) )
        # After inverting the logic above, we want points that are both
        # inside the real region, but within one padding of the boundary,
        # and this will do it.
        self.is_inside_annulus = na.bitwise_and(self.is_inside, inner)
        # Below we make a mapping of real particle index->local ID
        # Unf. this has to be a dict, because any task can have
        # particles of any particle_index, which means that if it were an
        # array every task would probably end up having this array be as long
        # as the full number of particles.
        # We can skip this the first time around.
        if round == 'second':
            self.rev_index = {}
            for i in xrange(int(self.size)):
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
        yt_counters("densestNN")
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
        yt_counters("densestNN")
        del chunk_NNtags, max_loc, n_dens
    
    def _build_chains(self):
        """
        Build the first round of particle chains. If the particle is too low in
        density, move on.
        """
        yt_counters("build_chains")
        chainIDmax = 0
        self.densest_in_chain = {} # chainID->part ID, one to one
        for i in xrange(int(self.size)):
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
        yt_counters("build_chains")
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
        yt_counters("globally_assign_chainIDs")
        # First find out the number of chains on each processor.
        self.mine, chain_info = self._mpi_info_dict(chain_count)
        self.nchains = sum(chain_info.itervalues())
        # Figure out our offset.
        self.my_first_id = sum([v for k,v in chain_info.iteritems() if k < self.mine])
        # Change particle IDs, -1 always means no chain assignment.
        select = (self.chainID != -1)
        select = select * self.my_first_id
        self.chainID += select
        del select
        yt_counters("globally_assign_chainIDs")

    def _create_global_densest_in_chain(self):
        """
        With the globally unique chainIDs, update densest_in_chain.
        """
        yt_counters("create_global_densest_in_chain")
        # Shift the values over.
        temp = self.densest_in_chain.copy()
        self.densest_in_chain = {}
        for dens in temp:
            self.densest_in_chain[dens + self.my_first_id] = temp[dens]
        # Distribute this.
        self.densest_in_chain = self._mpi_joindict_unpickled_float(self.densest_in_chain)
        del temp
        yt_counters("create_global_densest_in_chain")

    def _communicate_uphill_info(self):
        """
        Communicate the links to the correct neighbors from uphill_info.
        """
        yt_counters("communicate_uphill_info")
        # Find out how many particles we're going to receive, and make arrays
        # of the right size and type to store them.
        to_recv_count = 0
        temp_indices = {}
        temp_chainIDs = {}
        for opp_neighbor in self.neighbors:
            opp_size = self.global_padded_count[opp_neighbor]
            to_recv_count += opp_size
            temp_indices[opp_neighbor] = na.empty(opp_size, dtype='int64')
            temp_chainIDs[opp_neighbor] = na.empty(opp_size, dtype='int64')
        # The arrays we'll actually keep around...
        self.recv_real_indices = na.empty(to_recv_count, dtype='int64')
        self.recv_chainIDs = na.empty(to_recv_count, dtype='int64')
        # Set up the receives, but don't actually use them.
        hooks = []
        for opp_neighbor in self.neighbors:
            hooks.append(MPI.COMM_WORLD.Irecv([temp_indices[opp_neighbor], MPI.LONG], opp_neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Irecv([temp_chainIDs[opp_neighbor], MPI.LONG], opp_neighbor, 0))
        # Send padded particles to our neighbors.
        for neighbor in self.neighbors:
            hooks.append(MPI.COMM_WORLD.Isend([self.uphill_real_indices, MPI.LONG], neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Isend([self.uphill_chainIDs, MPI.LONG], neighbor, 0))
        # Now actually use the data once it's good to go.
        MPI.Request.Waitall(hooks)
        so_far = 0
        for opp_neighbor in self.neighbors:
            opp_size = self.global_padded_count[opp_neighbor]
            # Only save the part of the buffer that we want to the right places
            # in the full listing.
            self.recv_real_indices[so_far:(so_far + opp_size)] = \
                temp_indices[opp_neighbor][0:opp_size]
            self.recv_chainIDs[so_far:(so_far + opp_size)] = \
                temp_chainIDs[opp_neighbor][0:opp_size]
            so_far += opp_size
        # Clean up.
        del temp_indices, temp_chainIDs, hooks
        yt_counters("communicate_uphill_info")

    def _recurse_global_chain_links(self, chainID_translate_map_global, chainID, seen):
        """
        Step up the global chain links until we reach the self-densest chain,
        very similarly to the recursion of particles to densest nearest
        neighbors.
        """
        new_chainID = chainID_translate_map_global[chainID]
        if  new_chainID == chainID:
            return int(chainID)
        elif new_chainID in seen:
            mylog.info('seen %s' % str(seen))
            for s in seen:
                mylog.info('%d %d' % (s, chainID_translate_map_global[s]))
        else:
            seen.append(new_chainID)
            return self._recurse_global_chain_links(chainID_translate_map_global, new_chainID, seen)
        

    def _connect_chains_across_tasks(self):
        """
        Using the uphill links of chains, chains are linked across boundaries.
        Chains that link to a remote chain are recorded, and a complete dict
        of chain connections is created, globally. Then chainIDs are
        reassigned recursively, assigning the ID of the most dense chainID
        to every chain that links to it.
        """
        yt_counters("connect_chains_across_tasks")
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
            # -1 above. Also, since links are supposed to go only uphill,
            # ensure that they are being recorded that way below.
            if localID != -1:
                if self.recv_chainIDs[i] != self.chainID[localID]:# and \
                    #self.densest_in_chain[self.chainID[localID]] > self.densest_in_chain[self.recv_chainIDs[i]]:
                    chainID_translate_map_local[self.recv_chainIDs[i]] = \
                        self.chainID[localID]
        # In chainID_translate_map_local, chains may
        # 'point' to only one chain, but a chain may have many that point to
        # it. Therefore each key (a chain) in this dict is unique, but the items
        # the keys point to are not necessarily unique. Most chains do not
        # have an entry in the dict, and that is addressed after the 
        # communication step.
        chainID_translate_map_global = \
            self._mpi_joindict_unpickled_int(chainID_translate_map_local)
        for i in xrange(int(self.nchains)):
            try:
                target = chainID_translate_map_global[i]
            except KeyError:
                # The the chain is either self most-dense, or a singleton.
                chainID_translate_map_global[i] = int(i)
        # Build a list of chain densities, sorted smallest to largest.
        dens_temp = []
        for key in self.densest_in_chain:
            item = [self.densest_in_chain[key], key]
            dens_temp.append(item)
        dens_temp.sort()
        # Loop over chains, smallest to largest density, recursively until
        # we reach a self-assigned chain. Then we assign that final chainID to
        # the *current* one only.
        for chain in dens_temp:
            seen = []
            seen.append(int(chain[1]))
            new_chainID = \
                self._recurse_global_chain_links(chainID_translate_map_global, int(chain[1]), seen)
            chainID_translate_map_global[int(chain[1])] = new_chainID
            # At the same time, remove chains from densest_in_chain that have
            # been reassigned.
            if chain[1] != new_chainID:
                del self.densest_in_chain[chain[1]]
                # Also fix nchains to keep up.
                self.nchains -= 1
        # Convert local particles to their new chainID
        for i in xrange(int(self.size)):
            old_chainID = self.chainID[i]
            if old_chainID == -1: continue
            new_chainID = chainID_translate_map_global[old_chainID]
            self.chainID[i] = new_chainID
        del chainID_translate_map_local, dens_temp, self.recv_chainIDs
        del self.recv_real_indices, self.uphill_real_indices, self.uphill_chainIDs
        yt_counters("connect_chains_across_tasks")

    def _communicate_annulus_chainIDs(self):
        """
        Transmit all of our chainID-ed particles that are within self.padding
        of the boundaries to all of our neighbors. Tests show that this is
        faster than trying to figure out which of the neighbors to send the data
        to.
        """
        yt_counters("communicate_annulus_chainIDs")
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
        # Set up the receiving arrays.
        recv_real_indices = {}
        recv_chainIDs = {}
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            recv_real_indices[opp_neighbor] = na.empty(opp_size, dtype='int64')
            recv_chainIDs[opp_neighbor] = na.empty(opp_size, dtype='int64')
        # Set up the receving hooks.
        hooks = []
        for opp_neighbor in self.neighbors:
            hooks.append(MPI.COMM_WORLD.Irecv([recv_real_indices[opp_neighbor], MPI.LONG], opp_neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Irecv([recv_chainIDs[opp_neighbor], MPI.LONG], opp_neighbor, 0))
        # Now we send them.
        for neighbor in self.neighbors:
            hooks.append(MPI.COMM_WORLD.Isend([real_indices, MPI.LONG], neighbor, 0))
            hooks.append(MPI.COMM_WORLD.Isend([chainIDs, MPI.LONG], neighbor, 0))
        # Now we use them when they're nice and ripe.
        MPI.Request.Waitall(hooks)
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            # Update our local data.
            for i,real_index in enumerate(recv_real_indices[opp_neighbor][0:opp_size]):
                try:
                    localID = self.rev_index[real_index]
                    # We are only updating our particles that are in our
                    # padding, so to be rigorous we will skip particles
                    # that are in our real region.
                    if self.is_inside[localID]:
                        continue
                    self.chainID[localID] = recv_chainIDs[opp_neighbor][i]
                except KeyError:
                    # We ignore data that's not for us.
                    continue
        # Clean up.
        del recv_real_indices, recv_chainIDs, real_indices, chainIDs, select
        del hooks
        yt_counters("communicate_annulus_chainIDs")

    def _connect_chains(self):
        """
        With the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        yt_counters("connect_chains")
        self.chain_densest_n = {} # chainID -> {chainIDs->boundary dens}
        self.reverse_map = {} # chainID -> groupID, one to one
        for chainID in self.densest_in_chain:
            #self.chain_densest_n[int(chainID)] = {}
            self.reverse_map[int(chainID)] = -1
        for i in xrange(int(self.size)):
            # Don't consider this particle if it's not part of a chain.
            if self.chainID[i] < 0: continue
            # If this particle is in the padding, don't make a connection.
            if not self.is_inside[i]: continue
            # Find this particle's chain max_dens.
            part_max_dens = self.densest_in_chain[self.chainID[i]]
            # Loop over nMerge closest nearest neighbors.
            for j in xrange(int(self.nMerge+2)):
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
                    # Make sure that the higher density chain has an entry.
                    try:
                        test = self.chain_densest_n[higher_chain]
                    except KeyError:
                        self.chain_densest_n[higher_chain] = {}
                    # See if this boundary density is higher than
                    # previously recorded for this pair of chains.
                    # Links only go one direction.
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
        yt_counters("connect_chains")

    def _make_global_chain_densest_n(self):
        """
        We want to record the maximum boundary density between all chains on
        all tasks.
        """
        yt_counters("make_global_chain_densest_n")
        self.chain_densest_n = \
                self._mpi_maxdict_dict(self.chain_densest_n)
        yt_counters("make_global_chain_densest_n")
    
    def _build_groups(self):
        """
        With the collection of possible chain links, build groups.
        """
        yt_counters("build_groups")
        #if self.mine == 0:
        g_high = []
        g_low = []
        g_dens = []
        densestbound = {} # chainID -> boundary density
        for chainID in self.densest_in_chain:
            densestbound[chainID] = -1.0
        groupID = 0
        # First assign a group to all chains with max_dens above peakthresh.
        # The initial groupIDs will be assigned with decending peak density.
        # This guarantees that the group with the smaller groupID is the
        # higher chain, as in chain_high below.
        def ksort(d):
            temp = []
            for key in d:
                item = [key, d[key]]
                temp.append(item)
            temp.sort(lambda x,y: -cmp(x[1],y[1]))
            temp = na.array(temp)
            # return the first column, which are the dict keys sorted by
            # the second column, the densities.
            return temp[:,0]
        group_equivalancy_map = defaultdict(set)
        for chainID in ksort(self.densest_in_chain):
            if self.densest_in_chain[chainID] >= self.peakthresh:
                self.reverse_map[chainID] = groupID
                groupID += 1
        big_count = 0
        for chain_high in self.chain_densest_n:
            for chain_low in self.chain_densest_n[chain_high]:
                big_count += 1
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
                        #for chID in self.reverse_map:
                        #    if self.reverse_map[chID] == group_low:
                        #        self.reverse_map[chID] = group_high
                        if group_low != group_high:
                            group_equivalancy_map[group_low].add(group_high)
                            group_equivalancy_map[group_high].add(group_low)
                        continue
                # Else, one is above peakthresh, the other below
                # find out if this is the densest boundary seen so far for
                # the lower chain.
                group_high = self.reverse_map[chain_high]
                if dens > densestbound[chain_low]:
                    densestbound[chain_low] = dens
                    self.reverse_map[chain_low] = group_high
                # Done double loop over links.
        # Now refactor group_equivalancy_map back into reverse_map. The group
        # mapping may be more than one link long, so we need to do it
        # recursively. The best way to think about this is a field full of 
        # rabbit holes. The holes are connected at nexuses at the surface.
        # Each groupID (key) in group_equivalancy_map represents a hole, and
        # the values the nexuses are the tunnels lead to. The tunnels are two-way,
        # and when you go through it, you block the passage through that
        # tunnel in that direction, so you don't repeat yourself later. You can
        # go back through that tunnel, but your search ends there because all
        # the other tunnels have been closed at the old nexus. In this fashion your search 
        # spreads out like the water shooting out of the ground in 'Caddy
        # Shack.'
        mylog.info("here 1")
        Set_list = []
        # We only want the holes that are modulo mine.
        keys = na.array(group_equivalancy_map.keys(), dtype='int64')
        if self._mpi_get_size() == None:
            size = 1
        else:
            size = self._mpi_get_size()
        select = (keys % size == self.mine)
        groupIDs = keys[select]
        mine_groupIDs = set([]) # Records only ones modulo mine.
        not_mine_groupIDs = set([]) # All the others.
        for groupID in groupIDs:
            if groupID in mine_groupIDs:
                continue
            mine_groupIDs.add(groupID)
            current_sets = []
            new_set = group_equivalancy_map[groupID]
            final_set = new_set.copy()
            while len(new_set) > 0:
                to_add_set = set([])
                liter = new_set.difference(mine_groupIDs).difference(not_mine_groupIDs)
                new_mine, new_other = [], []
                for link_gID in liter:
                    to_add_set.update(group_equivalancy_map[link_gID])
                    if link_gID % size == self.mine:
                        new_mine.append(link_gID)
                    else:
                        new_other.append(link_gID)
                mine_groupIDs.update(new_mine)
                not_mine_groupIDs.update(new_other)
                final_set.update(to_add_set)
                new_set = to_add_set
            # Make sure it's not empty
            final_set.add(groupID)
            Set_list.append(final_set)
        # Convert this list of sets into a look-up table
        lookup = {}
        for i,item in enumerate(Set_list):
            item_min = min(item)
            for groupID in item:
                lookup[groupID] = item_min
        mylog.info("here 2")
        # To bring it all together, join the dicts.
        lookup = self._mpi_joindict_unpickled_int(lookup)
        mylog.info("here 3")
        # Now apply this to reverse_map
        for chainID in self.reverse_map:
            groupID = self.reverse_map[chainID]
            if groupID == -1:
                continue
            try:
                self.reverse_map[chainID] = lookup[groupID]
            except KeyError:
                continue
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
        temp = list(set(temp))
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
        #else:
        #    group_count = 0
        #group_count = MPI.COMM_WORLD.bcast(group_count, root=0)
        #self.reverse_map = MPI.COMM_WORLD.bcast(self.reverse_map, root=0)
        #self.reverse_map = self._mpi_bcast_int_dict_unpickled(self.reverse_map)
        mylog.info("here 4")
        yt_counters("build_groups")
        return group_count

    def _translate_groupIDs(self, group_count):
        """
        Using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        yt_counters("translate_groupIDs")
        for i in xrange(int(self.size)):
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
        for i in xrange(int(group_count)):
            self.densest_in_group[i] = 0.0
        for chainID in self.densest_in_chain:
            groupID = self.reverse_map[chainID]
            if groupID == -1: continue
            max_dens = self.densest_in_chain[chainID]
            if self.densest_in_group[groupID] < max_dens:
                self.densest_in_group[groupID] = max_dens
        yt_counters("translate_groupIDs")

    def _precompute_group_info(self):
        yt_counters("Precomp.")
        """
        For all groups, compute the various global properties, except bulk
        velocity, to save time in HaloFinding.py (fewer barriers!).
        """
        # First we need to find the maximum density point for all groups.
        # I think this will be faster than several vector operations that need
        # to pull the entire chainID array out of memory several times.
        max_dens_point = na.zeros((self.group_count,4),dtype='float64')
        for part in xrange(int(self.size)):
            if self.chainID[part] == -1: continue
            groupID = self.chainID[part]
            if self.density[part] == self.densest_in_group[groupID]:
                max_dens_point[groupID] = na.array([self.density[part], \
                self.xpos[part], self.ypos[part], self.zpos[part]], dtype='float64')
        # Now we broadcast this, effectively, with an allsum. Even though
        # some groups are on multiple tasks, there is only one densest_in_chain
        # and only that task contributed above.
        self.max_dens_point = self._mpi_Allsum_float(max_dens_point)
        # Now CoM.
        c_vec = self.max_dens_point[:,1:4] - na.array([0.5,0.5,0.5])
        size = na.zeros(self.group_count, dtype='int64')
        CoM_M = na.zeros((self.group_count,3),dtype='float64')
        Tot_M = na.zeros(self.group_count, dtype='float64')
        # This keeps track of what we own so we don't divide by zero later...
        I_own = []
        for part in xrange(int(self.size)):
            if self.chainID[part] == -1: continue
            groupID = self.chainID[part]
            size[groupID] += 1
            if groupID not in I_own:
                I_own.append(groupID)
            loc = na.array([self.xpos[part], self.ypos[part], \
                self.zpos[part]], dtype='float64') - c_vec[groupID]
            loc = loc - na.floor(loc)
            CoM_M[groupID] += (loc * self.mass[part])
            Tot_M[groupID] += self.mass[part]
        for groupID in xrange(int(self.group_count)):
            if groupID in I_own:
                CoM_M[groupID] /= Tot_M[groupID]
                CoM_M[groupID] += c_vec[groupID]
                CoM_M[groupID] *= Tot_M[groupID]
        # Now we find their global values
        self.group_sizes = self._mpi_Allsum_int(size)
        CoM_M = self._mpi_Allsum_float(CoM_M)
        self.Tot_M = self._mpi_Allsum_float(Tot_M)
        self.CoM = na.empty((self.group_count,3), dtype='float64')
        for groupID in xrange(int(self.group_count)):
            self.CoM[groupID] = CoM_M[groupID] / self.Tot_M[groupID]
        # Now we find the maximum radius for all groups.
        max_radius = na.zeros(self.group_count, dtype='float64')
        for part in xrange(int(self.size)):
            groupID = self.chainID[part]
            if groupID == -1: continue
            loc = na.array([self.xpos[part], self.ypos[part], self.zpos[part]])
            dist = na.sqrt(\
                na.minimum(na.abs(self.CoM[groupID][0] - loc[0]), self.period[0] - na.abs(self.CoM[groupID][0] - loc[0]))**2 \
                + na.minimum(na.abs(self.CoM[groupID][1] - loc[1]), self.period[1] - na.abs(self.CoM[groupID][1] - loc[1]))**2 \
                + na.minimum(na.abs(self.CoM[groupID][2] - loc[2]), self.period[2] - na.abs(self.CoM[groupID][2] - loc[2]))**2)
            if dist > max_radius[groupID]:
                max_radius[groupID] = dist
        # Find the maximum across all tasks.
        self.max_radius = self._mpi_float_array_max(max_radius)
        yt_counters("Precomp.")

    def _chain_hop(self):
        mylog.info("Running Python version.")
        self._global_padding()
        self._global_bounds_neighbors()
        mylog.info("Distributing padded particles...")
        self._is_inside('first')
        self._communicate_raw_padding_data()
        mylog.info('Building kd tree for %d particles...' % \
            self.size)
        self._init_kd_tree()
        # Mark particles in as being in/out of the domain.
        self._is_inside('second')
        # Loop over the particles to find NN for each.
        mylog.info('Finding nearest neighbors/density...')
        yt_counters("chainHOP_tags_dens")
        chainHOP_tags_dens()
        yt_counters("chainHOP_tags_dens")
        self.density = fKD.dens
        self.NNtags = (fKD.nn_tags - 1).transpose()
        # We can free these right now, the rest later.
        del fKD.dens, fKD.nn_tags, fKD.mass, fKD.dens
        #count = len(na.where(self.density >= self.threshold)[0])
        #print 'count above thresh', count
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
        mylog.info('Precomputing info for %d groups...' % group_count)
        self._precompute_group_info()
        # We need to fix chainID and density because HaloFinding is expecting
        # an array only as long as the real data.
        self.chainID = self.chainID[:self.real_size]
        self.density = self.density[:self.real_size]
        del self.mass, self.xpos, self.ypos, self.zpos
