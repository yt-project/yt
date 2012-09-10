"""
A implementation of the HOP algorithm that runs in parallel.

Author: Stephen Skory <s@skory.us>
Affiliation: UCSD/CASS
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Stephen Skory.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from collections import defaultdict
import itertools, sys
import numpy as np
import gc

from yt.funcs import *
from yt.utilities.performance_counters import yt_counters, time_function
try:
    from yt.utilities.kdtree import \
        chainHOP_tags_dens, \
        create_tree, fKD, find_nn_nearest_neighbors, \
        free_tree, find_chunk_nearest_neighbors
except ImportError:
    mylog.debug("The Fortran kD-Tree did not import correctly.")

from yt.utilities.spatial import cKDTree

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_blocking_call, \
    ParallelAnalysisInterface

class ParallelHOPHaloFinder(ParallelAnalysisInterface):
    def __init__(self,period, padding, num_neighbors, bounds,
            particle_fields, threshold=160.0, rearrange=True,
            premerge=True, tree='F'):
        ParallelAnalysisInterface.__init__(self)
        self.threshold = threshold
        self.rearrange = rearrange
        self.premerge = premerge
        self.saddlethresh = 2.5 * threshold
        self.peakthresh = 3 * threshold
        self.period = period
        self.padding = padding
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.xpos = particle_fields.pop("particle_position_x")
        self.ypos = particle_fields.pop("particle_position_y")
        self.zpos = particle_fields.pop("particle_position_z")
        self.real_size = len(self.xpos)
        self.index = particle_fields.pop("particle_index")
        self.mass = particle_fields.pop("ParticleMassMsun")
        self.padded_particles = []
        self.nMerge = 4
        self.tree = tree
        yt_counters("chainHOP")
        self.max_mem = 0
        self.__max_memory()
        self._chain_hop()
        yt_counters("chainHOP")

    def _global_bounds_neighbors(self):
        """
        Build a dict of the boundaries of all the tasks, and figure out which
        tasks are our geometric neighbors.
        """
        self.neighbors = set([])
        self.mine, global_bounds = self.comm.mpi_info_dict(self.bounds)
        my_LE, my_RE = self.bounds
        # Put the vertices into a big list, each row is
        # array[x,y,z, taskID]
        vertices = []
        my_vertices = []
        for taskID in global_bounds:
            thisLE, thisRE = global_bounds[taskID]
            if self.mine != taskID:
                vertices.append(np.array([thisLE[0], thisLE[1], thisLE[2], taskID]))
                vertices.append(np.array([thisLE[0], thisLE[1], thisRE[2], taskID]))
                vertices.append(np.array([thisLE[0], thisRE[1], thisLE[2], taskID]))
                vertices.append(np.array([thisRE[0], thisLE[1], thisLE[2], taskID]))
                vertices.append(np.array([thisLE[0], thisRE[1], thisRE[2], taskID]))
                vertices.append(np.array([thisRE[0], thisLE[1], thisRE[2], taskID]))
                vertices.append(np.array([thisRE[0], thisRE[1], thisLE[2], taskID]))
                vertices.append(np.array([thisRE[0], thisRE[1], thisRE[2], taskID]))
            if self.mine == taskID:
                my_vertices.append(np.array([thisLE[0], thisLE[1], thisLE[2]]))
                my_vertices.append(np.array([thisLE[0], thisLE[1], thisRE[2]]))
                my_vertices.append(np.array([thisLE[0], thisRE[1], thisLE[2]]))
                my_vertices.append(np.array([thisRE[0], thisLE[1], thisLE[2]]))
                my_vertices.append(np.array([thisLE[0], thisRE[1], thisRE[2]]))
                my_vertices.append(np.array([thisRE[0], thisLE[1], thisRE[2]]))
                my_vertices.append(np.array([thisRE[0], thisRE[1], thisLE[2]]))
                my_vertices.append(np.array([thisRE[0], thisRE[1], thisRE[2]]))
        # Find the neighbors we share corners with. Yes, this is lazy with
        # a double loop, but it works and this is definitely not a performance
        # bottleneck.
        for my_vertex in my_vertices:
            for vertex in vertices:
                if vertex[3] in self.neighbors: continue
                # If the corners touch, it's easy. This is the case if resizing
                # (load-balancing) is turned off.
                if (my_vertex % self.period == vertex[0:3] % self.period).all():
                    self.neighbors.add(int(vertex[3]))
                    continue
                # Also test to see if the distance to this corner is within
                # max_padding, which is more likely the case with load-balancing
                # turned on.
                dx = min( np.fabs(my_vertex[0] - vertex[0]), \
                    self.period[0] - np.fabs(my_vertex[0] - vertex[0]))
                dy = min( np.fabs(my_vertex[1] - vertex[1]), \
                    self.period[1] - np.fabs(my_vertex[1] - vertex[1]))
                dz = min( np.fabs(my_vertex[2] - vertex[2]), \
                    self.period[2] - np.fabs(my_vertex[2] - vertex[2]))
                d = np.sqrt(dx*dx + dy*dy + dz*dz)
                if d <= self.max_padding:
                    self.neighbors.add(int(vertex[3]))
        # Faces and edges.
        for dim in range(3):
            dim1 = (dim + 1) % 3
            dim2 = (dim + 2) % 3
            left_face = my_LE[dim]
            right_face = my_RE[dim]
            for taskID in global_bounds:
                if taskID == self.mine or taskID in self.neighbors: continue
                thisLE, thisRE = global_bounds[taskID]
                max1 = max(my_LE[dim1], thisLE[dim1])
                max2 = max(my_LE[dim2], thisLE[dim2])
                min1 = min(my_RE[dim1], thisRE[dim1])
                min2 = min(my_RE[dim2], thisRE[dim2])
                # Faces.
                # First, faces that touch directly.
                if (thisRE[dim] == left_face or thisRE[dim]%self.period[dim] == left_face) and \
                        max1 <= min1 and max2 <= min2:
                    self.neighbors.add(taskID)
                    continue
                elif (thisLE[dim] == right_face or thisLE[dim] == right_face%self.period[dim]) and \
                        max1 <= min1 and max2 <= min2:
                    self.neighbors.add(taskID)
                    continue
                # If an intervening subvolume has a width less than the padding
                # (rare, but possible), a neighbor may not actually touch, so
                # we need to account for that.
                if (abs(thisRE[dim] - left_face) <= self.max_padding or \
                        abs(thisRE[dim]%self.period[dim] - left_face) <= self.max_padding) and \
                        max1 <= min1 and max2 <= min2:
                    self.neighbors.add(taskID)
                    continue
                elif (abs(thisLE[dim] - right_face) <= self.max_padding or \
                        abs(thisLE[dim] - right_face%self.period[dim]) <= self.max_padding) and \
                        max1 <= min1 and max2 <= min2:
                    self.neighbors.add(taskID)
                    continue
                # Edges.
                # First, edges that touch.
                elif (my_LE[dim] == (thisRE[dim]%self.period[dim]) and \
                        my_LE[dim1] == (thisRE[dim1]%self.period[dim1]) and \
                        max2 <= min2) or \
                        (my_LE[dim] == (thisRE[dim]%self.period[dim]) and \
                        my_LE[dim2] == (thisRE[dim2]%self.period[dim2]) and \
                        max1 <= min1):
                    self.neighbors.add(taskID)
                    continue
                elif ((my_RE[dim]%self.period[dim]) == thisLE[dim] and \
                        (my_RE[dim1]%self.period[dim1]) == thisLE[dim1] and \
                        max2 <= min2) or \
                        ((my_RE[dim]%self.period[dim]) == thisLE[dim] and \
                        (my_RE[dim2]%self.period[dim2]) == thisLE[dim2] and \
                        max1 <= min1):
                    self.neighbors.add(taskID)
                    continue
                # Now edges that don't touch, but are close.
                if (abs(my_LE[dim] - thisRE[dim]%self.period[dim]) <= self.max_padding and \
                        abs(my_LE[dim1] - thisRE[dim1]%self.period[dim1]) <= self.max_padding and \
                        max2 <= min2) or \
                        (abs(my_LE[dim] - thisRE[dim]%self.period[dim]) <= self.max_padding and \
                        abs(my_LE[dim2] - thisRE[dim2]%self.period[dim2]) <= self.max_padding and \
                        max1 <= min1):
                    self.neighbors.add(taskID)
                    continue
                elif (abs(my_RE[dim]%self.period[dim] - thisLE[dim]) <= self.max_padding and \
                        abs(my_RE[dim1]%self.period[dim1] - thisLE[dim1]) <= self.max_padding and \
                        max2 <= min2) or \
                        (abs(my_RE[dim]%self.period[dim] - thisLE[dim]) <= self.max_padding and \
                        abs(my_RE[dim2]%self.period[dim2] - thisLE[dim2]) <= self.max_padding and \
                        max1 <= min1):
                    self.neighbors.add(taskID)
                    continue
        # Now we build a global dict of neighbor sets, and if a remote task
        # lists us as their neighbor, we add them as our neighbor. This is 
        # probably not needed because the stuff above should be symmetric,
        # but it isn't a big issue.
        self.mine, global_neighbors = self.comm.mpi_info_dict(self.neighbors)
        for taskID in global_neighbors:
            if taskID == self.mine: continue
            if self.mine in global_neighbors[taskID]:
                self.neighbors.add(taskID)
        # We can remove ourselves from the set if it got added somehow.
        self.neighbors.discard(self.mine)
        # Clean up.
        del global_neighbors, global_bounds, vertices, my_vertices
        
    def _global_padding(self, round):
        """
        Find the maximum padding of all our neighbors, used to send our
        annulus data.
        """
        if round == 'first':
            max_pad = np.max(self.padding)
            self.mine, self.global_padding = self.comm.mpi_info_dict(max_pad)
            self.max_padding = max(self.global_padding.itervalues())
        elif round == 'second':
            self.max_padding = 0.
            for neighbor in self.neighbors:
                self.max_padding = np.maximum(self.global_padding[neighbor], \
                    self.max_padding)

    def _communicate_padding_data(self):
        """
        Send the particles each of my neighbors need to build up their padding.
        """
        yt_counters("Communicate discriminated padding")
        # First build a global dict of the padded boundaries of all the tasks.
        (LE, RE) = self.bounds
        (LE_padding, RE_padding) = self.padding
        temp_LE = LE - LE_padding
        temp_RE = RE + RE_padding
        expanded_bounds = (temp_LE, temp_RE)
        self.mine, global_exp_bounds = self.comm.mpi_info_dict(expanded_bounds)
        send_real_indices = {}
        send_points = {}
        send_mass = {}
        send_size = {}
        # This will reduce the size of the loop over particles.
        yt_counters("Picking padding data to send.")
        send_count = self.is_inside_annulus.sum()
        points = np.empty((send_count, 3), dtype='float64')
        points[:,0] = self.xpos[self.is_inside_annulus]
        points[:,1] = self.ypos[self.is_inside_annulus]
        points[:,2] = self.zpos[self.is_inside_annulus]
        real_indices = self.index[self.is_inside_annulus].astype('int64')
        mass = self.mass[self.is_inside_annulus].astype('float64')
        # Make the arrays to send.
        shift_points = points.copy()
        for neighbor in self.neighbors:
            temp_LE, temp_RE = global_exp_bounds[neighbor]
            for i in xrange(3):
                left = ((points[:,i] < temp_LE[i]) * (points[:,i] < temp_RE[i])) * self.period[i]
                right = ((points[:,i] > temp_LE[i]) * (points[:,i] > temp_RE[i])) * self.period[i]
                shift_points[:,i] = points[:,i] + left - right
            is_inside = ( (shift_points >= temp_LE).all(axis=1) * \
                (shift_points < temp_RE).all(axis=1) )
            send_real_indices[neighbor] = real_indices[is_inside].copy()
            send_points[neighbor] = shift_points[is_inside].copy()
            send_mass[neighbor] = mass[is_inside].copy()
            send_size[neighbor] = is_inside.sum()
        del points, shift_points, mass, real_indices
        yt_counters("Picking padding data to send.")
        # Communicate the sizes to send.
        self.mine, global_send_count = self.comm.mpi_info_dict(send_size)
        del send_size
        # Initialize the arrays to receive data.
        yt_counters("Initalizing recv arrays.")
        recv_real_indices = {}
        recv_points = {}
        recv_mass = {}
        recv_size = 0
        for opp_neighbor in self.neighbors:
            opp_size = global_send_count[opp_neighbor][self.mine]
            recv_real_indices[opp_neighbor] = np.empty(opp_size, dtype='int64')
            recv_points[opp_neighbor] = np.empty((opp_size, 3), dtype='float64')
            recv_mass[opp_neighbor] = np.empty(opp_size, dtype='float64')
            recv_size += opp_size
        yt_counters("Initalizing recv arrays.")
        # Setup the receiving slots.
        yt_counters("MPI stuff.")
        hooks = []
        for opp_neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_recv(recv_real_indices[opp_neighbor], opp_neighbor))
            hooks.append(self.comm.mpi_nonblocking_recv(recv_points[opp_neighbor], opp_neighbor))
            hooks.append(self.comm.mpi_nonblocking_recv(recv_mass[opp_neighbor], opp_neighbor))
        # Let's wait here to be absolutely sure that all the receive buffers
        # have been created before any sending happens!
        self.comm.barrier()
        # Now we send the data.
        for neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_send(send_real_indices[neighbor], neighbor))
            hooks.append(self.comm.mpi_nonblocking_send(send_points[neighbor], neighbor))
            hooks.append(self.comm.mpi_nonblocking_send(send_mass[neighbor], neighbor))
        # Now we use the data, after all the comms are done.
        self.comm.mpi_Request_Waitall(hooks)
        yt_counters("MPI stuff.")
        yt_counters("Processing padded data.")
        del send_real_indices, send_points, send_mass
        # Now we add the data to ourselves.
        self.index_pad = np.empty(recv_size, dtype='int64')
        self.xpos_pad = np.empty(recv_size, dtype='float64')
        self.ypos_pad = np.empty(recv_size, dtype='float64')
        self.zpos_pad = np.empty(recv_size, dtype='float64')
        self.mass_pad = np.empty(recv_size, dtype='float64')
        so_far = 0
        for opp_neighbor in self.neighbors:
            opp_size = global_send_count[opp_neighbor][self.mine]
            self.index_pad[so_far:so_far+opp_size] = recv_real_indices[opp_neighbor]
            # Clean up immediately to reduce peak memory usage.
            del recv_real_indices[opp_neighbor]
            self.xpos_pad[so_far:so_far+opp_size] = recv_points[opp_neighbor][:,0]
            self.ypos_pad[so_far:so_far+opp_size] = recv_points[opp_neighbor][:,1]
            self.zpos_pad[so_far:so_far+opp_size] = recv_points[opp_neighbor][:,2]
            del recv_points[opp_neighbor]
            self.mass_pad[so_far:so_far+opp_size] = recv_mass[opp_neighbor]
            del recv_mass[opp_neighbor]
            so_far += opp_size
        yt_counters("Processing padded data.")
        # The KDtree node search wants the particles to be in the full box,
        # not the expanded dimensions of shifted (<0 or >1 generally) volume,
        # so we fix the positions of particles here.
        yt_counters("Flipping coordinates around the periodic boundary.")
        self.xpos_pad = self.xpos_pad % self.period[0]
        self.ypos_pad = self.ypos_pad % self.period[1]
        self.zpos_pad = self.zpos_pad % self.period[2]
        yt_counters("Flipping coordinates around the periodic boundary.")
        self.size = self.index.size + self.index_pad.size
        # Now that we have the full size, initialize the chainID array
        self.chainID = np.ones(self.size,dtype='int64') * -1
        # Clean up explicitly, but these should be empty dicts by now.
        del recv_real_indices, hooks, recv_points, recv_mass
        yt_counters("Communicate discriminated padding")

    def _init_kd_tree(self):
        """
        Set up the data objects that get passed to the kD-tree code.
        """
        yt_counters("init kd tree")
        if self.tree == 'F':
            # Yes, we really do need to initialize this many arrays.
            # They're deleted in _parallelHOP.
            fKD.dens = np.zeros(self.size, dtype='float64', order='F')
            fKD.mass = np.concatenate((self.mass, self.mass_pad))
            del self.mass
            fKD.pos = np.empty((3, self.size), dtype='float64', order='F')
            # This actually copies the data into the fortran space.
            self.psize = self.xpos.size
            fKD.pos[0, :self.psize] = self.xpos
            fKD.pos[1, :self.psize] = self.ypos
            fKD.pos[2, :self.psize] = self.zpos
            del self.xpos, self.ypos, self.zpos
            gc.collect()
            fKD.pos[0, self.psize:] = self.xpos_pad
            fKD.pos[1, self.psize:] = self.ypos_pad
            fKD.pos[2, self.psize:] = self.zpos_pad
            del self.xpos_pad, self.ypos_pad, self.zpos_pad
            gc.collect()
            fKD.qv = np.asfortranarray(np.empty(3, dtype='float64'))
            fKD.nn = self.num_neighbors
            # Plus 2 because we're looking for that neighbor, but only keeping 
            # nMerge + 1 neighbor tags, skipping ourselves.
            fKD.nMerge = self.nMerge + 2
            fKD.nparts = self.size
            fKD.sort = True # Slower, but needed in _connect_chains
            fKD.rearrange = self.rearrange # True is faster, but uses more memory
            # Now call the fortran.
            create_tree(0)
        elif self.tree == 'C':
            self.mass = np.concatenate((self.mass, self.mass_pad))
            self.pos = np.empty((self.size, 3), dtype='float64')
            self.psize = self.xpos.size
            self.pos[:self.psize, 0] = self.xpos
            self.pos[:self.psize, 1] = self.ypos
            self.pos[:self.psize, 2] = self.zpos
            del self.xpos, self.ypos, self.zpos
            gc.collect()
            self.pos[self.psize:, 0] = self.xpos_pad
            self.pos[self.psize:, 1] = self.ypos_pad
            self.pos[self.psize:, 2] = self.zpos_pad
            del self.xpos_pad, self.ypos_pad, self.zpos_pad
            gc.collect()
            self.kdtree = cKDTree(self.pos, leafsize = 64)
        self.__max_memory()
        yt_counters("init kd tree")

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
            points = np.empty((self.real_size, 3), dtype='float64')
            points[:,0] = self.xpos
            points[:,1] = self.ypos
            points[:,2] = self.zpos
            self.is_inside = ( (points >= LE).all(axis=1) * \
                (points < RE).all(axis=1) )
        elif round == 'second':
            if self.tree == 'F':
                self.is_inside = ( (fKD.pos.T >= LE).all(axis=1) * \
                    (fKD.pos.T < RE).all(axis=1) )
            elif self.tree == 'C':
                self.is_inside = ( (self.pos > LE).all(axis=1) * \
                    (self.pos < RE).all(axis=1) )
        # Below we find out which particles are in the `annulus', one padding
        # distance inside the boundaries. First we find the particles outside
        # this inner boundary.
        temp_LE = LE + self.max_padding
        temp_RE = RE - self.max_padding
        if round == 'first':
            inner = np.invert( (points >= temp_LE).all(axis=1) * \
                (points < temp_RE).all(axis=1) )
        elif round == 'second' or round == 'third':
            if self.tree == 'F':
                inner = np.invert( (fKD.pos.T >= temp_LE).all(axis=1) * \
                    (fKD.pos.T < temp_RE).all(axis=1) )
            elif self.tree == 'C':
                inner = np.invert( (self.pos >= temp_LE).all(axis=1) * \
                    (self.pos < temp_RE).all(axis=1) )
        if round == 'first':
            del points
        # After inverting the logic above, we want points that are both
        # inside the real region, but within one padding of the boundary,
        # and this will do it.
        self.is_inside_annulus = np.bitwise_and(self.is_inside, inner)
        del inner
        # Below we make a mapping of real particle index->local ID
        # Unf. this has to be a dict, because any task can have
        # particles of any particle_index, which means that if it were an
        # array every task would probably end up having this array be as long
        # as the full number of particles.
        # We can skip this the first two times around.
        if round == 'third':
            temp = np.arange(self.size)
            my_part = np.bitwise_or(np.invert(self.is_inside), self.is_inside_annulus)
            my_part = np.bitwise_and(my_part, (self.chainID != -1))
            catted_indices = np.concatenate(
                (self.index, self.index_pad))[my_part]
            self.rev_index = dict.fromkeys(catted_indices)
            self.rev_index.update(itertools.izip(catted_indices, temp[my_part]))
            del my_part, temp, catted_indices
        self.__max_memory()

    def _densestNN(self):
        """
        For all particles, find their densest nearest neighbor. It is done in
        chunks to keep the memory usage down.
        The first search of nearest neighbors (done earlier) did not return all 
        num_neighbor neighbors, so we need to do it again, but we're not
        keeping the all of this data, just using it.
        """
        yt_counters("densestNN")
        self.densestNN = np.empty(self.size,dtype='int64')
        # We find nearest neighbors in chunks.
        chunksize = 10000
        if self.tree == 'F':
            fKD.chunk_tags = np.asfortranarray(np.empty((self.num_neighbors, chunksize), dtype='int64'))
            start = 1 # Fortran counting!
            finish = 0
            while finish < self.size:
                finish = min(finish+chunksize,self.size)
                # Call the fortran. start and finish refer to the data locations
                # in fKD.pos, and specify the range of particles to find nearest
                # neighbors
                fKD.start = start
                fKD.finish = finish
                find_chunk_nearest_neighbors()
                chunk_NNtags = (fKD.chunk_tags[:,:finish-start+1] - 1).transpose()
                # Find the densest nearest neighbors by referencing the already
                # calculated density.
                n_dens = np.take(self.density,chunk_NNtags)
                max_loc = np.argmax(n_dens,axis=1)
                for i in xrange(finish - start + 1): # +1 for fortran counting.
                    j = start + i - 1 # -1 for fortran counting.
                    self.densestNN[j] = chunk_NNtags[i,max_loc[i]]
                start = finish + 1
        elif self.tree == 'C':
            start = 0
            finish = 0
            while finish < self.size - 1:
                finish = min(finish+chunksize, self.size)
                # Unlike above, this function returns a new chunk_NNtags
                # that is the right size every time. But this may not actually
                # be as memory efficient - fragmenting?
                chunk_NNtags = self.kdtree.find_chunk_nearest_neighbors(start, \
                    finish, num_neighbors=self.num_neighbors)
                n_dens = np.take(self.density, chunk_NNtags)
                max_loc = np.argmax(n_dens, axis=1)
                max_loc = np.argmax(n_dens,axis=1)
                for i in xrange(finish - start):
                    j = start + i
                    self.densestNN[j] = chunk_NNtags[i,max_loc[i]]
                start = finish
        yt_counters("densestNN")
        self.__max_memory()
        del chunk_NNtags, max_loc, n_dens
    
    def _build_chains(self):
        """
        Build the first round of particle chains. If the particle is too low in
        density, move on.
        """
        yt_counters("build_chains")
        chainIDmax = 0
        self.densest_in_chain = np.ones(10000, dtype='float64') * -1 # chainID->density, one to one
        self.densest_in_chain_real_index = np.ones(10000, dtype='int64') * -1 # chainID->real_index, one to one
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
        self.padded_particles = np.array(self.padded_particles, dtype='int64')
        self.densest_in_chain = self.__clean_up_array(self.densest_in_chain)
        self.densest_in_chain_real_index = self.__clean_up_array(self.densest_in_chain_real_index)
        yt_counters("build_chains")
        self.__max_memory()
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
            self.densest_in_chain = self.__add_to_array(self.densest_in_chain,
                chainIDmax, self.density[pi], 'float64')
            if pi < self.real_size:
                self.densest_in_chain_real_index = self.__add_to_array(self.densest_in_chain_real_index,
                chainIDmax, self.index[pi], 'int64')
            else:
                self.densest_in_chain_real_index = self.__add_to_array(self.densest_in_chain_real_index,
                chainIDmax, self.index_pad[pi-self.real_size], 'int64')
            # if this is a padded particle, record it for later
            if not inside:
                self.padded_particles.append(pi)
            return chainIDmax
        # Otherwise, recursively link to nearest neighbors.
        else:
            chainIDnew = self._recurse_links(nn, chainIDmax)
            self.chainID[pi] = chainIDnew
            return chainIDnew

    def _recurse_preconnected_links(self, chain_map, thisID):
        if min(thisID, min(chain_map[thisID])) == thisID:
            return thisID
        else:
            return self._recurse_preconnected_links(chain_map, min(chain_map[thisID]))

    def _preconnect_chains(self, chain_count):
        """
        In each subvolume, chains that share a boundary that both have high
        enough peak densities are prelinked in order to reduce the size of the
        global chain objects. This is very similar to _connect_chains().
        """
        # First we'll sort them, which will be used below.
        mylog.info("Locally sorting chains...")
        yt_counters("preconnect_chains")
        yt_counters("local chain sorting.")
        sort = self.densest_in_chain.argsort()
        sort = np.flipud(sort)
        map = np.empty(sort.size,dtype='int64')
        map[sort] = np.arange(sort.size)
        self.densest_in_chain = self.densest_in_chain[sort]
        self.densest_in_chain_real_index = self.densest_in_chain_real_index[sort]
        del sort
        for i,chID in enumerate(self.chainID):
            if chID == -1: continue
            self.chainID[i] = map[chID]
        del map
        yt_counters("local chain sorting.")
        mylog.info("Preconnecting %d chains..." % chain_count)
        chain_map = defaultdict(set)
        for i in xrange(max(self.chainID)+1):
            chain_map[i].add(i)
        yt_counters("preconnect kd tree search.")
        if self.tree == 'C':
            nn = self.nMerge + 2
            rv = self.kdtree.chainHOP_preconnect(
                self.chainID, self.density, self.densest_in_chain,
                self.is_inside, self.search_again,
                self.peakthresh, self.saddlethresh, nn, self.nMerge,
                chain_map)
            self.search_again = rv.astype("bool")
            yt_counters("preconnect kd tree search.")
        elif self.tree == 'F':
            # Plus 2 because we're looking for that neighbor, but only keeping 
            # nMerge + 1 neighbor tags, skipping ourselves.
            fKD.dist = np.empty(self.nMerge+2, dtype='float64')
            fKD.tags = np.empty(self.nMerge+2, dtype='int64')
            # We can change this here to make the searches faster.
            fKD.nn = self.nMerge + 2
            for i in xrange(self.size):
                # Don't consider this particle if it's not part of a chain.
                if self.chainID[i] < 0: continue
                chainID_i = self.chainID[i]
                # If this particle is in the padding, don't make a connection.
                if not self.is_inside[i]: continue
                # Find this particle's chain max_dens.
                part_max_dens = self.densest_in_chain[chainID_i]
                # We're only connecting >= peakthresh chains now.
                if part_max_dens < self.peakthresh: continue
                # Loop over nMerge closest nearest neighbors.
                if self.tree == 'F':
                    fKD.qv = fKD.pos[:, i]
                    find_nn_nearest_neighbors()
                    NNtags = fKD.tags[:] - 1
                elif self.tree == 'C':
                    qv = self.pos[i, :]
                    NNtags = self.kdtree.query(qv, nn)[1]
                same_count = 0
                for j in xrange(int(self.nMerge+1)):
                    thisNN = NNtags[j+1] # Don't consider ourselves at NNtags[0]
                    thisNN_chainID = self.chainID[thisNN]
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
                        thisNN_max_dens = self.densest_in_chain[thisNN_chainID]
                        # We're only linking peakthresh chains
                        if thisNN_max_dens < self.peakthresh: continue
                        # Calculate the two groups boundary density.
                        boundary_density = (self.density[thisNN] + self.density[i]) / 2.
                        # Don't connect if the boundary is too low.
                        if boundary_density < self.saddlethresh: continue
                        # Mark these chains as related.
                        chain_map[thisNN_chainID].add(chainID_i)
                        chain_map[chainID_i].add(thisNN_chainID)
                if same_count == self.nMerge + 1:
                    # All our neighbors are in the same chain already, so 
                    # we don't need to search again.
                    self.search_again[i] = False
            try:
                del NNtags
            except UnboundLocalError:
                pass
        yt_counters("preconnect kd tree search.")
        # Recursively jump links until we get to a chain whose densest
        # link is to itself. At that point we've found the densest chain
        # in this set of sets and we keep a record of that.
        yt_counters("preconnect pregrouping.")
        final_chain_map = np.empty(max(self.chainID)+1, dtype='int64')
        removed = 0
        for i in xrange(self.chainID.max()+1):
            j = chain_count - i - 1
            densest_link = self._recurse_preconnected_links(chain_map, j)
            final_chain_map[j] = densest_link
            if j != densest_link:
                removed += 1
                self.densest_in_chain[j] = -1
                self.densest_in_chain_real_index[j] = -1
        del chain_map
        for i in xrange(self.size):
            if self.chainID[i] != -1:
                self.chainID[i] = final_chain_map[self.chainID[i]]
        del final_chain_map
        # Now make the chainID assignments consecutive.
        map = np.empty(self.densest_in_chain.size, dtype='int64')
        dic_new = np.empty(chain_count - removed, dtype='float64')
        dicri_new = np.empty(chain_count - removed, dtype='int64')
        new = 0
        for i,dic in enumerate(self.densest_in_chain):
            if dic > 0:
                map[i] = new
                dic_new[new] = dic
                dicri_new[new] = self.densest_in_chain_real_index[i]
                new += 1
            else:
                map[i] = -1
        for i in range(self.size):
            if self.chainID[i] != -1:
                self.chainID[i] = map[self.chainID[i]]
        del map
        self.densest_in_chain = dic_new
        self.densest_in_chain_real_index = dicri_new
        self.__max_memory()
        yt_counters("preconnect pregrouping.")
        mylog.info("Preconnected %d chains." % removed)
        yt_counters("preconnect_chains")

        return chain_count - removed

    def _globally_assign_chainIDs(self, chain_count):
        """
        Convert local chainIDs into globally unique chainIDs.
        """
        yt_counters("globally_assign_chainIDs")
        # First find out the number of chains on each processor.
        self.mine, chain_info = self.comm.mpi_info_dict(chain_count)
        self.nchains = sum(chain_info.values())
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
        # Shift the values over effectively by concatenating them in the same
        # order as the values have been shifted in _globally_assign_chainIDs()
        yt_counters("global chain MPI stuff.")
        self.densest_in_chain = self.comm.par_combine_object(self.densest_in_chain,
                datatype="array", op="cat")
        self.densest_in_chain_real_index = self.comm.par_combine_object(
                self.densest_in_chain_real_index,
                datatype="array", op="cat")
        yt_counters("global chain MPI stuff.")
        # Sort the chains by density here. This is an attempt to make it such
        # that the merging stuff in a few steps happens in the same order
        # all the time.
        mylog.info("Sorting chains...")
        yt_counters("global chain sorting.")
        sort = self.densest_in_chain.argsort()
        sort = np.flipud(sort)
        map = np.empty(sort.size,dtype='int64')
        map[sort] =np.arange(sort.size)
        self.densest_in_chain = self.densest_in_chain[sort]
        self.densest_in_chain_real_index = self.densest_in_chain_real_index[sort]
        del sort
        for i,chID in enumerate(self.chainID):
            if chID == -1: continue
            self.chainID[i] = map[chID]
        del map
        yt_counters("global chain sorting.")
        # For some reason chains that share the most-dense particle are not
        # being linked, so we link them 'by hand' here.
        mylog.info("Pre-linking chains 'by hand'...")
        yt_counters("global chain hand-linking.")
        # If there are no repeats, we can skip this mess entirely.
        uniq = np.unique(self.densest_in_chain_real_index)
        if uniq.size != self.densest_in_chain_real_index.size:
            # Find only the real particle indices that are repeated to reduce
            # the dict workload below.
            dicri = self.densest_in_chain_real_index[self.densest_in_chain_real_index.argsort()]
            diff = np.ediff1d(dicri)
            diff = (diff == 0) # Picks out the places where the ids are equal
            diff = np.concatenate((diff, [False])) # Makes it the same length
            # This has only the repeated IDs. Sets are faster at searches than
            # arrays.
            dicri = set(dicri[diff])
            reverse = defaultdict(set)
            # Here we find a reverse mapping of real particle ID to chainID
            for chainID, real_index in enumerate(self.densest_in_chain_real_index):
                if real_index in dicri:
                    reverse[real_index].add(chainID)
            del dicri, diff
            # If the real index has len(set)>1, there are multiple chains that need
            # to be linked
            tolink = defaultdict(set)
            for real in reverse:
                if len(reverse[real]) > 1:
                    # Unf. can't slice a set, so this will have to do.
                    tolink[min(reverse[real])] = reverse[real]
                    tolink[min(reverse[real])].discard(min(reverse[real]))
            del reverse
            # Now we will remove the other chains from the dicts and re-assign
            # particles to their new chainID.
            fix_map = {}
            for tokeep in tolink:
                for remove in tolink[tokeep]:
                    fix_map[remove] = tokeep
                    self.densest_in_chain[remove] = -1.0
                    self.densest_in_chain_real_index[remove] = -1
            for i, chainID in enumerate(self.chainID):
                try:
                    new = fix_map[chainID]
                except KeyError:
                    continue
                self.chainID[i] = new
            del tolink, fix_map
        del uniq
        yt_counters("global chain hand-linking.")
        yt_counters("create_global_densest_in_chain")

    def _communicate_uphill_info(self):
        """
        Communicate the links to the correct neighbors from uphill_info.
        """
        yt_counters("communicate_uphill_info")
        # Find out how many particles we're going to receive, and make arrays
        # of the right size and type to store them.
        to_recv_count = 0
        temp_indices = dict.fromkeys(self.neighbors)
        temp_chainIDs = dict.fromkeys(self.neighbors)
        for opp_neighbor in self.neighbors:
            opp_size = self.global_padded_count[opp_neighbor]
            to_recv_count += opp_size
            temp_indices[opp_neighbor] = np.empty(opp_size, dtype='int64')
            temp_chainIDs[opp_neighbor] = np.empty(opp_size, dtype='int64')
        # The arrays we'll actually keep around...
        self.recv_real_indices = np.empty(to_recv_count, dtype='int64')
        self.recv_chainIDs = np.empty(to_recv_count, dtype='int64')
        # Set up the receives, but don't actually use them.
        hooks = []
        for opp_neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_recv(temp_indices[opp_neighbor], opp_neighbor))
            hooks.append(self.comm.mpi_nonblocking_recv(temp_chainIDs[opp_neighbor], opp_neighbor))
        # Make sure all the receive buffers are set before continuing.
        self.comm.barrier()
        # Send padded particles to our neighbors.
        for neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_send(self.uphill_real_indices, neighbor))
            hooks.append(self.comm.mpi_nonblocking_send(self.uphill_chainIDs, neighbor))
        # Now actually use the data once it's good to go.
        self.comm.mpi_Request_Waitall(hooks)
        self.__max_memory()
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
            # Bad things are about to happen if this condition is met! The
            # padding probably needs to be increased (using the safety factor).
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
        chainID_translate_map_local = np.arange(self.nchains, dtype='int64')
        # Build the stuff to send.
        self.uphill_real_indices = np.concatenate((
            self.index, self.index_pad))[self.padded_particles]
        self.uphill_chainIDs = self.chainID[self.padded_particles]
        del self.padded_particles
        # Now we make a global dict of how many particles each task is
        # sending.
        self.global_padded_count = {self.mine:self.uphill_chainIDs.size}
        self.global_padded_count = self.comm.par_combine_object(
                self.global_padded_count, datatype = "dict", op = "join")
        # Send/receive 'em.
        self._communicate_uphill_info()
        del self.global_padded_count
        self.__max_memory()
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
            if localID != -1 and self.chainID[localID] != -1:
                if self.recv_chainIDs[i] != self.chainID[localID] and \
                        self.densest_in_chain[self.chainID[localID]] >= self.densest_in_chain[self.recv_chainIDs[i]] and \
                        self.densest_in_chain[self.chainID[localID]] != -1.0 and \
                        self.densest_in_chain[self.recv_chainIDs[i]] != -1.0:
                    chainID_translate_map_local[self.recv_chainIDs[i]] = \
                        self.chainID[localID]
        self.__max_memory()
        # In chainID_translate_map_local, chains may
        # 'point' to only one chain, but a chain may have many that point to
        # it. Therefore each key (a chain) in this dict is unique, but the items
        # the keys point to are not necessarily unique.
        chainID_translate_map_global = \
            self.comm.mpi_allreduce(chainID_translate_map_local, op='min',
            dtype='int64')
        # Loop over chains, smallest to largest density, recursively until
        # we reach a self-assigned chain. Then we assign that final chainID to
        # the *current* one only.
        seen = []
        for key, density in enumerate(self.densest_in_chain):
            if density == -1: continue # Skip 'deleted' chains
            seen = []
            seen.append(key)
            new_chainID = \
                self._recurse_global_chain_links(chainID_translate_map_global, key, seen)
            chainID_translate_map_global[key] = new_chainID
            # At the same time, remove chains from densest_in_chain that have
            # been reassigned.
            if key != new_chainID:
                self.densest_in_chain[key] = -1.0
                self.densest_in_chain_real_index[key] = -1
                # Also fix nchains to keep up.
                self.nchains -= 1
        self.__max_memory()
        # Convert local particles to their new chainID
        for i in xrange(int(self.size)):
            old_chainID = self.chainID[i]
            if old_chainID == -1: continue
            new_chainID = chainID_translate_map_global[old_chainID]
            self.chainID[i] = new_chainID
        del chainID_translate_map_local, self.recv_chainIDs
        del self.recv_real_indices, self.uphill_real_indices, self.uphill_chainIDs
        del seen, chainID_translate_map_global
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
        real_indices = np.concatenate(
            (self.index, self.index_pad))[self.is_inside_annulus]
        chainIDs = self.chainID[self.is_inside_annulus]
        # We're done with this here.
        del self.is_inside_annulus
        # Eliminate un-assigned particles.
        select = (chainIDs != -1)
        real_indices = real_indices[select]
        chainIDs = chainIDs[select]
        send_count = real_indices.size
        # Here distribute the counts globally. Unfortunately, it's a barrier(), 
        # but there's so many places in this that need to be globally synched
        # that it's not worth the effort right now to make this one spot better.
        global_annulus_count = {self.mine:send_count}
        global_annulus_count = self.comm.par_combine_object(
                global_annulus_count, datatype = "dict", op = "join")
        # Set up the receiving arrays.
        recv_real_indices = dict.fromkeys(self.neighbors)
        recv_chainIDs = dict.fromkeys(self.neighbors)
        for opp_neighbor in self.neighbors:
            opp_size = global_annulus_count[opp_neighbor]
            recv_real_indices[opp_neighbor] = np.empty(opp_size, dtype='int64')
            recv_chainIDs[opp_neighbor] = np.empty(opp_size, dtype='int64')
        # Set up the receving hooks.
        hooks = []
        for opp_neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_recv(recv_real_indices[opp_neighbor], opp_neighbor))
            hooks.append(self.comm.mpi_nonblocking_recv(recv_chainIDs[opp_neighbor], opp_neighbor))
        # Make sure the recv buffers are set before continuing.
        self.comm.barrier()
        # Now we send them.
        for neighbor in self.neighbors:
            hooks.append(self.comm.mpi_nonblocking_send(real_indices, neighbor))
            hooks.append(self.comm.mpi_nonblocking_send(chainIDs, neighbor))
        # Now we use them when they're nice and ripe.
        self.comm.mpi_Request_Waitall(hooks)
        self.__max_memory()
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
        del hooks, global_annulus_count
        # We're done with this here.
        del self.rev_index
        yt_counters("communicate_annulus_chainIDs")


    def _connect_chains(self):
        """
        With the set of particle chains, build a mapping of connected chainIDs
        by finding the highest boundary density neighbor for each chain. Some
        chains will have no neighbors!
        """
        yt_counters("connect_chains")
        self.chain_densest_n = {} # chainID -> {chainIDs->boundary dens}
        # Plus 2 because we're looking for that neighbor, but only keeping 
        # nMerge + 1 neighbor tags, skipping ourselves.
        if self.tree == 'F':
            fKD.dist = np.empty(self.nMerge+2, dtype='float64')
            fKD.tags = np.empty(self.nMerge+2, dtype='int64')
            # We can change this here to make the searches faster.
            fKD.nn = self.nMerge+2
        elif self.tree == 'C':
            nn = self.nMerge + 2
        for i in xrange(int(self.size)):
            # Don't consider this particle if it's not part of a chain.
            if self.chainID[i] < 0: continue
            # If this particle is in the padding, don't make a connection.
            if not self.is_inside[i]: continue
            # Make sure that we should search this particle again.
            if not self.search_again[i]: continue
            # Find this particle's chain max_dens.
            part_max_dens = self.densest_in_chain[self.chainID[i]]
            # Make sure we're skipping deleted chains.
            if part_max_dens == -1.0: continue
            # Loop over nMerge closest nearest neighbors.
            if self.tree == 'F':
                fKD.qv = fKD.pos[:, i]
                find_nn_nearest_neighbors()
                NNtags = fKD.tags[:] - 1
            elif self.tree == 'C':
                qv = self.pos[i, :]
                NNtags = self.kdtree.query(qv, nn)[1]
            for j in xrange(int(self.nMerge+1)):
                thisNN = NNtags[j+1] # Don't consider ourselves at NNtags[0]
                thisNN_chainID = self.chainID[thisNN]
                # If our neighbor is in the same chain, move on.
                if self.chainID[i] == thisNN_chainID: continue
                # Everything immediately below is for
                # neighboring particles with a chainID. 
                if thisNN_chainID >= 0:
                    # Find thisNN's chain's max_dens.
                    thisNN_max_dens = self.densest_in_chain[thisNN_chainID]
                    if thisNN_max_dens == -1.0: continue
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
                        test = self.chain_densest_n[int(higher_chain)]
                    except KeyError:
                        self.chain_densest_n[int(higher_chain)] = {}
                    # See if this boundary density is higher than
                    # previously recorded for this pair of chains.
                    # Links only go one direction.
                    try:
                        old = self.chain_densest_n[int(higher_chain)][int(lower_chain)]
                        if old < boundary_density:
                            # make this the new densest boundary between this pair
                            self.chain_densest_n[int(higher_chain)][int(lower_chain)] = \
                                boundary_density
                    except KeyError:
                        # we haven't seen this pairing before, record this as the
                        # new densest boundary between chains
                        self.chain_densest_n[int(higher_chain)][int(lower_chain)] = \
                            boundary_density
                else:
                    continue
        try:
            del point, NNtags, results
        except UnboundLocalError:
            pass
        self.__max_memory()
        yt_counters("connect_chains")

    def _make_global_chain_densest_n(self):
        """
        We want to record the maximum boundary density between all chains on
        all tasks.
        """
        yt_counters("make_global_chain_densest_n")
        (self.top_keys, self.bot_keys, self.vals) = \
            self.linearize_chain_dict(self.chain_densest_n)
        self.__max_memory()
        del self.chain_densest_n
        yt_counters("make_global_chain_densest_n")

    def linearize_chain_dict(self, data):
        """
        Similar to above, but finds maximums for dicts of dicts. This is
        specificaly for a part of chainHOP.
        """
        top_keys = []
        bot_keys = []
        vals = []
        for top_key in data:
            for bot_key in data[top_key]:
                top_keys.append(top_key)
                bot_keys.append(bot_key)
                vals.append(data[top_key][bot_key])
        top_keys = np.array(top_keys, dtype='int64')
        bot_keys = np.array(bot_keys, dtype='int64')
        vals = np.array(vals, dtype='float64')

        data.clear()

        top_keys = self.comm.par_combine_object(top_keys, datatype='array', op='cat')
        bot_keys = self.comm.par_combine_object(bot_keys, datatype='array', op='cat')
        vals     = self.comm.par_combine_object(vals, datatype='array', op='cat')
        return (top_keys, bot_keys, vals)

    def _build_groups(self):
        """
        With the collection of possible chain links, build groups.
        """
        yt_counters("build_groups")
        # We need to find out which pairs of self.top_keys, self.bot_keys are
        # both < self.peakthresh, and create arrays that will store this
        # relationship.
        both = np.bitwise_and((self.densest_in_chain[self.top_keys] < self.peakthresh),
            (self.densest_in_chain[self.bot_keys] < self.peakthresh))
        g_high = self.top_keys[both]
        g_low = self.bot_keys[both]
        g_dens = self.vals[both]
        del both
        self.reverse_map = np.ones(self.densest_in_chain.size) * -1
        densestbound = np.ones(self.densest_in_chain.size) * -1.0
        for i, gl in enumerate(g_low):
            if g_dens[i] > densestbound[gl]:
                densestbound[gl] = g_dens[i]
        groupID = 0
        # First assign a group to all chains with max_dens above peakthresh.
        # The initial groupIDs will be assigned with decending peak density.
        # This guarantees that the group with the smaller groupID is the
        # higher chain, as in chain_high below.
        for chainID,density in enumerate(self.densest_in_chain):
            if density == -1.0: continue
            if self.densest_in_chain[chainID] >= self.peakthresh:
                self.reverse_map[chainID] = groupID
                groupID += 1
        group_equivalancy_map = np.empty(groupID, dtype='object')
        for i in xrange(groupID):
            group_equivalancy_map[i] = set([])
        # Loop over all of the chain linkages.
        for i,chain_high in enumerate(self.top_keys):
            chain_low = self.bot_keys[i]
            dens = self.vals[i]
            max_dens_high = self.densest_in_chain[chain_high]
            max_dens_low = self.densest_in_chain[chain_low]
            if max_dens_high == -1.0 or max_dens_low == -1.0: continue
            # If neither are peak density groups, mark them for later
            # consideration.
            if max_dens_high < self.peakthresh and \
                max_dens_low < self.peakthresh:
                    # This step is now done vectorized above, with the g_dens
                    # stuff.
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
                    if group_high == -1 or group_low == -1: continue
                    # Both are already identified as groups, so we need
                    # to re-assign the less dense group to the denser
                    # groupID.
                    if group_low != group_high:
                        group_equivalancy_map[group_low].add(group_high)
                        group_equivalancy_map[group_high].add(group_low)
                    continue
            # Else, one is above peakthresh, the other below
            # find out if this is the densest boundary seen so far for
            # the lower chain.
            group_high = self.reverse_map[chain_high]
            if group_high == -1: continue
            if dens >= densestbound[chain_low]:
                densestbound[chain_low] = dens
                self.reverse_map[chain_low] = group_high
        self.__max_memory()
        del self.top_keys, self.bot_keys, self.vals
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
        Set_list = []
        # We only want the holes that are modulo mine.
        keys = np.arange(groupID, dtype='int64')
        size = self.comm.size
        select = (keys % size == self.mine)
        groupIDs = keys[select]
        mine_groupIDs = set([]) # Records only ones modulo mine.
        not_mine_groupIDs = set([]) # All the others.
        # Declare these to prevent Errors when they're del-ed below, in case
        # this task doesn't create them in the loop, for whatever reason.
        current_sets, new_mine, new_other = [], [], []
        new_set, final_set, to_add_set, liter = set([]), set([]), set([]), set([])
        to_add_set = set([])
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
        self.__max_memory()
        del group_equivalancy_map, final_set, keys, select, groupIDs, current_sets
        del mine_groupIDs, not_mine_groupIDs, new_set, to_add_set, liter
        # Convert this list of sets into a look-up table
        lookup = np.ones(self.densest_in_chain.size, dtype='int64') * (self.densest_in_chain.size + 2)
        for i,item in enumerate(Set_list):
            item_min = min(item)
            for groupID in item:
                lookup[groupID] = item_min
        self.__max_memory()
        del Set_list
        # To bring it all together, find the minimum values at each entry
        # globally.
        lookup = self.comm.mpi_allreduce(lookup, op='min')
        # Now apply this to reverse_map
        for chainID,groupID in enumerate(self.reverse_map):
            if groupID == -1:
                continue
            if lookup[groupID] != (self.densest_in_chain.size + 2):
                self.reverse_map[chainID] = lookup[groupID]
        del lookup
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
                # replace it. We also don't want to link to un-assigned 
                # neighbors, and we can skip neighbors we're already assigned to.
                if dens >= densestbound[chain_low] and \
                        densestbound[chain_high] > densestbound[chain_low] and \
                        self.reverse_map[chain_high] != -1 and \
                        self.reverse_map[chain_low] != self.reverse_map[chain_high]:
                    changes += 1
                    if dens < densestbound[chain_high]:
                        densestbound[chain_low] = dens
                    else:
                        densestbound[chain_low] = densestbound[chain_high]
                    self.reverse_map[chain_low] = self.reverse_map[chain_high]
        self.__max_memory()
        del g_high, g_low, g_dens, densestbound
        # Now we have to find the unique groupIDs, since they may have been
        # merged.
        temp = list(set(self.reverse_map))
        # Remove -1 from the list.
        try:
            temp.pop(temp.index(-1))
        except ValueError:
            # There are no groups, probably.
            pass
        # Make a secondary map to make the IDs consecutive.
        values = np.arange(len(temp))
        secondary_map = dict(itertools.izip(temp, values))
        del values
        # Update reverse_map
        for chain, map in enumerate(self.reverse_map):
            # Don't attempt to fix non-assigned chains.
            if map == -1: continue
            self.reverse_map[chain] = secondary_map[map]
        group_count = len(temp)
        del secondary_map, temp
        yt_counters("build_groups")
        self.__max_memory()
        return group_count

    def _translate_groupIDs(self, group_count):
        """
        Using the maps, convert the particle chainIDs into their locally-final
        groupIDs.
        """
        yt_counters("translate_groupIDs")
        self.I_own = set([])
        for i in xrange(int(self.size)):
            # Don't translate non-affiliated particles.
            if self.chainID[i] == -1: continue
            # We want to remove the group tag from padded particles,
            # so when we return it to HaloFinding, there is no duplication.
            if self.is_inside[i]:
                self.chainID[i] = self.reverse_map[self.chainID[i]]
                self.I_own.add(self.chainID[i])
            else:
                self.chainID[i] = -1
        del self.is_inside
        # Create a densest_in_group, analogous to densest_in_chain.
        keys = np.arange(group_count)
        vals = np.zeros(group_count)
        self.densest_in_group = dict(itertools.izip(keys,vals))
        self.densest_in_group_real_index = self.densest_in_group.copy()
        del keys, vals
        for chainID,max_dens in enumerate(self.densest_in_chain):
            if max_dens == -1.0: continue
            groupID = self.reverse_map[chainID]
            if groupID == -1: continue
            if self.densest_in_group[groupID] < max_dens:
                self.densest_in_group[groupID] = max_dens
                self.densest_in_group_real_index[groupID] = self.densest_in_chain_real_index[chainID]
        del self.densest_in_chain, self.densest_in_chain_real_index, self.reverse_map
        del self.densest_in_group
        yt_counters("translate_groupIDs")

    def _precompute_group_info(self):
        yt_counters("Precomp.")
        """
        For all groups, compute the various global properties, except bulk
        velocity, to save time in HaloFinding.py (fewer barriers!).
        """
        select = (self.chainID != -1)
        calc = len(np.where(select == True)[0])
        loc = np.empty((calc, 3), dtype='float64')
        if self.tree == 'F':
            loc[:, 0] = np.concatenate((self.xpos, self.xpos_pad))[select]
            loc[:, 1] = np.concatenate((self.ypos, self.ypos_pad))[select]
            loc[:, 2] = np.concatenate((self.zpos, self.zpos_pad))[select]
            self.__max_memory()
            del self.xpos_pad, self.ypos_pad, self.zpos_pad
        elif self.tree == 'C':
            loc = self.pos[select]
        subchain = self.chainID[select]
        # First we need to find the maximum density point for all groups.
        # I think this will be faster than several vector operations that need
        # to pull the entire chainID array out of memory several times.
        yt_counters("max dens point")
        max_dens_point = np.zeros((self.group_count,4),dtype='float64')
        for i,part in enumerate(np.arange(self.size)[select]):
            groupID = self.chainID[part]
            if part < self.real_size:
                real_index = self.index[part]
            else:
                real_index = self.index_pad[part - self.real_size]
            if real_index == self.densest_in_group_real_index[groupID]:
                max_dens_point[groupID] = np.array([self.density[part], \
                loc[i, 0], loc[i, 1], loc[i, 2]])
        del self.index, self.index_pad, self.densest_in_group_real_index
        # Now we broadcast this, effectively, with an allsum. Even though
        # some groups are on multiple tasks, there is only one densest_in_chain
        # and only that task contributed above.
        self.max_dens_point = self.comm.mpi_allreduce(max_dens_point, op='sum')
        del max_dens_point
        yt_counters("max dens point")
        # Now CoM.
        yt_counters("CoM")
        CoM_M = np.zeros((self.group_count,3),dtype='float64')
        Tot_M = np.zeros(self.group_count, dtype='float64')
        #c_vec = self.max_dens_point[:,1:4][subchain] - np.array([0.5,0.5,0.5])
        if calc:
            c_vec = self.max_dens_point[:,1:4][subchain] - np.array([0.5,0.5,0.5])
            size = np.bincount(self.chainID[select]).astype('int64')
        else:
            # This task has no particles in groups!
            size = np.zeros(self.group_count, dtype='int64')
        # In case this task doesn't have all the groups, add trailing zeros.
        if size.size != self.group_count:
            size = np.concatenate((size, np.zeros(self.group_count - size.size, dtype='int64')))
        if calc:
            cc = loc - c_vec
            cc = cc - np.floor(cc)
            ms = np.concatenate((self.mass, self.mass_pad))[select]
            # Most of the time, the masses will be all the same, and we can try
            # to save some effort.
            ms_u = np.unique(ms)
            if ms_u.size == 1:
                single = True
                Tot_M = size.astype('float64') * ms_u
                del ms_u
            else:
                single = False
                del ms_u
            cc[:,0] = cc[:,0] * ms
            cc[:,1] = cc[:,1] * ms
            cc[:,2] = cc[:,2] * ms
            sort = subchain.argsort()
            cc = cc[sort]
            sort_subchain = subchain[sort]
            uniq_subchain = np.unique(sort_subchain)
            diff_subchain = np.ediff1d(sort_subchain)
            marks = (diff_subchain > 0)
            marks = np.arange(calc)[marks] + 1
            marks = np.concatenate(([0], marks, [calc]))
            for i, u in enumerate(uniq_subchain):
                CoM_M[u] = np.sum(cc[marks[i]:marks[i+1]], axis=0)
            if not single:
                for i,groupID in enumerate(subchain):
                    Tot_M[groupID] += ms[i]
            del cc, ms
            for groupID in xrange(int(self.group_count)):
                # Don't divide by zero.
                if groupID in self.I_own:
                    CoM_M[groupID] /= Tot_M[groupID]
                    CoM_M[groupID] += self.max_dens_point[groupID,1:4] - np.array([0.5,0.5,0.5])
                    CoM_M[groupID] *= Tot_M[groupID]
        # Now we find their global values
        self.group_sizes = self.comm.mpi_allreduce(size, op='sum')
        CoM_M = self.comm.mpi_allreduce(CoM_M, op='sum')
        self.Tot_M = self.comm.mpi_allreduce(Tot_M, op='sum')
        self.CoM = np.empty((self.group_count,3), dtype='float64')
        for groupID in xrange(int(self.group_count)):
            self.CoM[groupID] = CoM_M[groupID] / self.Tot_M[groupID]
        yt_counters("CoM")
        self.__max_memory()
        # Now we find the maximum radius for all groups.
        yt_counters("max radius")
        max_radius = np.zeros(self.group_count, dtype='float64')
        if calc:
            com = self.CoM[subchain]
            rad = np.fabs(com - loc)
            dist = (np.minimum(rad, self.period - rad)**2.).sum(axis=1)
            dist = dist[sort]
            for i, u in enumerate(uniq_subchain):
                max_radius[u] = np.max(dist[marks[i]:marks[i+1]])
        # Find the maximum across all tasks.
        mylog.info('Fraction of particles in this region in groups: %f' % (float(calc)/self.size))
        self.max_radius = self.comm.mpi_allreduce(max_radius, op='max')
        self.max_radius = np.sqrt(self.max_radius)
        yt_counters("max radius")
        yt_counters("Precomp.")
        self.__max_memory()
        del select, loc, subchain, CoM_M, Tot_M, size, max_radius
        if calc:
            del c_vec
            del sort_subchain, uniq_subchain, diff_subchain, marks, dist, sort
            del rad, com

    def _chain_hop(self):
        self._global_padding('first')
        self._global_bounds_neighbors()
        self._global_padding('second')
        self._is_inside('first')
        mylog.info('Distributing padded particles...')
        self._communicate_padding_data()
        mylog.info('Building kd tree for %d particles...' % \
            self.size)
        self._init_kd_tree()
        # Mark particles in as being in/out of the domain.
        self._is_inside('second')
        # Loop over the particles to find NN for each.
        mylog.info('Finding nearest neighbors/density...')
        yt_counters("chainHOP_tags_dens")
        if self.tree == 'F':
            chainHOP_tags_dens()
        elif self.tree == 'C':
            self.density = self.kdtree.chainHOP_get_dens(self.mass, \
            num_neighbors = self.num_neighbors, nMerge = self.nMerge + 2)
        yt_counters("chainHOP_tags_dens")
        if self.tree == 'F':
            self.density = fKD.dens.copy()
        elif self.tree == 'C':
            pass
        # Now each particle a local self density.
        # Let's find densest NN
        mylog.info('Finding densest nearest neighbors...')
        self._densestNN()
        # Build the chain of links.
        mylog.info('Building particle chains...')
        chain_count = self._build_chains()
        # This array tracks whether or not relationships for this particle
        # need to be examined twice, in preconnect_chains and in connect_chains
        self.search_again = np.ones(self.size, dtype='bool')
        if self.premerge:
            chain_count = self._preconnect_chains(chain_count)
        mylog.info('Gobally assigning chainIDs...')
        self._globally_assign_chainIDs(chain_count)
        mylog.info('Globally finding densest in chains...')
        self._create_global_densest_in_chain()
        mylog.info('Building chain connections across tasks...')
        self._is_inside('third')
        self._connect_chains_across_tasks()
        mylog.info('Communicating connected chains...')
        self._communicate_annulus_chainIDs()
        mylog.info('Connecting %d chains into groups...' % self.nchains)
        self._connect_chains()
        if self.tree == 'F':
            self.mass = fKD.mass[:self.psize]
            self.mass_pad = fKD.mass[self.psize:]
            del fKD.dens, fKD.mass, fKD.dens
            self.xpos = fKD.pos[0, :self.psize]
            self.ypos = fKD.pos[1, :self.psize]
            self.zpos = fKD.pos[2, :self.psize]
            self.xpos_pad = fKD.pos[0, self.psize:]
            self.ypos_pad = fKD.pos[1, self.psize:]
            self.zpos_pad = fKD.pos[2, self.psize:]
            del fKD.pos, fKD.chunk_tags
            free_tree(0) # Frees the kdtree object.
            gc.collect()
        elif self.tree == 'C':
            del self.kdtree
            gc.collect()
        del self.densestNN
        mylog.info('Communicating group links globally...')
        self._make_global_chain_densest_n()
        mylog.info('Building final groups...')
        group_count = self._build_groups()
        self.group_count = group_count
        mylog.info('Remapping particles to final groups...')
        self._translate_groupIDs(group_count)
        mylog.info('Precomputing info for %d groups...' % group_count)
        self._precompute_group_info()
        mylog.info("All done! Max Memory = %d MB" % self.max_mem)
        # We need to fix chainID and density because HaloFinding is expecting
        # an array only as long as the real data.
        self.chainID = self.chainID[:self.real_size]
        self.density = self.density[:self.real_size]
        # We'll make this a global object, which can be used to write a text
        # file giving the names of hdf5 files the particles for each halo.
        self.mine, self.I_own = self.comm.mpi_info_dict(self.I_own)
        self.halo_taskmap = defaultdict(set)
        for taskID in self.I_own:
            for groupID in self.I_own[taskID]:
                self.halo_taskmap[groupID].add(taskID)
        del self.I_own
        if self.tree == 'F':
            del self.xpos, self.ypos, self.zpos
        elif self.tree == 'C':
            pass

    def __add_to_array(self, arr, key, value, type):
        """
        In an effort to replace the functionality of a dict with an array, in
        order to save memory, this function adds items to an array. If the
        array is not long enough, it is resized and filled with 'bad' values."""
        
        try:
            arr[key] = value
        except IndexError:
            arr = np.concatenate((arr, np.ones(10000, dtype=type)*-1))
            arr[key] = value
        return arr
    
    def __clean_up_array(self, arr):
        good = (arr != -1)
        return arr[good]
    
    def __max_memory(self):
        my_mem = get_memory_usage()
        self.max_mem = max(my_mem, self.max_mem)
