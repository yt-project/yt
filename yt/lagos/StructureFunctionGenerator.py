"""
Structure Function Generator

Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Stephen Skory.  All Rights Reserved.

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

from yt.lagos import *
from yt.math_utils import *
try:
    from yt.extensions.kdtree import *
except ImportError:
    raise RuntimeError("The Fortran kD-Tree did not import correctly.")

import math, sys, itertools, inspect, types, random
from collections import defaultdict

class StructFcnGen(ParallelAnalysisInterface):
    def __init__(self, pf, left_edge=None, right_edge=None,
            total_values=1000000, comm_size=10000, length_type="lin",
            length_number=10, length_range=None, vol_ratio = 1):
        """
        *total_values* (int) How many total (global) pair calculations to run
        for each of the structure functions specified.
        A single integer. Default: 1,000,000.
        *comm_size* (int) How entries are sent during communication.
        Default: 10,000.
        Set the parameters used to search the simulational volume with randomly
        placed 'rulers'.
        *length_type* (str) controls the even spacing of the rulers lengths in
        logarithmic or linear space, set by "log" or "lin", respectively.
        A single string. Default: "lin"
        *length_number* (int) sets how many lengths to create, evenly spaced by the above
        parameter.
        A single integer. Default: "10"
        *length_range* (float) a min/max pair for the range of values to search the over
        the simulational volume.
        A single pair (a list or array). Default: [0., 1/2*shortest box edge].
        """
        self.fsets = []
        self.fields = set([])
        # MPI stuff.
        self.size = self._mpi_get_size()
        self.mine = self._mpi_get_rank()
        self.vol_ratio = vol_ratio
        self.total_values = total_values / self.size
        # For communication.
        self.recv_hooks = []
        self.send_hooks = []
        self.done_hooks = []
        self.comm_size = comm_size
        self.pf = pf
        self.period = self.pf['DomainRightEdge'] - self.pf['DomainLeftEdge']
        self.min_edge = min(self.period)
        self.hierarchy = pf.h
        self.center = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0
        # Figure out the range of ruler lengths.
        if length_range == None:
            length_range = [0, self.min_edge/2.]
        else:
            if len(length_range) != 2:
                raise ValueError("length_range must have two values.")
            if length_range[1] <= length_range[0]:
                raise ValueError("length_range[1] must be larger than length_range[0]")
            if length_range[1] > self.min_edge/2.:
                length_range[1] = self.min_edge/2.
                mylog.info("Automatically adjusting length_range[1] to half the shortest box edge.")
        # Make the list of ruler lengths.
        if length_type == "lin":
            self.lengths = na.linspace(length_range[0], length_range[1],
                length_number)
        elif length_type == "log":
            self.lengths = na.logspace(math.log10(length_range[0]),
                math.log10(length_range[1]), length_number)
        else:
            # Something went wrong.
            raise SyntaxError("length_type is either \"lin\" or \"log\".")
        # Subdivide the volume.
        if not left_edge or not right_edge:
            self.left_edge = self.pf['DomainLeftEdge']
            self.right_edge = self.pf['DomainRightEdge']
            padded, self.LE, self.RE, self.ds = self._partition_hierarchy_3d(padding=0.,
                rank_ratio = self.vol_ratio)
        else:
            self.left_edge = left_edge
            self.right_edge = right_edge
            # We do this twice, first with no 'buffer' to get the unbuffered
            # self.LE/RE, and then second to get a buffered self.ds.
            padded, self.LE, self.RE, temp = \
                self._partition_region_3d(left_edge, right_edge,
                    rank_ratio=self.vol_ratio)
            padded, temp, temp, self.ds = \
                self._partition_region_3d(left_edge - self.lengths[-1], \
                right_edge + self.lengths[-1], rank_ratio=self.vol_ratio)
        mylog.info("LE %s RE %s %s" % (str(self.LE), str(self.RE), str(self.ds)))
        random.seed(-1 * self.mine)
    
    def add_function(self, function, fields, out_labels):
        """
        Add a structure function to the list that will be evaluated at the
        generated pairs of points.
        """
        fargs = inspect.getargspec(function)
        if len(fargs.args) != 4:
            raise SyntaxError("The structure function %s needs four arguments." %\
                function.__name__)
        fields = list(fields)
        if len(fields) < 1:
            raise SyntaxError("Please specify at least one field for function %s." %\
                function.__name__)
        out_labels = list(out_labels)
        if len(out_labels) < 1:
            raise SyntaxError("Please specify at least one out_labels for function %s." %\
                function.__name__)
        self.fsets.append(StructSet(self, function, self.min_edge, fields,
            out_labels))
        self.fields.update(fields)
        return self.fsets[-1]
    
    def run_generator(self):
        """
        After all the structure functions have been added, run the generator.
        """
        # We need a function!
        if len(self.fsets) == 0:
            mylog.error("You need to add at least one structure function!")
            return None
        # Do all the startup tasks.
        self._init_kd_tree()
        self._build_fields_vals()
        self._build_points_array()
        for bigloop, length in enumerate(self.lengths):
            #mylog.info("Doing length %1.5e" % length)
            # Things stop when this value below equals total_values.
            self.generated_points = 0
            self.gen_array = na.zeros(self.size, dtype='int64')
            self.comm_cycle_count = 0
            self.final_comm_cycle_count = 0
            self.sent_done = False
            self._setup_done_hooks()
            # While everyone else isn't done or I'm not done, we loop.
            # self.points[0][1] = -0.3 * (self.mine + 1) # for testing.
            uni = na.unique(self.points)
            #while (self.gen_array < self.total_values).any() or uni.size > 1:
            while self._should_cycle():
                self._setup_recv_arrays()
                self._send_arrays()
                self._mpi_Request_Waitall(self.send_hooks)
                self._mpi_Request_Waitall(self.recv_hooks)
                #self._barrier()
                # print self.mine, 'p', self.points
                #self._barrier()
                # print self.mine, 'r', self.recv_points
                #self._barrier()
                if (self.recv_points < -1.).any() or (self.recv_points > 1.).any(): # or \
                        #(na.abs(na.log10(na.abs(self.recv_points))) > 20).any():
                    raise ValueError("self.recv_points is no good!")
                self.points = self.recv_points.copy()
                self.fields_vals = self.recv_fields_vals.copy()
                self.gen_array = self.recv_gen_array.copy()
                self._eval_points(length)
                self.gen_array[self.mine] = self.generated_points
                uni = na.unique(self.points)
                #print 'yy', self.mine, uni.size, self.gen_array, self.comm_cycle_count
                self.comm_cycle_count += 1
                if self.generated_points == self.total_values:
                    self._send_done_toall()
            #print 'done!', self.mine
            #self._barrier()
            mylog.info("Length (%d of %d) %1.5e took %d communication cycles to complete." % \
                (bigloop+1, len(self.lengths), length, self.comm_cycle_count))
        self._allsum_bin_hits()
    
    def _init_kd_tree(self):
        """
        Builds the kd tree of grid center points.
        """
        # Grid cell centers.
        mylog.info("Building kD-Tree.")
        xp = self.ds["x"]
        yp = self.ds["y"]
        zp = self.ds["z"]
        fKD.pos = na.asfortranarray(na.empty((3,xp.size), dtype='float64'))
        fKD.pos[0, :] = xp[:]
        fKD.pos[1, :] = yp[:]
        fKD.pos[2, :] = zp[:]
        fKD.qv = na.empty(3, dtype='float64')
        fKD.nn = 1
        fKD.sort = False
        fKD.rearrange = True
        fKD.tags = na.empty(1, dtype='int64')
        fKD.dist = na.empty(1, dtype='float64')
        create_tree()

    def _build_fields_vals(self):
        """
        Builds an array to store the field values array.
        """
        self.fields_vals = na.empty((self.comm_size, len(self.fields)), \
            dtype='float64')
        # At the same time build a dict to label the columns.
        self.fields_columns = {}
        for i,field in enumerate(self.fields):
            self.fields_columns[field] = i

    def _build_points_array(self):
        """
        Initializes the array that contains the random points as all negatives
        to start with.
        """
        self.points = na.ones((self.comm_size, 6), dtype='float64') * -1.0

    def _setup_done_hooks(self):
        """
        Opens non-blocking receives pointing to all the other tasks.
        """
        self.recv_done = {}
        for task in xrange(self.size):
            if task == self.mine: continue
            self.recv_done[task] = na.zeros(1, dtype='int64')
            self.done_hooks.append(self._mpi_Irecv_long(self.recv_done[task], \
                task, tag=15))

    def _send_done_toall(self):
        """
        Signal all the other tasks that this task has created all the points
        it needs to.
        """
        # If I've already done this, don't do it again.
        if self.sent_done: return
        self.send_done = {}
        for task in xrange(self.size):
            if task == self.mine: continue
            # We'll send the cycle at which I think things should stop.
            self.send_done[task] = na.ones(1, dtype='int64') * \
                (self.size / self.vol_ratio -1) + self.comm_cycle_count
            self.done_hooks.append(self._mpi_Isend_long(self.send_done[task], \
                task, tag=15))
        # I need to mark my *own*, too!
        self.recv_done[self.mine] = (self.size / self.vol_ratio -1) + \
            self.comm_cycle_count
        self.sent_done = True

    def _should_cycle(self):
        """
        Depending on several factors, it returns whether or not the
        communication cycle should continue.
        """
        # We need to continue if:
        # We haven't seen that all the points have been created.
        #if (self.gen_array < self.total_values).any(): return True
        # We own unprocessed points.
        #if uni.size > 1: return True
        # If I'm not finished.
        if self.generated_points < self.total_values: return True
        # If other tasks aren't finished
        if not self._mpi_Request_Testall(self.done_hooks): return True
        # If they are all finished, meaning Testall returns True, we find
        # the biggest value in self.recv_done and stop there.
        #if self.final_comm_cycle_count==0: # For testing.
        #    print self.mine, "its cycle",self.comm_cycle_count,"and I think everyone is done",self.recv_done
        stop = max(self.recv_done.values())
        if self.comm_cycle_count < stop:
            self.final_comm_cycle_count += 1
            return True
        # If we've come this far, we're done.
        return False

    def _setup_recv_arrays(self):
        """
        Creates the recv buffers and calls a non-blocking MPI receive pointing
        to the left-hand neighbor.
        """
        self.recv_points = na.ones((self.comm_size, 6), dtype='float64') * -1.
        self.recv_fields_vals = na.zeros((self.comm_size, len(self.fields)), \
            dtype='float64')
        self.recv_gen_array = na.zeros(self.size, dtype='int64')
        self.recv_hooks.append(self._mpi_Irecv_double(self.recv_points, \
            (self.mine-1)%self.size, tag=10))
        self.recv_hooks.append(self._mpi_Irecv_double(self.recv_fields_vals, \
            (self.mine-1)%self.size, tag=20))
        self.recv_hooks.append(self._mpi_Irecv_long(self.recv_gen_array, \
            (self.mine-1)%self.size, tag=40))

    def _send_arrays(self):
        """
        Send the data arrays to the right-hand neighbor.
        """
        self.send_hooks.append(self._mpi_Isend_double(self.points,\
            (self.mine+1)%self.size, tag=10))
        self.send_hooks.append(self._mpi_Isend_double(self.fields_vals,\
            (self.mine+1)%self.size, tag=20))
        self.send_hooks.append(self._mpi_Isend_long(self.gen_array, \
            (self.mine+1)%self.size, tag=40))

    def _allsum_bin_hits(self):
        """
        Add up the hits to all the bins globally for all functions.
        """
        print 'here'
        for fset in self.fsets:
            fset.too_low = self._mpi_allsum(fset.too_low)
            fset.too_high = self._mpi_allsum(fset.too_high)
            fset.binned = {}
            if self.mine == 0:
                mylog.info("Function %s had %d values too high and %d too low that were not binned." % \
                    (fset.function.__name__, fset.too_high, fset.too_low))
            for length in self.lengths:
                fset.length_bin_hits[length] = \
                    self._mpi_Allsum_long(fset.length_bin_hits[length])
                # Find out how many were successfully binned.
                fset.binned[length] = fset.length_bin_hits[length].sum()
                # Normalize the counts.
                fset.length_bin_hits[length] = \
                    fset.length_bin_hits[length].astype('float64') / \
                    fset.binned[length]
                # Return it to its original shape.
                fset.length_bin_hits[length] = \
                    fset.length_bin_hits[length].reshape(fset.bin_number)
    
    def _pick_random_points(self, length):
        """
        Picks out two random points separated by *length*.
        """
        # A random point inside this subvolume.
        r1 = na.empty(3, dtype='float64')
        for i in xrange(3):
            r1[i] = random.uniform(self.LE[i], self.RE[i]) # 3 randoms.
        # Pick our theta, phi random pair.
        theta = random.uniform(0, 2.*math.pi) # 1 random.
        phi = random.uniform(-math.pi/2., math.pi/2.) # 1 random.
        # Find our second point.
        r2 = na.array([r1[0] + length * math.cos(theta) * math.cos(phi),
            r1[1] + length * math.sin(theta) * math.cos(phi),
            r1[2] + length * math.sin(phi)], dtype='float64')
        # Reflect it so it's inside the (full) volume.
        r2 = r2 % self.period
        return (r1, na.array(r2))

    def _find_nearest_cell(self, r):
        """
        Perform the search for the closest 'grid' point.
        """
        fKD.qv = r
        find_nn_nearest_neighbors()
        # The -1 is because the kdtree is using fortran ordering which starts
        # at 1
        n = fKD.tags[0] - 1
        return n

    def _get_fields_vals(self, r):
        """
        Given a point r, return the values for the fields we need for that
        point.
        """
        # First find the grid data index field.
        index = self._find_nearest_cell(r)
        results = {}
        for field in self.fields:
            results[field] = self.ds[field][index]
        return results

    def _eval_points(self, length):
        """
        The main work loop. This loops over the self.points array. It attempts
        to evaluate the structure functions at the points, and if it can't,
        moves on. If it comes upon an empty (negative entries) entry, it 
        creates a new pair of random points, and attempts to evaluate them.
        If this pair are both local, the structure function results are binned
        and the loop over self.points is not advanced. If it cannot, we do
        advance. Things stop when we reach the end of the self.points array, 
        meaning it's full of points this task cannot evaluate. Things also
        stop when this task reaches self.total_values, and it doesn't have to
        make any more points.
        """
        #print 'starting',self.mine, self.generated_points
        broken = 0
        evaled = 0
        for i, points in enumerate(self.points):
            if (points < 0).any():
                # Make new points for this slot and don't advance the points
                # array as long as they are local.
                while self.generated_points < self.total_values:
                    r1, r2 = self._pick_random_points(length)
                    self.generated_points += 1
                    # Get the values for r1 because it's always local.
                    r1_results = self._get_fields_vals(r1)
                    # Store them.
                    for field in r1_results:
                        col = self.fields_columns[field]
                        self.fields_vals[i, col] = r1_results[field]
                    # if r2 is not local, we can skip out of this while loop
                    # now.
                    if (r2 < self.ds.left_edge).any() or (r2 >= self.ds.right_edge).any():
                        broken += 1
                        self.points[i][:3] = r1
                        self.points[i][3:] = r2
                        break
                    r2_results = self._get_fields_vals(r2)
                    # Evaluate the functions and bin the results.
                    for fcn_set in self.fsets:
                        fcn_results = fcn_set._eval_st_fcn(self.fields_vals[i],
                            r2_results, r1, r2)
                        fcn_set._bin_results(length, fcn_results)
                        evaled += 1
            else:
                # Now we've come to a point I've received and we'll attempt
                # to evaluate it.
                if (points[3:] < self.ds.left_edge).any() or (points[3:] >= self.ds.right_edge).any():
                    # I don't own this r2 point, move on.
                    #print 'moving on', points, self.LE, self.RE
                    continue
                r2 = points[3:]
                r2_results = self._get_fields_vals(r2)
                for fcn_set in self.fsets:
                    fcn_results = fcn_set._eval_st_fcn(self.fields_vals[i],
                        r2_results, points[:3], r2)
                    fcn_set._bin_results(length, fcn_results)
                    evaled += 1
                # reset points
                self.points[i] = na.array([-1.]*6)
                # If we need to, we add our own points here as long as we can.
                while self.generated_points < self.total_values:
                    r1, r2 = self._pick_random_points(length)
                    #self.points[i][:3] = r1
                    #self.points[i][3:] = r2
                    self.generated_points += 1
                    # Get the values for r1 because it's always local.
                    r1_results = self._get_fields_vals(r1)
                    # Store them.
                    for field in r1_results:
                        col = self.fields_columns[field]
                        self.fields_vals[i, col] = r1_results[field]
                    # if r2 is not local, we can skip out of this while loop
                    # now.
                    if (r2 < self.ds.left_edge).any() or (r2 >= self.ds.left_edge).any():
                        broken += 1
                        self.points[i][:3] = r1
                        self.points[i][3:] = r2
                        break
                    r2_results = self._get_fields_vals(r2)
                    # Evaluate the functions and bin the results.
                    for fcn_set in self.fsets:
                        fcn_results = fcn_set._eval_st_fcn(self.fields_vals[i],
                            r2_results, r1, r2)
                        fcn_set._bin_results(length, fcn_results)
                        evaled += 1
        #print 'done here',self.mine, self.generated_points, broken, evaled

    @parallel_blocking_call
    def write_out_means(self, sqrt=False):
        """
        Writes out the weighted-average value for each structure function for
        each dimension for each ruler length to a text file. The data is written
        to files of the name 'function_name.txt' in the current working
        directory.
        *sqrt* (bool) Sets whether or not to square-root the averages. For
        example, for RMS velocity differences set this to True, otherwise the
        output will be just MS velocity. Default: False.
        """
        sep = 12
        if type(sqrt) != types.ListType:
            sqrt = [sqrt]
        for fset in self.fsets:
            fp = self._write_on_root("%s.txt" % fset.function.__name__)
            fset._avg_bin_hits()
            line = "# length".ljust(sep)
            line += "count".ljust(sep)
            for dim in fset.dims:
                line += ("%s" % fset.out_labels[dim]).ljust(sep)
            fp.write(line + "\n")
            for length in self.lengths:
                line = ("%1.5e" % length).ljust(sep)
                line += ("%d" % fset.binned[length]).ljust(sep)
                for dim in fset.dims:
                    if sqrt[dim]:
                        line += ("%1.5e" % \
                            math.sqrt(fset.length_avgs[length][dim])).ljust(sep)
                    else:
                        line += ("%1.5e" % \
                            fset.length_avgs[length][dim]).ljust(sep)
                line += "\n"
                fp.write(line)
            fp.close()
    
    @parallel_root_only
    def write_out_arrays(self):
        """
        Writes out the raw probability bins and the bin edges to an HDF5 file
        for each of the structure functions. The files are named 
        'function_name.txt' and saved in the current working directory.
        """
        if self.mine == 0:
            for fset in self.fsets:
                f = h5py.File("%s.h5" % fset.function.__name__, "w")
                bin_names = []
                prob_names = []
                bin_counts = []
                for dim in fset.dims:
                    f.create_dataset("/bin_edges_%02d_%s" % \
                        (dim, fset.out_labels[dim]), \
                        data=fset.bin_edges[dim])
                    bin_names.append("/bin_edges_%02d_%s" % \
                    (dim, fset.out_labels[dim]))
                for i,length in enumerate(self.lengths):
                    f.create_dataset("/prob_bins_%05d" % i, \
                        data=fset.length_bin_hits[length])
                    prob_names.append("/prob_bins_%05d" % i)
                    bin_counts.append(fset.binned[length])
                f.create_dataset("/bin_edges_names", data=bin_names)
                #f.create_dataset("/prob_names", data=prob_names)
                f.create_dataset("/lengths", data=self.lengths)
                f.create_dataset("/counts", data=bin_counts)
                f.close()

class StructSet(StructFcnGen):
    def __init__(self,sfg, function, min_edge, fields, out_labels):
        self.sfg = sfg # The overarching SFG class
        self.function = function # Function to eval between the two points.
        self.min_edge = min_edge # The length of the minimum edge of the box.
        self.function_fields = fields # The fields for this fcn in order.
        self.out_labels = out_labels # For output.
        # These below are used to track how many times the function returns
        # unbinned results.
        self.too_low = 0
        self.too_high = 0
        
    def set_pdf_params(self, bin_type="lin", bin_number=1000, bin_range=None):
        """
        Set the parameters used to build the Probability Distribution Function
        for each ruler length for this function. The values output by the
        function are slotted into the bins described here.
        *bin_type* (str) controls the edges of the bins spaced evenly in
        logarithmic or linear space, set by "log" or "lin", respectively.
        A single string, or list of strings for N-dim binning. Default: "lin".
        *bin_number* (int) sets how many bins to create, evenly spaced by the above
        parameter.
        A single integer, or a list of integers for N-dim binning. Default: 1000
        *bin_range* (float) A pair of values giving the range for the bins.
        A pair of floats (a list), or a list of pairs for N-dim binning.
        Default: None.
        """
        # This should be called after setSearchParams.
        if not hasattr(self.sfg, "lengths"):
            mylog.error("Please call setSearchParams() before calling setPDFParams().")
            return None
        # Make sure they're either all lists or only one is.
        input = [bin_type, bin_number, bin_range]
        lists = 0
        for thing in input:
            if type(thing) == types.ListType:
                lists += 1
        if lists > 1 and lists < 3:
            mylog.error("Either all the inputs need to be lists, or only one.")
            return None
        # Make sure they're all the same length if they're lists.
        if lists == 3:
            first_len = 0
            for thing in input:
                if first_len == 0:
                    first_len = len(thing)
                    if first_len == 0:
                        mylog.error("Input cannot be an empty list.")
                        return None
                    continue
                if first_len != len(thing):
                    mylog.error("All the inputs need to have the same length.")
                    return None
        # If they are not all lists, put the input into lists for convenience.
        if lists == 1:
            bin_type, bin_number = [bin_type], [bin_number]
            bin_range = [bin_range]
        self.bin_type = bin_type
        self.bin_number = na.array(bin_number) - 1
        self.dims = range(len(bin_type))
        # Create the dict that stores the arrays to store the bin hits, and
        # the arrays themselves.
        self.length_bin_hits = {}
        for length in self.sfg.lengths:
            # It's easier to index flattened, but will be unflattened later.
            self.length_bin_hits[length] = na.zeros(self.bin_number,
                dtype='int64').flatten()
        # Create the bin edges for each dimension.
        # self.bins is indexed by dimension
        self.bin_edges = {}
        for dim in self.dims:
            # Error check.
            if len(bin_range[dim]) != 2:
                raise ValueError("bin_range must have two values.")
            if bin_range[dim][1] <= bin_range[dim][0]:
                raise ValueError("bin_range[1] must be larger than bin_range[0]")
            # Make the edges for this dimension.
            if bin_type[dim] == "lin":
                self.bin_edges[dim] = na.linspace(bin_range[dim][0], bin_range[dim][1],
                    bin_number[dim])
            elif bin_type[dim] == "log":
                self.bin_edges[dim] = na.logspace(math.log10(bin_range[dim][0]),
                    math.log10(bin_range[dim][1]), bin_number[dim])
            else:
                raise SyntaxError("bin_edges is either \"lin\" or \"log\".")

    def _eval_st_fcn(self, r1_results, r2_results, r1, r2):
        """
        Return the value of the structure function using the provided results.
        """
        # We need to pick out the field values this function needs and put them in
        # correct order. r1_results is an array, and r2_results is a dict,
        # which although it's a bit confusing, is actually the simplest way
        # to go.
        fcn_r1_res, fcn_r2_res = [], []
        for field in self.function_fields:
            fcn_r1_res.append(r1_results[self.sfg.fields_columns[field]])
            fcn_r2_res.append(r2_results[field])
        return self.function(fcn_r1_res, fcn_r2_res, r1, r2)
        """
        NOTE - A function looks like:
        def stuff(a,b):
            return [(a[0] - b[0])/(a[1] + b[1])]
        where a and b refer to different points in space and the indices
        are for the different fields, which are given when the function is
        added. The results need to be a list or array even if it's only one
        item.
        """

    def _bin_results(self, length, results):
        """
        Add hits to the bins corresponding to these results. length_hit_bins
        is flattened, so we need to figure out the offset for this hit by
        factoring the sizes of the other dimensions.
        """
        hit_bin = 0
        multi = 1
        for dim, res in enumerate(results):
            if res < self.bin_edges[dim][0]:
                self.too_low += 1
                return
            if res >= self.bin_edges[dim][-1]:
                self.too_high += 1
                return
            for d1 in range(dim):
                multi *= self.bin_edges[d1].size
            hit_bin += (na.digitize([res], self.bin_edges[dim])[0] - 1) * multi
            if hit_bin >= self.length_bin_hits[length].size:
                self.too_high += 1
                return
        self.length_bin_hits[length][hit_bin] += 1

    def _dim_sum(self, a, dim):
        """
        Given a multidimensional array a, this finds the sum over all the
        elements leaving the dimension dim untouched.
        """
        dims = na.arange(len(a.shape))
        dims = na.flipud(dims)
        gt_dims = dims[dims > dim]
        lt_dims = dims[dims < dim]
        iter_dims = na.concatenate((gt_dims, lt_dims))
        for this_dim in iter_dims:
            a = a.sum(axis=this_dim)
        return a

    def _avg_bin_hits(self):
        """
        For each dimension and length of bin_hits return the weighted average.
        """
        self.length_avgs = defaultdict(dict)
        for length in self.sfg.lengths:
            for dim in self.dims:
                self.length_avgs[length][dim] = \
                    (self._dim_sum(self.length_bin_hits[length], dim) * \
                    ((self.bin_edges[dim][:-1] + self.bin_edges[dim][1:]) / 2.)).sum()
