"""
A simple distributed object mechanism, for storing array-heavy objects.
Meant to be subclassed.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

from itertools import izip

import numpy as na

from yt.funcs import *

from .parallel_analysis_interface import ParallelAnalysisInterface

class DistributedObjectCollection(ParallelAnalysisInterface):
    valid = True

    def _get_object_info(self):
        pass

    def _set_object_info(self):
        pass

    def join_lists(self):
        info_dict = self._get_object_info()
        info_dict = self._mpi_catdict(info_dict)
        self._set_object_info(info_dict)

    def _collect_objects(self, desired_indices):
        # We figure out which indices belong to which processor,
        # then we pack them up, and we send a list to each processor.
        request_count = []
        owners = self._object_owners[desired_indices]
        mylog.debug("Owner list: %s", na.unique1d(owners))
        # Even if we have a million bricks, this should not take long.
        s = self._mpi_get_size()
        m = self._mpi_get_rank()
        requests = dict( ( (i, []) for i in xrange(s) ) )
        for i, p in izip(desired_indices, owners):
            requests[p].append(i)
        for p in sorted(requests):
            requests[p] = na.array(requests[p], dtype='int64')
            request_count.append(len(requests[p]))
        size = len(request_count)
        mylog.debug("Requesting: %s", request_count)
        request_count = na.array(request_count, dtype='int64')
        # Now we distribute our requests to all the processors.
        # This is two-pass.  One to get the length of the arrays.  The second
        # pass is to get the actual indices themselves.
        request_count = self._mpi_joindict({m : request_count})
        # Now we have our final array of requests, with arrangement
        # (Nproc,Nproc).  First index corresponds to requesting proc, second to
        # sending.  So [them,us] = 5 means we owe 5, whereas [us, them] means
        # we are owed.
        send_hooks = []
        dsend_buffers, dsend_hooks = [], []
        recv_hooks, recv_buffers = [], []
        drecv_buffers, drecv_hooks = [], []
        # We post our index-list and data receives from each processor.
        mylog.debug("Posting data buffer receives")
        proc_hooks = {}
        for p, request_from in request_count.items():
            if p == m: continue
            size = request_from[m]
            #if size == 0: continue
            # We post receives of the grids we *asked* for.
            # Note that indices into this are not necessarily processor ids.
            # So we store.  This has to go before the appends or it's an
            # off-by-one.
            mylog.debug("Setting up index buffer of size %s for receive from %s",
                        size, p)
            proc_hooks[len(drecv_buffers)] = p
            drecv_buffers.append(self._create_buffer(requests[p]))
            drecv_hooks.append(self._mpi_Irecv_double(drecv_buffers[-1], p, 1))
            recv_buffers.append(na.zeros(size, dtype='int64'))
            # Our index list goes on 0, our buffer goes on 1.  We know how big
            # the index list will be, now.
            recv_hooks.append(self._mpi_Irecv_long(recv_buffers[-1], p, 0))
        # Send our index lists into hte waiting buffers
        mylog.debug("Sending index lists")
        for p, ind_list in requests.items():
            if p == m: continue
            if len(ind_list) == 0: continue
            # Now, we actually send our index lists.
            send_hooks.append(self._mpi_Isend_long(ind_list, p, 0))
        # Now we post receives for all of the data buffers.
        mylog.debug("Sending data")
        for i in self._mpi_Request_Waititer(recv_hooks):
            # We get back the index, which here is identical to the processor
            # number doing the send.  At this point, we can post our receives.
            p = proc_hooks[i]
            mylog.debug("Processing from %s", p)
            ind_list = recv_buffers[i]
            dsend_buffers.append(self._create_buffer(ind_list))
            self._pack_buffer(ind_list, dsend_buffers[-1])
            dsend_hooks.append(self._mpi_Isend_double(
                dsend_buffers[-1], p, 1))
        mylog.debug("Waiting on data receives: %s", len(drecv_hooks))
        for i in self._mpi_Request_Waititer(drecv_hooks):
            mylog.debug("Unpacking from %s", proc_hooks[i])
            # Now we have to unpack our buffers
            # Our key into this is actually the request for the processor
            # number.
            p = proc_hooks[i]
            self._unpack_buffer(requests[p], drecv_buffers[i])
        mylog.debug("Finalizing sends: %s", len(dsend_hooks))
        for i in self._mpi_Request_Waititer(dsend_hooks):
            continue

    def _create_buffer(self, ind_list):
        pass

    def _pack_buffer(self, ind_list):
        pass

    def _unpack_buffer(self, ind_list, my_buffer):
        pass

