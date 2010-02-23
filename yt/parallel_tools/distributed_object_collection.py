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

from yt.funcs import *
from yt.parallel_tools import ParallelAnalysisInterface
from itertools import izip

class DistributedObjectCollection(ParallelAnalysisInterface):
    valid = True

    def _get_object_info(self):
        pass

    def _set_object_info(self):
        pass

    def join_lists(self):
        info_dict = self._get_object_list_info()
        info_dict = self._mpi_catdict()
        self._set_object_info(info_dict)

    def _collect_objects(self, desired_indices):
        # We figure out which indices belong to which processor,
        # then we pack them up, and we send a list to each processor.
        requests = defaultdict(lambda: list)
        parents = self._object_parents[desired_indices]
        # Even if we have a million bricks, this should not take long.
        s = self._mpi_get_size()
        for i, p in izip(desired_indices, parents):
            requests[p].append(i)
        for p in sorted(requests):
            requests[p] = na.array(requests[p], dtype='int64')
            request_count.append(len(requests[p]))
        request_count = na.array(request_count, dtype='int64').reshape((s,1))
        # Now we distribute our requests to all the processors.
        # This is two-pass.  One to get the length of the arrays.  The second
        # pass is to get the actual indices themselves.
        requests = self._mpi_catdict(dict(request_count=request_count))
        # Now we have our final array of requests, with arrangement
        # (Nproc,Nproc).  First index corresponds to requesting proc, second to
        # sending.  So [them,us] = 5 means we owe 5, whereas [us, them] means
        # we are owed.

    def _pack_object(self, index):
        pass

    def _unpack_object(self, index):
        pass
