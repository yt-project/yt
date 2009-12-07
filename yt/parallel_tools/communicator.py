"""
This module handles communicators and sub-communicators.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

from yt.parallel_tools import parallel_capable, mylog

if parallel_capable:
    from mpi4py import MPI

# This space reserved.  self.group should not be called directly, but
# instead a communicator group should be used.

__tocast = 'c'

def get_mpi_type(self, arr):
    # A dictionary won't do this, because
    # the dtypes won't hash to the same thing.
    dd = ((na.float, MPI.FLOAT),
          (na.int, MPI.INT),
          (na.double, MPI.DOUBLE),
          (na.long, MPI.LONG))
    for a,m in dd:
        if arr.dtype == a: return m
    return None

def __gen_mpi_op(opname):
    def pass_through(self, obj):
        return obj
    if not parallel_capable:
        return pass_through
    op = getattr(MPI, opname)
    def operation(self, obj):
        tt = get_mpi_type(obj)
        if tt is not None:
            return self.group.Allreduce(obj, op=op)
        return self.group.allreduce(obj, op=op)
    return operation

class ProcessCommunciator(object):

    def __init__(self):
        self.group = self.group

    def barrier(self):
        if not self._distributed: return
        mylog.debug("Opening MPI Barrier on %s", self.rank)
        self.group.Barrier()

    def send_array(arr, dest, tag = 0):
        if not isinstance(arr, na.ndarray):
            self.group.send((None,None), dest=dest, tag=tag)
            self.group.send(arr, dest=dest, tag=tag)
            return
        tmp = arr.view(__tocast) # Cast to CHAR
        # communicate type and shape
        self.group.send((arr.dtype.str, arr.shape), dest=dest, tag=tag)
        self.group.Send([arr, MPI.CHAR], dest=dest, tag=tag)
        del tmp

    def recv_array(source, tag = 0):
        dt, ne = self.group.recv(source=source, tag=tag)
        if dt is None and ne is None:
            return self.group.recv(source=source, tag=tag)
        arr = na.empty(ne, dtype=dt)
        tmp = arr.view(__tocast)
        self.group.Recv([tmp, MPI.CHAR], source=source, tag=tag)
        return arr

    def bcast_array(arr, root = 0):
        if self.rank == root:
            tmp = arr.view(__tocast) # Cast to CHAR
            self.group.bcast((arr.dtype.str, arr.shape), root=root)
        else:
            dt, ne = self.group.bcast(None, root=root)
            arr = na.empty(ne, dtype=dt)
            tmp = arr.view(__tocast)
        self.group.Bcast([tmp, MPI.CHAR], root=root)
        return arr

    def alltoallv_array(send, total_size, offsets, sizes):
        if len(send.shape) > 1:
            recv = []
            for i in range(send.shape[0]):
                recv.append(_alltoallv_array(send[i,:].copy(), total_size, offsets, sizes))
            recv = na.array(recv)
            return recv
        offset = offsets[self.rank]
        tmp_send = send.view(__tocast)
        recv = na.empty(total_size, dtype=send.dtype)
        recv[offset:offset+send.size] = send[:]
        dtr = send.dtype.itemsize / tmp_send.dtype.itemsize # > 1
        soff = [0] * self.size
        ssize = [tmp_send.size] * self.size
        roff = [off * dtr for off in offsets]
        rsize = [siz * dtr for siz in sizes]
        tmp_recv = recv.view(__tocast)
        self.group.Alltoallv((tmp_send, (ssize, soff), MPI.CHAR),
                                 (tmp_recv, (rsize, roff), MPI.CHAR))
        return recv

    def concat_dictionary(self, data):
        field_keys = data.keys()
        field_keys.sort()
        size = data[field_keys[0]].shape[-1]
        # We do an all-to-all, but the components end up getting replaced.
        # This is a python object; because it's an int that's being broadcast,
        # the pickling overhead is minimal.
        sizes = self.group.alltoall( [size]*self.size )
        offsets = na.add.accumulate([0] + sizes)[:-1]
        arr_size = self.group.allreduce(size, op=MPI.SUM)
        for key in field_keys:
            dd = data[key]
            rv = self.alltoallv_array(dd, arr_size, offsets, sizes)
            data[key] = rv
        return data

    def concat_list(self, data):
        # Very simple!  We just fill up an array with Nones and our data.
        to_send = [None] * self.size
        if data is not None:
            to_send[self.rank] = ensure_list(data)
        else:
            # If it's none, we just pickle an empty list so as not to break the
            # concatenation at the end
            to_send[self.rank] = []
        to_send = self.group.alltoall(to_send)
        to_return = []
        # Now we pop off the list from the front
        for p in xrange(self.size):
            to_return += to_send.pop(0)
        return to_return

    def concat_array(self, data):
        size = data.shape[-1]
        sizes = self.group.alltoall( [size]*self.size )
        offsets = na.add.accumulate([0] + sizes)[:-1]
        arr_size = self.group.allreduce(size, op=MPI.SUM)
        rv = self.alltoallv_array(data, arr_size, offsets, sizes)
        return rv

    min = __gen_mpi_op("MIN")
    max = __gen_mpi_op("MAX")
    sum = __gen_mpi_op("SUM")
    prod = __gen_mpi_op("PROD")

    @property
    def rank(self):
        return self.group.rank

    @property
    def size(self):
        return self.group.size
