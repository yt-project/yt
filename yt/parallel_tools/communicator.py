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

from yt.parallel_tools import parallel_capable

if parallel_capable:
    from mpi4py import MPI

# This space reserved.  self._comm should not be called directly, but
# instead a communicator group should be used.

__tocast = 'c'

class ProcessCommunciator(object):
    def __init__(self):
        self._comm = MPI.COMM_WORLD

    def send_array(arr, dest, tag = 0):
        if not isinstance(arr, na.ndarray):
            self._comm.send((None,None), dest=dest, tag=tag)
            self._comm.send(arr, dest=dest, tag=tag)
            return
        tmp = arr.view(__tocast) # Cast to CHAR
        # communicate type and shape
        self._comm.send((arr.dtype.str, arr.shape), dest=dest, tag=tag)
        self._comm.Send([arr, MPI.CHAR], dest=dest, tag=tag)
        del tmp

    def recv_array(source, tag = 0):
        dt, ne = self._comm.recv(source=source, tag=tag)
        if dt is None and ne is None:
            return self._comm.recv(source=source, tag=tag)
        arr = na.empty(ne, dtype=dt)
        tmp = arr.view(__tocast)
        self._comm.Recv([tmp, MPI.CHAR], source=source, tag=tag)
        return arr

    def bcast_array(arr, root = 0):
        if self._comm.rank == root:
            tmp = arr.view(__tocast) # Cast to CHAR
            self._comm.bcast((arr.dtype.str, arr.shape), root=root)
        else:
            dt, ne = self._comm.bcast(None, root=root)
            arr = na.empty(ne, dtype=dt)
            tmp = arr.view(__tocast)
        self._comm.Bcast([tmp, MPI.CHAR], root=root)
        return arr

    def alltoallv_array(send, total_size, offsets, sizes):
        if len(send.shape) > 1:
            recv = []
            for i in range(send.shape[0]):
                recv.append(_alltoallv_array(send[i,:].copy(), total_size, offsets, sizes))
            recv = na.array(recv)
            return recv
        offset = offsets[self._comm.rank]
        tmp_send = send.view(__tocast)
        recv = na.empty(total_size, dtype=send.dtype)
        recv[offset:offset+send.size] = send[:]
        dtr = send.dtype.itemsize / tmp_send.dtype.itemsize # > 1
        soff = [0] * self._comm.size
        ssize = [tmp_send.size] * self._comm.size
        roff = [off * dtr for off in offsets]
        rsize = [siz * dtr for siz in sizes]
        tmp_recv = recv.view(__tocast)
        self._comm.Alltoallv((tmp_send, (ssize, soff), MPI.CHAR),
                                 (tmp_recv, (rsize, roff), MPI.CHAR))
        return recv
        
