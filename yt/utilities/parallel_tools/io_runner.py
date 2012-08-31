"""
A simple IO staging mechanism

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

import os
from .parallel_analysis_interface import ProcessorPool
from yt.utilities.io_handler import BaseIOHandler
from contextlib import contextmanager
import time

try:
    from .parallel_analysis_interface import MPI
except ImportError:
    pass

YT_TAG_MESSAGE = 317 # Cell 317 knows where to go

class IOCommunicator(BaseIOHandler):
    def __init__(self, pf, wg, pool):
        mylog.info("Initializing IOCommunicator")
        self.pf = pf
        self.wg = wg # We don't need to use this!
        self.pool = pool
        self.comm = pool.comm
        # We read our grids here
        self.grids = []
        storage = {}
        grids = pf.h.grids.tolist()
        grids.sort(key=lambda a:a.filename)
        for sto, g in parallel_objects(grids, storage = storage):
            sto.result = self.comm.rank
            sto.result_id = g.id
            self.grids.append(g)
        self._id_offset = pf.h.grids[0]._id_offset
        mylog.info("Reading from disk ...")
        self.initialize_data()
        mylog.info("Broadcasting ...")
        self.comm.comm.bcast(storage, root = wg.ranks[0])
        mylog.info("Done.")
        self.hooks = []

    def initialize_data(self):
        pf = self.pf
        fields = [f for f in pf.h.field_list
                  if not pf.field_info[f].particle_type]
        pfields = [f for f in pf.h.field_list
                   if pf.field_info[f].particle_type]
        # Preload is only defined for Enzo ...
        if pf.h.io._data_style == "enzo_packed_3d":
            self.queue = pf.h.io.queue
            pf.h.io.preload(self.grids, fields)
            for g in self.grids:
                for f in fields:
                    if f not in self.queue[g.id]:
                        d = np.zeros(g.ActiveDimensions, dtype='float64')
                        self.queue[g.id][f] = d
                for f in pfields:
                    self.queue[g.id][f] = self._read(g, f)
        else:
            self.queue = {}
            for g in self.grids:
                for f in fields + pfields:
                    self.queue[g.id][f] = pf.h.io._read(g, f)

    def _read(self, g, f):
        fi = self.pf.field_info[f]
        if fi.particle_type and g.NumberOfParticles == 0:
            # because this gets upcast to float
            return np.array([],dtype='float64')
        try:
            temp = self.pf.h.io._read_data_set(g, f)
        except:# self.pf.hierarchy.io._read_exception as exc:
            if fi.not_in_all:
                temp = np.zeros(g.ActiveDimensions, dtype='float64')
            else:
                raise
        return temp

    def wait(self):
        status = MPI.Status()
        while 1:
            if self.comm.comm.Iprobe(MPI.ANY_SOURCE,
                                YT_TAG_MESSAGE,
                                status = status):
                msg = self.comm.comm.recv(
                        source = status.source, tag = YT_TAG_MESSAGE)
                if msg['op'] == "end":
                    mylog.debug("Shutting down IO.")
                    break
                self._send_data(msg, status.source)
                status = MPI.Status()
            else:
                time.sleep(1e-2)

    def _send_data(self, msg, dest):
        grid_id = msg['grid_id']
        field = msg['field']
        ts = self.queue[grid_id][field].astype("float64")
        mylog.debug("Opening send to %s (%s)", dest, ts.shape)
        self.hooks.append(self.comm.comm.Isend([ts, MPI.DOUBLE], dest = dest))

class IOHandlerRemote(BaseIOHandler):
    _data_style = "remote"

    def __init__(self, pf, wg, pool):
        self.pf = pf
        self.wg = wg # probably won't need
        self.pool = pool
        self.comm = pool.comm
        self.proc_map = self.comm.comm.bcast(None,
                root = pool['io'].ranks[0])
        super(IOHandlerRemote, self).__init__()

    def _read_data_set(self, grid, field):
        dest = self.proc_map[grid.id]
        msg = dict(grid_id = grid.id, field = field, op="read")
        mylog.debug("Requesting %s for %s from %s", field, grid, dest)
        if self.pf.field_info[field].particle_type:
            data = np.empty(grid.NumberOfParticles, 'float64')
        else:
            data = np.empty(grid.ActiveDimensions, 'float64')
        hook = self.comm.comm.Irecv([data, MPI.DOUBLE], source = dest)
        self.comm.comm.send(msg, dest = dest, tag = YT_TAG_MESSAGE)
        mylog.debug("Waiting for data.")
        MPI.Request.Wait(hook)
        return data

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        #sl = tuple(reversed(sl))
        return self._read_data_set(grid,field)[sl]

    def terminate(self):
        msg = dict(op='end')
        if self.wg.comm.rank == 0:
            for rank in self.pool['io'].ranks:
                mylog.debug("Sending termination message to %s", rank)
                self.comm.comm.send(msg, dest=rank, tag=YT_TAG_MESSAGE)

@contextmanager
def remote_io(pf, wg, pool):
    original_io = pf.h.io
    pf.h.io = IOHandlerRemote(pf, wg, pool)
    yield
    pf.h.io.terminate()
    pf.h.io = original_io

def io_nodes(fn, n_io, n_work, func, *args, **kwargs):
    pool, wg = ProcessorPool.from_sizes([(n_io, "io"), (n_work, "work")])
    rv = None
    if wg.name == "work":
        pf = load(fn)
        with remote_io(pf, wg, pool):
            rv = func(pf, *args, **kwargs)
    elif wg.name == "io":
        pf = load(fn)
        io = IOCommunicator(pf, wg, pool)
        io.wait()
    # We should broadcast the result
    rv = pool.comm.mpi_bcast(rv, root=pool['work'].ranks[0])
    pool.free_all()
    mylog.debug("Return value: %s", rv)
    return rv

# Here is an example of how to use this functionality.
if __name__ == "__main__":
    def gq(pf):
        dd = pf.h.all_data()
        return dd.quantities["TotalQuantity"]("CellMassMsun")
    q = io_nodes("DD0087/DD0087", 8, 24, gq)
    mylog.info(q)


