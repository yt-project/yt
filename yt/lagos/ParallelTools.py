"""
Parallel data mapping techniques for yt

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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
from yt.funcs import *
import yt.logger
import itertools, sys

if os.path.basename(sys.argv[0]) == "mpi4py":
    from mpi4py import MPI
    parallel_capable = (MPI.COMM_WORLD.size > 1)
    if parallel_capable:
        mylog.info("Parallel computation enabled: %s / %s",
                   MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size)
        ytcfg["yt","__parallel_rank"] = str(MPI.COMM_WORLD.rank)
        ytcfg["yt","__parallel_size"] = str(MPI.COMM_WORLD.size)
        ytcfg["yt","__parallel"] = "True"
        # Now let's make sure we have the right options set.
        if MPI.COMM_WORLD.rank > 0:
            if ytcfg.getboolean("lagos","serialize"):
                ytcfg["lagos","onlydeserialize"] = "True"
            if ytcfg.getboolean("yt","LogFile"):
                ytcfg["yt","LogFile"] = "False"
                yt.logger.disable_file_logging()
else:
    parallel_capable = False

class GridIterator(object):
    def __init__(self, pobj, just_list = False):
        self.pobj = pobj
        if hasattr(pobj, '_grids') and pobj._grids is not None:
            gs = pobj._grids
        else:
            gs = pobj._data_source._grids
        if hasattr(gs[0], 'proc_num'):
            # This one sort of knows about MPI, but not quite
            self._grids = [g for g in gs if g.proc_num ==
                            ytcfg.getint('yt','__parallel_rank')]
            self._use_all = True
        else:
            self._grids = sorted(gs, key = lambda g: g.filename)
            self._use_all = False
        self.ng = len(self._grids)
        self.just_list = just_list

    def __iter__(self):
        self.pos = 0
        return self

    def next(self):
        # We do this manually in case
        # something else asks for us.pos
        if self.pos < len(self._grids):
            self.pos += 1
            return self._grids[self.pos - 1]
        raise StopIteration

class ParallelGridIterator(GridIterator):
    """
    This takes an object, pobj, that implements ParallelAnalysisInterface,
    and then does its thing.
    """
    def __init__(self, pobj, just_list = False):
        GridIterator.__init__(self, pobj, just_list)
        self._offset = MPI.COMM_WORLD.rank
        self._skip = MPI.COMM_WORLD.size
        # Note that we're doing this in advance, and with a simple means
        # of choosing them; more advanced methods will be explored later.
        if self._use_all:
            self.my_grid_ids = na.arange(len(self._grids))
        else:
            #upper, lower = na.mgrid[0:self.ng:(self._skip+1)*1j][self._offset:self._offset+2]
            #self.my_grid_ids = na.mgrid[upper:lower-1].astype("int64")
            self.my_grid_ids = na.array_split(
                            na.arange(len(self._grids)), self._skip)[self._offset]
        
    def __iter__(self):
        self.pos = 0
        return self

    def next(self):
        if self.pos < len(self.my_grid_ids):
            gid = self.my_grid_ids[self.pos]
            self.pos += 1
            return self._grids[gid]
        if not self.just_list: self.pobj._finalize_parallel()
        raise StopIteration

def parallel_passthrough(func):
    @wraps(func)
    def passage(self, data):
        if not parallel_capable: return data
        return func(self, data)
    return passage

class ParallelAnalysisInterface(object):
    _grids = None

    def _get_grids(self, *args, **kwargs):
        if parallel_capable:
            self._initialize_parallel(*args, **kwargs)
            return ParallelGridIterator(self)
        return GridIterator(self)

    def _get_grid_objs(self):
        if parallel_capable:
            return ParallelGridIterator(self, True)
        return GridIterator(self, True)

    def _initialize_parallel(self):
        pass

    def _finalize_parallel(self):
        pass

    def _partition_hierarchy_2d(self, axis):
        if not parallel_capable:
           return False, self.hierarchy.grid_collection(self.center, self.hierarchy.grids)

        cc = MPI.Compute_dims(MPI.COMM_WORLD.size, 2)
        mi = MPI.COMM_WORLD.rank
        cx, cy = na.unravel_index(mi, cc)
        x = na.mgrid[0:1:(cc[0]+1)*1j][cx:cx+2]
        y = na.mgrid[0:1:(cc[1]+1)*1j][cy:cy+2]

        LE = na.zeros(3, dtype='float64')
        RE = na.ones(3, dtype='float64')
        LE[x_dict[axis]] = x[0]  # It actually doesn't matter if this is x or y
        RE[x_dict[axis]] = x[1]
        LE[y_dict[axis]] = y[0]
        RE[y_dict[axis]] = y[1]

        return True, self.hierarchy.region(self.center, LE, RE)

    def _partition_hierarchy_3d(self, padding=0.0):
        if not parallel_capable:
           return False, self.hierarchy.grid_collection(self.center, self.hierarchy.grids)

        cc = MPI.Compute_dims(MPI.COMM_WORLD.size, 3)
        mi = MPI.COMM_WORLD.rank
        cx, cy, cz = na.unravel_index(mi, cc)
        x = na.mgrid[0:1:(cc[0]+1)*1j][cx:cx+2]
        y = na.mgrid[0:1:(cc[1]+1)*1j][cy:cy+2]
        z = na.mgrid[0:1:(cc[2]+1)*1j][cz:cz+2]

        LE = na.array([x[0], y[0], z[0]], dtype='float64')
        RE = na.array([x[1], y[1], z[1]], dtype='float64')

        if padding > 0:
            return True, \
                LE, RE, self.hierarchy.periodic_region(self.center, LE-padding, RE+padding)

        return LE, RE, self.hierarchy.region(self.center, LE, RE)
        

    @parallel_passthrough
    def _mpi_catdict(self, data):
        mylog.debug("Opening MPI Barrier on %s", MPI.COMM_WORLD.rank)
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.rank == 0:
            data = self.__mpi_recvdict(data)
        else:
            MPI.COMM_WORLD.Send(data, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.Bcast(data, root=0)
        MPI.COMM_WORLD.Barrier()
        return data

    @parallel_passthrough
    def __mpi_recvdict(self, data):
        # First we receive, then we make a new dict.
        for i in range(1,MPI.COMM_WORLD.size):
            buf = MPI.COMM_WORLD.Recv(source=i, tag=0)
            for j in buf: data[j] = na.concatenate([data[j],buf[j]], axis=-1)
        return data

    @parallel_passthrough
    def __mpi_recvlist(self, data):
        # First we receive, then we make a new list.
        for i in range(1,MPI.COMM_WORLD.size):
            buf = MPI.COMM_WORLD.Recv(source=i, tag=0)
            data += buf
        return na.array(data)

    @parallel_passthrough
    def _mpi_catlist(self, data):
        mylog.debug("Opening MPI Barrier on %s", MPI.COMM_WORLD.rank)
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.rank == 0:
            data = self.__mpi_recvlist(data)
        else:
            MPI.COMM_WORLD.Send(data, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.Bcast(data, root=0)
        MPI.COMM_WORLD.Barrier()
        return data

    @parallel_passthrough
    def __mpi_recvarray(self, data):
        # First we receive, then we make a new list.
        for i in range(1,MPI.COMM_WORLD.size):
            buf = MPI.COMM_WORLD.Recv(source=i, tag=0)
            data = na.concatenate([data, buf])
        return data

    @parallel_passthrough
    def _mpi_catarray(self, data):
        mylog.debug("Opening MPI Barrier on %s", MPI.COMM_WORLD.rank)
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.rank == 0:
            data = self.__mpi_recvarray(data)
        else:
            MPI.COMM_WORLD.Send(data, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.Bcast(data, root=0)
        MPI.COMM_WORLD.Barrier()

    def _should_i_write(self):
        if not parallel_capable: return True
        return (MPI.COMM_WORLD == 0)

    def _preload(self, grids, fields, queue):
        # This will preload if it detects we are parallel capable and
        # if so, we load *everything* that we need.  Use with some care.
        if not parallel_capable: return
        queue.preload(grids, fields)

    @parallel_passthrough
    def _mpi_allsum(self, data):
        MPI.COMM_WORLD.Barrier()
        return MPI.COMM_WORLD.Allreduce(data, op=MPI.SUM)


    def _get_dependencies(self, fields):
        deps = []
        for field in fields:
            deps += ensure_list(fieldInfo[field].get_dependencies().requested)
        return list(set(deps))
