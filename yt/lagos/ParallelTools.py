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
import itertools, sys, cStringIO

if os.path.basename(sys.executable) in ["mpi4py", "embed_enzo"] \
    or "--parallel" in sys.argv or '_parallel' in dir(sys):
    from mpi4py import MPI
    parallel_capable = (MPI.COMM_WORLD.size > 1)
    if parallel_capable:
        mylog.info("Parallel computation enabled: %s / %s",
                   MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size)
        ytcfg["yt","__parallel_rank"] = str(MPI.COMM_WORLD.rank)
        ytcfg["yt","__parallel_size"] = str(MPI.COMM_WORLD.size)
        ytcfg["yt","__parallel"] = "True"
        # I believe we do not need to turn this off manually
        #ytcfg["yt","StoreParameterFiles"] = "False"
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

def parallel_simple_proxy(func):
    if not parallel_capable: return func
    @wraps(func)
    def single_proc_results(self, *args, **kwargs):
        retval = None
        if self._processing or not self._distributed:
            return func(self, *args, **kwargs)
        if self._owner == MPI.COMM_WORLD.rank:
            self._processing = True
            retval = func(self, *args, **kwargs)
            self._processing = False
        retval = MPI.COMM_WORLD.Bcast(retval, root=self._owner)
        MPI.COMM_WORLD.Barrier()
        return retval
    return single_proc_results

class ParallelDummy(type):
    # All attributes that don't start with _ get replaced with
    # parallel_simple_proxy attributes.
    def __init__(cls, name, bases, d):
        super(ParallelDummy, cls).__init__(name, bases, d)
        skip = d.pop("dont_wrap", [])
        for attrname in d:
            if attrname.startswith("_") or attrname in skip: continue
            attr = getattr(cls, attrname)
            if type(attr) == types.MethodType:
                setattr(cls, attrname, parallel_simple_proxy(attr))

def parallel_passthrough(func):
    @wraps(func)
    def passage(self, data):
        if not parallel_capable: return data
        return func(self, data)
    return passage

def parallel_blocking_call(func):
    @wraps(func)
    def barrierize(*args, **kwargs):
        mylog.info("Entering barrier before %s", func.func_name)
        MPI.COMM_WORLD.Barrier()
        retval = func(*args, **kwargs)
        mylog.info("Entering barrier after %s", func.func_name)
        MPI.COMM_WORLD.Barrier()
        return retval
    if parallel_capable:
        return barrierize
    else:
        return func

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

        reg = self.hierarchy.region(self.center, LE, RE)
        return True, reg

    def _partition_hierarchy_3d(self, padding=0.0):
        LE, RE = self.pf["DomainLeftEdge"], self.pf["DomainRightEdge"]
        if not parallel_capable:
           return False, LE, RE, self.hierarchy.grid_collection(self.center, self.hierarchy.grids)

        cc = MPI.Compute_dims(MPI.COMM_WORLD.size, 3)
        mi = MPI.COMM_WORLD.rank
        cx, cy, cz = na.unravel_index(mi, cc)
        x = na.mgrid[LE[0]:RE[0]:(cc[0]+1)*1j][cx:cx+2]
        y = na.mgrid[LE[1]:RE[1]:(cc[1]+1)*1j][cy:cy+2]
        z = na.mgrid[LE[2]:RE[2]:(cc[2]+1)*1j][cz:cz+2]

        LE = na.array([x[0], y[0], z[0]], dtype='float64')
        RE = na.array([x[1], y[1], z[1]], dtype='float64')

        if padding > 0:
            return True, \
                LE, RE, self.hierarchy.periodic_region(self.center, LE-padding, RE+padding)

        return False, LE, RE, self.hierarchy.region(self.center, LE, RE)
        
    def _barrier(self):
        if not parallel_capable: return
        mylog.debug("Opening MPI Barrier on %s", MPI.COMM_WORLD.rank)
        MPI.COMM_WORLD.Barrier()

    @parallel_passthrough
    def _mpi_catdict(self, data):
        self._barrier()
        field_keys = data.keys()
        field_keys.sort()
        np = MPI.COMM_WORLD.size
        for key in field_keys:
            mylog.debug("Joining %s (%s) on %s", key, type(data[key]),
                        MPI.COMM_WORLD.rank)
            if MPI.COMM_WORLD.rank == 0:
                temp_data = []
                for i in range(1,np):
                    temp_data.append(_recv_array(source=i, tag=0))
                data[key] = na.concatenate([data[key]] + temp_data, axis=-1)
            else:
                _send_array(data[key], dest=0, tag=0)
            self._barrier()
            data[key] = _bcast_array(data[key])
        self._barrier()
        return data

    @parallel_passthrough
    def __mpi_recvlist(self, data):
        # First we receive, then we make a new list.
        data = ensure_list(data)
        for i in range(1,MPI.COMM_WORLD.size):
            buf = ensure_list(MPI.COMM_WORLD.recv(source=i, tag=0))
            data += buf
        return na.array(data)

    @parallel_passthrough
    def _mpi_catlist(self, data):
        self._barrier()
        if MPI.COMM_WORLD.rank == 0:
            data = self.__mpi_recvlist(data)
        else:
            MPI.COMM_WORLD.send(data, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.bcast(data, root=0)
        self._barrier()
        return data

    @parallel_passthrough
    def __mpi_recvarrays(self, data):
        # First we receive, then we make a new list.
        for i in range(1,MPI.COMM_WORLD.size):
            buf = _recv_array(source=i, tag=0)
            data = na.concatenate([data, buf])
        return data

    @parallel_passthrough
    def _mpi_catarray(self, data):
        self._barrier()
        if MPI.COMM_WORLD.rank == 0:
            data = self.__mpi_recvarrays(data)
        else:
            _send_array(data, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.Bcast(data, root=0)
        self._barrier()
        return data

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
        self._barrier()
        # We use old-school pickling here on the assumption the arrays are
        # relatively small ( < 1e7 elements )
        return MPI.COMM_WORLD.allreduce(data, op=MPI.SUM)

    def _mpi_info_dict(self, info):
        mylog.info("Parallel capable: %s", parallel_capable)
        if not parallel_capable: return 0, {0:info}
        self._barrier()
        data = None
        if MPI.COMM_WORLD.rank == 0:
            data = {0:info}
            for i in range(1, MPI.COMM_WORLD.size):
                data[i] = MPI.COMM_WORLD.recv(source=i, tag=0)
        else:
            MPI.COMM_WORLD.send(info, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", MPI.COMM_WORLD.rank)
        data = MPI.COMM_WORLD.bcast(data, root=0)
        self._barrier()
        return MPI.COMM_WORLD.rank, data

    def _get_dependencies(self, fields):
        deps = []
        fi = self.pf.field_info
        for field in fields:
            deps += ensure_list(fi[field].get_dependencies(pf=self.pf).requested)
        return list(set(deps))

    def _claim_object(self, obj):
        if not parallel_capable: return
        obj._owner = MPI.COMM_WORLD.rank
        obj._distributed = True

    def _write_on_root(self, fn):
        if not parallel_capable: return open(fn, "w")
        if MPI.COMM_WORLD.rank == 0:
            return open(fn, "w")
        else:
            return cStringIO.StringIO()

__tocast = 'c'

def _send_array(arr, dest, tag = 0):
    if not isinstance(arr, na.ndarray):
        MPI.COMM_WORLD.send((None,None), dest=dest, tag=tag)
        MPI.COMM_WORLD.send(arr, dest=dest, tag=tag)
        return
    tmp = arr.view(__tocast) # Cast to CHAR
    # communicate type and shape
    MPI.COMM_WORLD.send((arr.dtype.str, arr.shape), dest=dest, tag=tag)
    MPI.COMM_WORLD.Send([arr, MPI.CHAR], dest=dest, tag=tag)
    del tmp

def _recv_array(source, tag = 0):
    dt, ne = MPI.COMM_WORLD.recv(source=source, tag=tag)
    if dt is None and ne is None:
        return MPI.COMM_WORLD.recv(source=source, tag=tag)
    arr = na.empty(ne, dtype=dt)
    tmp = arr.view(__tocast)
    MPI.COMM_WORLD.Recv([tmp, MPI.CHAR], source=source, tag=tag)
    return arr

def _bcast_array(arr, root = 0):
    if MPI.COMM_WORLD.rank == root:
        tmp = arr.view(__tocast) # Cast to CHAR
        MPI.COMM_WORLD.bcast((arr.dtype.str, arr.shape), root=root)
    else:
        dt, ne = MPI.COMM_WORLD.bcast(None, root=root)
        arr = na.empty(ne, dtype=dt)
        tmp = arr.view(__tocast)
    MPI.COMM_WORLD.Bcast([tmp, MPI.CHAR], root=root)
    return arr
