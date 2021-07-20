import time
from contextlib import contextmanager

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

from .parallel_analysis_interface import ProcessorPool, parallel_objects

try:
    from .parallel_analysis_interface import MPI
except ImportError:
    pass

YT_TAG_MESSAGE = 317  # Cell 317 knows where to go


class IOCommunicator(BaseIOHandler):
    def __init__(self, ds, wg, pool):
        mylog.info("Initializing IOCommunicator")
        self.ds = ds
        self.wg = wg  # We don't need to use this!
        self.pool = pool
        self.comm = pool.comm
        # We read our grids here
        self.grids = []
        storage = {}
        grids = ds.index.grids.tolist()
        grids.sort(key=lambda a: a.filename)
        for sto, g in parallel_objects(grids, storage=storage):
            sto.result = self.comm.rank
            sto.result_id = g.id
            self.grids.append(g)
        self._id_offset = ds.index.grids[0]._id_offset
        mylog.info("Reading from disk ...")
        self.initialize_data()
        mylog.info("Broadcasting ...")
        self.comm.comm.bcast(storage, root=wg.ranks[0])
        mylog.info("Done.")
        self.hooks = []

    def initialize_data(self):
        ds = self.ds
        fields = [
            f for f in ds.field_list if not ds.field_info[f].sampling_type == "particle"
        ]
        pfields = [
            f for f in ds.field_list if ds.field_info[f].sampling_type == "particle"
        ]
        # Preload is only defined for Enzo ...
        if ds.index.io._dataset_type == "enzo_packed_3d":
            self.queue = ds.index.io.queue
            ds.index.io.preload(self.grids, fields)
            for g in self.grids:
                for f in fields:
                    if f not in self.queue[g.id]:
                        d = np.zeros(g.ActiveDimensions, dtype="float64")
                        self.queue[g.id][f] = d
                for f in pfields:
                    self.queue[g.id][f] = self._read(g, f)
        else:
            self.queue = {}
            for g in self.grids:
                for f in fields + pfields:
                    self.queue[g.id][f] = ds.index.io._read(g, f)

    def _read(self, g, f):
        fi = self.ds.field_info[f]
        if fi.sampling_type == "particle" and g.NumberOfParticles == 0:
            # because this gets upcast to float
            return np.array([], dtype="float64")
        try:
            temp = self.ds.index.io._read_data_set(g, f)
        except Exception:  # self.ds.index.io._read_exception as exc:
            if fi.not_in_all:
                temp = np.zeros(g.ActiveDimensions, dtype="float64")
            else:
                raise
        return temp

    def wait(self):
        status = MPI.Status()
        while True:
            if self.comm.comm.Iprobe(MPI.ANY_SOURCE, YT_TAG_MESSAGE, status=status):
                msg = self.comm.comm.recv(source=status.source, tag=YT_TAG_MESSAGE)
                if msg["op"] == "end":
                    mylog.debug("Shutting down IO.")
                    break
                self._send_data(msg, status.source)
                status = MPI.Status()
            else:
                time.sleep(1e-2)

    def _send_data(self, msg, dest):
        grid_id = msg["grid_id"]
        field = msg["field"]
        ts = self.queue[grid_id][field].astype("float64")
        mylog.debug("Opening send to %s (%s)", dest, ts.shape)
        self.hooks.append(self.comm.comm.Isend([ts, MPI.DOUBLE], dest=dest))


class IOHandlerRemote(BaseIOHandler):
    _dataset_type = "remote"

    def __init__(self, ds, wg, pool):
        self.ds = ds
        self.wg = wg  # probably won't need
        self.pool = pool
        self.comm = pool.comm
        self.proc_map = self.comm.comm.bcast(None, root=pool["io"].ranks[0])
        super().__init__()

    def _read_data_set(self, grid, field):
        dest = self.proc_map[grid.id]
        msg = dict(grid_id=grid.id, field=field, op="read")
        mylog.debug("Requesting %s for %s from %s", field, grid, dest)
        if self.ds.field_info[field].sampling_type == "particle":
            data = np.empty(grid.NumberOfParticles, "float64")
        else:
            data = np.empty(grid.ActiveDimensions, "float64")
        hook = self.comm.comm.Irecv([data, MPI.DOUBLE], source=dest)
        self.comm.comm.send(msg, dest=dest, tag=YT_TAG_MESSAGE)
        mylog.debug("Waiting for data.")
        MPI.Request.Wait(hook)
        return data

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        # sl = tuple(reversed(sl))
        return self._read_data_set(grid, field)[tuple(sl)]

    def terminate(self):
        msg = dict(op="end")
        if self.wg.comm.rank == 0:
            for rank in self.pool["io"].ranks:
                mylog.debug("Sending termination message to %s", rank)
                self.comm.comm.send(msg, dest=rank, tag=YT_TAG_MESSAGE)


@contextmanager
def remote_io(ds, wg, pool):
    original_io = ds.index.io
    ds.index.io = IOHandlerRemote(ds, wg, pool)
    yield
    ds.index.io.terminate()
    ds.index.io = original_io


def io_nodes(fn, n_io, n_work, func, *args, **kwargs):
    from yt.loaders import load

    pool, wg = ProcessorPool.from_sizes([(n_io, "io"), (n_work, "work")])
    rv = None
    if wg.name == "work":
        ds = load(fn)
        with remote_io(ds, wg, pool):
            rv = func(ds, *args, **kwargs)
    elif wg.name == "io":
        ds = load(fn)
        io = IOCommunicator(ds, wg, pool)
        io.wait()
    # We should broadcast the result
    rv = pool.comm.mpi_bcast(rv, root=pool["work"].ranks[0])
    pool.free_all()
    mylog.debug("Return value: %s", rv)
    return rv


# Here is an example of how to use this functionality.
if __name__ == "__main__":

    def gq(ds):
        dd = ds.all_data()
        return dd.quantities["TotalQuantity"](("gas", "cell_mass"))

    q = io_nodes("DD0087/DD0087", 8, 24, gq)
    mylog.info(q)
