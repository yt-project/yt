import itertools
import logging
import os
import sys
import traceback
from functools import wraps
from io import StringIO
from typing import List

import numpy as np
from more_itertools import always_iterable

import yt.utilities.logger
from yt.config import ytcfg
from yt.data_objects.image_array import ImageArray
from yt.funcs import is_sequence
from yt.units.unit_registry import UnitRegistry  # type: ignore
from yt.units.yt_array import YTArray
from yt.utilities.exceptions import YTNoDataInObjectError
from yt.utilities.lib.quad_tree import QuadTree, merge_quadtrees
from yt.utilities.logger import ytLogger as mylog

# We default to *no* parallelism unless it gets turned on, in which case this
# will be changed.
MPI = None
parallel_capable = False

dtype_names = dict(
    float32="MPI.FLOAT",
    float64="MPI.DOUBLE",
    int32="MPI.INT",
    int64="MPI.LONG",
    c="MPI.CHAR",
)
op_names = dict(sum="MPI.SUM", min="MPI.MIN", max="MPI.MAX")


class FilterAllMessages(logging.Filter):
    """
    This is a simple filter for logging.Logger's that won't let any
    messages pass.
    """

    def filter(self, record):
        return 0


# Set up translation table and import things


def traceback_writer_hook(file_suffix=""):
    def write_to_file(exc_type, exc, tb):
        sys.__excepthook__(exc_type, exc, tb)
        fn = f"yt_traceback{file_suffix}"
        with open(fn, "w") as fhandle:
            traceback.print_exception(exc_type, exc, tb, file=fhandle)
            print(f"Wrote traceback to {fn}")
        MPI.COMM_WORLD.Abort(1)

    return write_to_file


def default_mpi_excepthook(exception_type, exception_value, tb):
    traceback.print_tb(tb)
    mylog.error("%s: %s", exception_type.__name__, exception_value)
    comm = yt.communication_system.communicators[-1]
    if comm.size > 1:
        mylog.error("Error occurred on rank %d.", comm.rank)
    MPI.COMM_WORLD.Abort(1)


def enable_parallelism(suppress_logging: bool = False, communicator=None) -> bool:
    """
    This method is used inside a script to turn on MPI parallelism, via
    mpi4py.  More information about running yt in parallel can be found
    here: https://yt-project.org/docs/3.0/analyzing/parallel_computation.html

    Parameters
    ----------
    suppress_logging : bool
       If set to True, only rank 0 will log information after the initial
       setup of MPI.

    communicator : mpi4py.MPI.Comm
        The MPI communicator to use. This controls which processes yt can see.
        If not specified, will be set to COMM_WORLD.

    Returns
    -------
    parallel_capable: bool
        True if the call was successful. False otherwise.
    """
    global parallel_capable, MPI
    try:
        from mpi4py import MPI as _MPI
    except ImportError:
        mylog.error("Could not enable parallelism: mpi4py is not installed")
        return False
    MPI = _MPI
    exe_name = os.path.basename(sys.executable)

    # if no communicator specified, set to COMM_WORLD
    if communicator is None:
        communicator = MPI.COMM_WORLD

    parallel_capable = communicator.size > 1
    if not parallel_capable:
        mylog.error(
            "Could not enable parallelism: only one mpi process is running. "
            "To remedy this, launch the Python interpreter as\n"
            "  mpirun -n <X> python3 <yourscript>.py  # with X > 1 ",
        )
        return False

    mylog.info(
        "Global parallel computation enabled: %s / %s",
        communicator.rank,
        communicator.size,
    )
    communication_system.push(communicator)
    ytcfg["yt", "internals", "global_parallel_rank"] = communicator.rank
    ytcfg["yt", "internals", "global_parallel_size"] = communicator.size
    ytcfg["yt", "internals", "parallel"] = True
    if exe_name == "embed_enzo" or ("_parallel" in dir(sys) and sys._parallel):  # type: ignore
        ytcfg["yt", "inline"] = True
    yt.utilities.logger.uncolorize_logging()
    # Even though the uncolorize function already resets the format string,
    # we reset it again so that it includes the processor.
    f = logging.Formatter(
        "P%03i %s" % (communicator.rank, yt.utilities.logger.ufstring)
    )
    if len(yt.utilities.logger.ytLogger.handlers) > 0:
        yt.utilities.logger.ytLogger.handlers[0].setFormatter(f)

    if ytcfg.get("yt", "parallel_traceback"):
        sys.excepthook = traceback_writer_hook("_%03i" % communicator.rank)
    else:
        sys.excepthook = default_mpi_excepthook

    if ytcfg.get("yt", "log_level") < 20:
        yt.utilities.logger.ytLogger.warning(
            "Log Level is set low -- this could affect parallel performance!"
        )
    dtype_names.update(
        dict(
            float32=MPI.FLOAT,
            float64=MPI.DOUBLE,
            int32=MPI.INT,
            int64=MPI.LONG,
            c=MPI.CHAR,
        )
    )
    op_names.update(dict(sum=MPI.SUM, min=MPI.MIN, max=MPI.MAX))
    # Turn off logging on all but the root rank, if specified.
    if suppress_logging:
        if communicator.rank > 0:
            mylog.addFilter(FilterAllMessages())
    return True


# Because the dtypes will == correctly but do not hash the same, we need this
# function for dictionary access.
def get_mpi_type(dtype):
    for dt, val in dtype_names.items():
        if dt == dtype:
            return val


class ObjectIterator:
    """
    This is a generalized class that accepts a list of objects and then
    attempts to intelligently iterate over them.
    """

    def __init__(self, pobj, just_list=False, attr="_grids"):
        self.pobj = pobj
        if hasattr(pobj, attr) and getattr(pobj, attr) is not None:
            gs = getattr(pobj, attr)
        else:
            gs = getattr(pobj._data_source, attr)
        if len(gs) == 0:
            raise YTNoDataInObjectError(pobj)
        if hasattr(gs[0], "proc_num"):
            # This one sort of knows about MPI, but not quite
            self._objs = [
                g
                for g in gs
                if g.proc_num == ytcfg.get("yt", "internals", "topcomm_parallel_rank")
            ]
            self._use_all = True
        else:
            self._objs = gs
            if hasattr(self._objs[0], "filename"):
                self._objs = sorted(self._objs, key=lambda g: g.filename)
            self._use_all = False
        self.ng = len(self._objs)
        self.just_list = just_list

    def __iter__(self):
        yield from self._objs


class ParallelObjectIterator(ObjectIterator):
    """
    This takes an object, *pobj*, that implements ParallelAnalysisInterface,
    and then does its thing, calling initialize and finalize on the object.
    """

    def __init__(self, pobj, just_list=False, attr="_grids", round_robin=False):
        ObjectIterator.__init__(self, pobj, just_list, attr=attr)
        # pobj has to be a ParallelAnalysisInterface, so it must have a .comm
        # object.
        self._offset = pobj.comm.rank
        self._skip = pobj.comm.size
        # Note that we're doing this in advance, and with a simple means
        # of choosing them; more advanced methods will be explored later.
        if self._use_all:
            self.my_obj_ids = np.arange(len(self._objs))
        else:
            if not round_robin:
                self.my_obj_ids = np.array_split(
                    np.arange(len(self._objs)), self._skip
                )[self._offset]
            else:
                self.my_obj_ids = np.arange(len(self._objs))[self._offset :: self._skip]

    def __iter__(self):
        for gid in self.my_obj_ids:
            yield self._objs[gid]
        if not self.just_list:
            self.pobj._finalize_parallel()


def parallel_simple_proxy(func):
    """
    This is a decorator that broadcasts the result of computation on a single
    processor to all other processors.  To do so, it uses the _processing and
    _distributed flags in the object to check for blocks.  Meant only to be
    used on objects that subclass
    :class:`~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelAnalysisInterface`.
    """

    @wraps(func)
    def single_proc_results(self, *args, **kwargs):
        retval = None
        if hasattr(self, "dont_wrap"):
            if func.__name__ in self.dont_wrap:
                return func(self, *args, **kwargs)
        if not parallel_capable or self._processing or not self._distributed:
            return func(self, *args, **kwargs)
        comm = _get_comm((self,))
        if self._owner == comm.rank:
            self._processing = True
            retval = func(self, *args, **kwargs)
            self._processing = False
        # To be sure we utilize the root= kwarg, we manually access the .comm
        # attribute, which must be an instance of MPI.Intracomm, and call bcast
        # on that.
        retval = comm.comm.bcast(retval, root=self._owner)
        return retval

    return single_proc_results


class ParallelDummy(type):
    """
    This is a base class that, on instantiation, replaces all attributes that
    don't start with ``_`` with
    :func:`~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_simple_proxy`-wrapped
    attributes.  Used as a metaclass.
    """

    def __init__(cls, name, bases, d):
        super().__init__(name, bases, d)
        skip = d.pop("dont_wrap", [])
        extra = d.pop("extra_wrap", [])
        for attrname in d:
            if attrname.startswith("_") or attrname in skip:
                if attrname not in extra:
                    continue
            attr = getattr(cls, attrname)
            if callable(attr):
                setattr(cls, attrname, parallel_simple_proxy(attr))


def parallel_passthrough(func):
    """
    If we are not run in parallel, this function passes the input back as
    output; otherwise, the function gets called.  Used as a decorator.
    """

    @wraps(func)
    def passage(self, *args, **kwargs):
        if not self._distributed:
            return args[0]
        return func(self, *args, **kwargs)

    return passage


def _get_comm(args):
    if len(args) > 0 and hasattr(args[0], "comm"):
        comm = args[0].comm
    else:
        comm = communication_system.communicators[-1]
    return comm


def parallel_blocking_call(func):
    """
    This decorator blocks on entry and exit of a function.
    """

    @wraps(func)
    def barrierize(*args, **kwargs):
        if not parallel_capable:
            return func(*args, **kwargs)
        mylog.debug("Entering barrier before %s", func.__name__)
        comm = _get_comm(args)
        comm.barrier()
        retval = func(*args, **kwargs)
        mylog.debug("Entering barrier after %s", func.__name__)
        comm.barrier()
        return retval

    return barrierize


def parallel_root_only(func):
    """
    This decorator blocks and calls the function on the root processor,
    but does not broadcast results to the other processors.
    """

    @wraps(func)
    def root_only(*args, **kwargs):
        if not parallel_capable:
            return func(*args, **kwargs)
        comm = _get_comm(args)
        rv = None
        if comm.rank == 0:
            try:
                rv = func(*args, **kwargs)
                all_clear = 1
            except Exception:
                traceback.print_last()
                all_clear = 0
        else:
            all_clear = None
        all_clear = comm.mpi_bcast(all_clear)
        if not all_clear:
            raise RuntimeError
        return rv

    return root_only


class Workgroup:
    def __init__(self, size, ranks, comm, name):
        self.size = size
        self.ranks = ranks
        self.comm = comm
        self.name = name


class ProcessorPool:
    comm = None
    size = None
    ranks = None
    available_ranks = None
    tasks = None

    def __init__(self):
        self.comm = communication_system.communicators[-1]
        self.size = self.comm.size
        self.ranks = list(range(self.size))
        self.available_ranks = list(range(self.size))
        self.workgroups = []

    def add_workgroup(self, size=None, ranks=None, name=None):
        if size is None:
            size = len(self.available_ranks)
        if len(self.available_ranks) < size:
            mylog.error(
                "Not enough resources available, asked for %d have %d",
                size,
                self.available_ranks,
            )
            raise RuntimeError
        if ranks is None:
            ranks = [self.available_ranks.pop(0) for i in range(size)]
        # Default name to the workgroup number.
        if name is None:
            name = str(len(self.workgroups))
        group = self.comm.comm.Get_group().Incl(ranks)
        new_comm = self.comm.comm.Create(group)
        if self.comm.rank in ranks:
            communication_system.communicators.append(Communicator(new_comm))
        self.workgroups.append(Workgroup(len(ranks), ranks, new_comm, name))

    def free_workgroup(self, workgroup):
        # If you want to actually delete the workgroup you will need to
        # pop it out of the self.workgroups list so you don't have references
        # that are left dangling, e.g. see free_all() below.
        for i in workgroup.ranks:
            if self.comm.rank == i:
                communication_system.communicators.pop()
            self.available_ranks.append(i)
        self.available_ranks.sort()

    def free_all(self):
        for wg in self.workgroups:
            self.free_workgroup(wg)
        while self.workgroups:
            self.workgroups.pop(0)

    @classmethod
    def from_sizes(cls, sizes):
        pool = cls()
        rank = pool.comm.rank
        for i, size in enumerate(always_iterable(sizes)):
            if is_sequence(size):
                size, name = size
            else:
                name = "workgroup_%02i" % i
            pool.add_workgroup(size, name=name)
        for wg in pool.workgroups:
            if rank in wg.ranks:
                workgroup = wg
        return pool, workgroup

    def __getitem__(self, key):
        for wg in self.workgroups:
            if wg.name == key:
                return wg
        raise KeyError(key)


class ResultsStorage:
    slots = ["result", "result_id"]
    result = None
    result_id = None


def parallel_objects(objects, njobs=0, storage=None, barrier=True, dynamic=False):
    r"""This function dispatches components of an iterable to different
    processors.

    The parallel_objects function accepts an iterable, *objects*, and based on
    the number of jobs requested and number of available processors, decides
    how to dispatch individual objects to processors or sets of processors.
    This can implicitly include multi-level parallelism, such that the
    processor groups assigned each object can be composed of several or even
    hundreds of processors.  *storage* is also available, for collation of
    results at the end of the iteration loop.

    Calls to this function can be nested.

    This should not be used to iterate over datasets --
    :class:`~yt.data_objects.time_series.DatasetSeries` provides a much nicer
    interface for that.

    Parameters
    ----------
    objects : Iterable
        The list of objects to dispatch to different processors.
    njobs : int
        How many jobs to spawn.  By default, one job will be dispatched for
        each available processor.
    storage : dict
        This is a dictionary, which will be filled with results during the
        course of the iteration.  The keys will be the dataset
        indices and the values will be whatever is assigned to the *result*
        attribute on the storage during iteration.
    barrier : bool
        Should a barier be placed at the end of iteration?
    dynamic : bool
        This governs whether or not dynamic load balancing will be enabled.
        This requires one dedicated processor; if this is enabled with a set of
        128 processors available, only 127 will be available to iterate over
        objects as one will be load balancing the rest.


    Examples
    --------
    Here is a simple example of iterating over a set of centers and making
    slice plots centered at each.

    >>> for c in parallel_objects(centers):
    ...     SlicePlot(ds, "x", "Density", center=c).save()
    ...

    Here's an example of calculating the angular momentum vector of a set of
    spheres, but with a set of four jobs of multiple processors each.  Note
    that we also store the results.

    >>> storage = {}
    >>> for sto, c in parallel_objects(centers, njobs=4, storage=storage):
    ...     sp = ds.sphere(c, (100, "kpc"))
    ...     sto.result = sp.quantities["AngularMomentumVector"]()
    ...
    >>> for sphere_id, L in sorted(storage.items()):
    ...     print(centers[sphere_id], L)
    ...

    """
    if dynamic:
        from .task_queue import dynamic_parallel_objects

        yield from dynamic_parallel_objects(objects, njobs=njobs, storage=storage)
        return

    if not parallel_capable:
        njobs = 1
    my_communicator = communication_system.communicators[-1]
    my_size = my_communicator.size
    if njobs <= 0:
        njobs = my_size
    if njobs > my_size:
        mylog.error(
            "You have asked for %s jobs, but you only have %s processors.",
            njobs,
            my_size,
        )
        raise RuntimeError
    my_rank = my_communicator.rank
    all_new_comms = np.array_split(np.arange(my_size), njobs)
    for i, comm_set in enumerate(all_new_comms):
        if my_rank in comm_set:
            my_new_id = i
            break
    if parallel_capable:
        communication_system.push_with_ids(all_new_comms[my_new_id].tolist())
    to_share = {}
    # If our objects object is slice-aware, like time series data objects are,
    # this will prevent intermediate objects from being created.
    oiter = itertools.islice(enumerate(objects), my_new_id, None, njobs)
    for result_id, obj in oiter:
        if storage is not None:
            rstore = ResultsStorage()
            rstore.result_id = result_id
            yield rstore, obj
            to_share[rstore.result_id] = rstore.result
        else:
            yield obj
    if parallel_capable:
        communication_system.pop()
    if storage is not None:
        # Now we have to broadcast it
        new_storage = my_communicator.par_combine_object(
            to_share, datatype="dict", op="join"
        )
        storage.update(new_storage)
    if barrier:
        my_communicator.barrier()


def parallel_ring(objects, generator_func, mutable=False):
    r"""This function loops in a ring around a set of objects, yielding the
    results of generator_func and passing from one processor to another to
    avoid IO or expensive computation.

    This function is designed to operate in sequence on a set of objects, where
    the creation of those objects might be expensive.  For instance, this could
    be a set of particles that are costly to read from disk.  Processor N will
    run generator_func on an object, and the results of that will both be
    yielded and passed to processor N-1.  If the length of the objects is not
    equal to the number of processors, then the final processor in the top
    communicator will re-generate the data as needed.

    In all likelihood, this function will only be useful internally to yt.

    Parameters
    ----------
    objects : Iterable
        The list of objects to operate on.
    generator_func : callable
        This function will be called on each object, and the results yielded.
        It must return a single NumPy array; for multiple values, it needs to
        have a custom dtype.
    mutable : bool
        Should the arrays be considered mutable?  Currently, this will only
        work if the number of processors equals the number of objects.

    Examples
    --------
    Here is a simple example of a ring loop around a set of integers, with a
    custom dtype.

    >>> dt = np.dtype([("x", "float64"), ("y", "float64"), ("z", "float64")])
    >>> def gfunc(o):
    ...     np.random.seed(o)
    ...     rv = np.empty(1000, dtype=dt)
    ...     rv["x"] = np.random.random(1000)
    ...     rv["y"] = np.random.random(1000)
    ...     rv["z"] = np.random.random(1000)
    ...     return rv
    ...
    >>> obj = range(8)
    >>> for obj, arr in parallel_ring(obj, gfunc):
    ...     print(arr["x"].sum(), arr["y"].sum(), arr["z"].sum())
    ...

    """
    if mutable:
        raise NotImplementedError
    my_comm = communication_system.communicators[-1]
    my_size = my_comm.size
    my_rank = my_comm.rank  # This will also be the first object we access
    if not parallel_capable and not mutable:
        for obj in objects:
            yield obj, generator_func(obj)
        return
    generate_endpoints = len(objects) != my_size
    # gback False: send the object backwards
    # gforw False: receive an object from forwards
    if len(objects) == my_size:
        generate_endpoints = False
        gback = False
        gforw = False
    else:
        # In this case, the first processor (my_rank == 0) will generate.
        generate_endpoints = True
        gback = my_rank == 0
        gforw = my_rank == my_size - 1
    if generate_endpoints and mutable:
        raise NotImplementedError
    # Now we need to do pairwise sends
    source = (my_rank + 1) % my_size
    dest = (my_rank - 1) % my_size
    oiter = itertools.islice(itertools.cycle(objects), my_rank, my_rank + len(objects))
    idata = None
    isize = np.zeros((1,), dtype="int64")
    osize = np.zeros((1,), dtype="int64")
    for obj in oiter:
        if idata is None or gforw:
            idata = generator_func(obj)
            idtype = odtype = idata.dtype
            if get_mpi_type(idtype) is None:
                idtype = "c"
        yield obj, idata
        # We first send to the previous processor
        tags = []
        if not gforw:
            tags.append(my_comm.mpi_nonblocking_recv(isize, source))
        if not gback:
            osize[0] = idata.size
            tags.append(my_comm.mpi_nonblocking_send(osize, dest))
        my_comm.mpi_Request_Waitall(tags)
        odata = idata
        tags = []
        if not gforw:
            idata = np.empty(isize[0], dtype=odtype)
            tags.append(
                my_comm.mpi_nonblocking_recv(idata.view(idtype), source, dtype=idtype)
            )
        if not gback:
            tags.append(
                my_comm.mpi_nonblocking_send(odata.view(idtype), dest, dtype=idtype)
            )
        my_comm.mpi_Request_Waitall(tags)
        del odata


class CommunicationSystem:
    communicators: List["Communicator"] = []

    def __init__(self):
        self.communicators.append(Communicator(None))

    def push(self, new_comm):
        if not isinstance(new_comm, Communicator):
            new_comm = Communicator(new_comm)
        self.communicators.append(new_comm)
        self._update_parallel_state(new_comm)

    def push_with_ids(self, ids):
        group = self.communicators[-1].comm.Get_group().Incl(ids)
        new_comm = self.communicators[-1].comm.Create(group)
        self.push(new_comm)
        return new_comm

    def _update_parallel_state(self, new_comm):
        ytcfg["yt", "internals", "topcomm_parallel_size"] = new_comm.size
        ytcfg["yt", "internals", "topcomm_parallel_rank"] = new_comm.rank
        if new_comm.rank > 0 and ytcfg.get("yt", "serialize"):
            ytcfg["yt", "only_deserialize"] = True

    def pop(self):
        self.communicators.pop()
        self._update_parallel_state(self.communicators[-1])


def _reconstruct_communicator():
    return communication_system.communicators[-1]


class Communicator:
    comm = None
    _grids = None
    _distributed = None
    __tocast = "c"

    def __init__(self, comm=None):
        self.comm = comm
        self._distributed = comm is not None and self.comm.size > 1

    def __del__(self):
        if self.comm is not None:
            self.comm.Free()

    """
    This is an interface specification providing several useful utility
    functions for analyzing something in parallel.
    """

    def __reduce__(self):
        # We don't try to reconstruct any of the properties of the communicator
        # or the processors.  In general, we don't want to.
        return (_reconstruct_communicator, ())

    def barrier(self):
        if not self._distributed:
            return
        mylog.debug("Opening MPI Barrier on %s", self.comm.rank)
        self.comm.Barrier()

    def mpi_exit_test(self, data=False):
        # data==True -> exit. data==False -> no exit
        mine, statuses = self.mpi_info_dict(data)
        if True in statuses.values():
            raise RuntimeError("Fatal error. Exiting.")
        return None

    @parallel_passthrough
    def par_combine_object(self, data, op, datatype=None):
        # op can be chosen from:
        #   cat
        #   join
        # data is selected to be of types:
        #   np.ndarray
        #   dict
        #   data field dict
        if datatype is not None:
            pass
        elif isinstance(data, dict):
            datatype = "dict"
        elif isinstance(data, np.ndarray):
            datatype = "array"
        elif isinstance(data, list):
            datatype = "list"
        # Now we have our datatype, and we conduct our operation
        if datatype == "dict" and op == "join":
            if self.comm.rank == 0:
                for i in range(1, self.comm.size):
                    data.update(self.comm.recv(source=i, tag=0))
            else:
                self.comm.send(data, dest=0, tag=0)

            # Send the keys first, then each item one by one
            # This is to prevent MPI from crashing when sending more
            # than 2GiB of data over the network.
            keys = self.comm.bcast(list(data.keys()), root=0)
            for key in keys:
                tmp = data.get(key, None)
                data[key] = self.comm.bcast(tmp, root=0)
            return data
        elif datatype == "dict" and op == "cat":
            field_keys = sorted(data.keys())
            size = data[field_keys[0]].shape[-1]
            sizes = np.zeros(self.comm.size, dtype="int64")
            outsize = np.array(size, dtype="int64")
            self.comm.Allgather([outsize, 1, MPI.LONG], [sizes, 1, MPI.LONG])
            # This nested concatenate is to get the shapes to work out correctly;
            # if we just add [0] to sizes, it will broadcast a summation, not a
            # concatenation.
            offsets = np.add.accumulate(np.concatenate([[0], sizes]))[:-1]
            arr_size = self.comm.allreduce(size, op=MPI.SUM)
            for key in field_keys:
                dd = data[key]
                rv = self.alltoallv_array(dd, arr_size, offsets, sizes)
                data[key] = rv
            return data
        elif datatype == "array" and op == "cat":
            if data is None:
                ncols = -1
                size = 0
                dtype = "float64"
                mylog.warning(
                    "Array passed to par_combine_object was None. "
                    "Setting dtype to float64. This may break things!"
                )
            else:
                dtype = data.dtype
                if len(data) == 0:
                    ncols = -1
                    size = 0
                elif len(data.shape) == 1:
                    ncols = 1
                    size = data.shape[0]
                else:
                    ncols, size = data.shape
            ncols = self.comm.allreduce(ncols, op=MPI.MAX)
            if ncols == 0:
                data = np.zeros(0, dtype=dtype)  # This only works for
            elif data is None:
                data = np.zeros((ncols, 0), dtype=dtype)
            size = data.shape[-1]
            sizes = np.zeros(self.comm.size, dtype="int64")
            outsize = np.array(size, dtype="int64")
            self.comm.Allgather([outsize, 1, MPI.LONG], [sizes, 1, MPI.LONG])
            # This nested concatenate is to get the shapes to work out correctly;
            # if we just add [0] to sizes, it will broadcast a summation, not a
            # concatenation.
            offsets = np.add.accumulate(np.concatenate([[0], sizes]))[:-1]
            arr_size = self.comm.allreduce(size, op=MPI.SUM)
            data = self.alltoallv_array(data, arr_size, offsets, sizes)
            return data
        elif datatype == "list" and op == "cat":
            recv_data = self.comm.allgather(data)
            # Now flatten into a single list, since this
            # returns us a list of lists.
            data = []
            while recv_data:
                data.extend(recv_data.pop(0))
            return data
        raise NotImplementedError

    @parallel_passthrough
    def mpi_bcast(self, data, root=0):
        # The second check below makes sure that we know how to communicate
        # this type of array. Otherwise, we'll pickle it.
        if isinstance(data, np.ndarray) and get_mpi_type(data.dtype) is not None:
            if self.comm.rank == root:
                if isinstance(data, YTArray):
                    info = (
                        data.shape,
                        data.dtype,
                        str(data.units),
                        data.units.registry.lut,
                    )
                    if isinstance(data, ImageArray):
                        info += ("ImageArray",)
                    else:
                        info += ("YTArray",)
                else:
                    info = (data.shape, data.dtype)
            else:
                info = ()
            info = self.comm.bcast(info, root=root)
            if self.comm.rank != root:
                if len(info) == 5:
                    registry = UnitRegistry(lut=info[3], add_default_symbols=False)
                    if info[-1] == "ImageArray":
                        data = ImageArray(
                            np.empty(info[0], dtype=info[1]),
                            units=info[2],
                            registry=registry,
                        )
                    else:
                        data = YTArray(
                            np.empty(info[0], dtype=info[1]), info[2], registry=registry
                        )
                else:
                    data = np.empty(info[0], dtype=info[1])
            mpi_type = get_mpi_type(info[1])
            self.comm.Bcast([data, mpi_type], root=root)
            return data
        else:
            # Use pickled methods.
            data = self.comm.bcast(data, root=root)
            return data

    def preload(self, grids, fields, io_handler):
        # This is non-functional.
        return

    @parallel_passthrough
    def mpi_allreduce(self, data, dtype=None, op="sum"):
        op = op_names[op]
        if isinstance(data, np.ndarray) and data.dtype != bool:
            if dtype is None:
                dtype = data.dtype
            if dtype != data.dtype:
                data = data.astype(dtype)
            temp = data.copy()
            self.comm.Allreduce(
                [temp, get_mpi_type(dtype)], [data, get_mpi_type(dtype)], op
            )
            return data
        else:
            # We use old-school pickling here on the assumption the arrays are
            # relatively small ( < 1e7 elements )
            return self.comm.allreduce(data, op)

    ###
    # Non-blocking stuff.
    ###

    def mpi_nonblocking_recv(self, data, source, tag=0, dtype=None):
        if not self._distributed:
            return -1
        if dtype is None:
            dtype = data.dtype
        mpi_type = get_mpi_type(dtype)
        return self.comm.Irecv([data, mpi_type], source, tag)

    def mpi_nonblocking_send(self, data, dest, tag=0, dtype=None):
        if not self._distributed:
            return -1
        if dtype is None:
            dtype = data.dtype
        mpi_type = get_mpi_type(dtype)
        return self.comm.Isend([data, mpi_type], dest, tag)

    def mpi_Request_Waitall(self, hooks):
        if not self._distributed:
            return
        MPI.Request.Waitall(hooks)

    def mpi_Request_Waititer(self, hooks):
        for _hook in hooks:
            req = MPI.Request.Waitany(hooks)
            yield req

    def mpi_Request_Testall(self, hooks):
        """
        This returns False if any of the request hooks are un-finished,
        and True if they are all finished.
        """
        if not self._distributed:
            return True
        return MPI.Request.Testall(hooks)

    ###
    # End non-blocking stuff.
    ###

    ###
    # Parallel rank and size properties.
    ###

    @property
    def size(self):
        if not self._distributed:
            return 1
        return self.comm.size

    @property
    def rank(self):
        if not self._distributed:
            return 0
        return self.comm.rank

    def mpi_info_dict(self, info):
        if not self._distributed:
            return 0, {0: info}
        data = None
        if self.comm.rank == 0:
            data = {0: info}
            for i in range(1, self.comm.size):
                data[i] = self.comm.recv(source=i, tag=0)
        else:
            self.comm.send(info, dest=0, tag=0)
        mylog.debug("Opening MPI Broadcast on %s", self.comm.rank)
        data = self.comm.bcast(data, root=0)
        return self.comm.rank, data

    def claim_object(self, obj):
        if not self._distributed:
            return
        obj._owner = self.comm.rank
        obj._distributed = True

    def do_not_claim_object(self, obj):
        if not self._distributed:
            return
        obj._owner = -1
        obj._distributed = True

    def write_on_root(self, fn):
        if not self._distributed:
            return open(fn, "w")
        if self.comm.rank == 0:
            return open(fn, "w")
        else:
            return StringIO()

    def get_filename(self, prefix, rank=None):
        if not self._distributed:
            return prefix
        if rank is None:
            return "%s_%04i" % (prefix, self.comm.rank)
        else:
            return "%s_%04i" % (prefix, rank)

    def is_mine(self, obj):
        if not obj._distributed:
            return True
        return obj._owner == self.comm.rank

    def send_quadtree(self, target, buf, tgd, args):
        sizebuf = np.zeros(1, "int64")
        sizebuf[0] = buf[0].size
        self.comm.Send([sizebuf, MPI.LONG], dest=target)
        self.comm.Send([buf[0], MPI.INT], dest=target)
        self.comm.Send([buf[1], MPI.DOUBLE], dest=target)
        self.comm.Send([buf[2], MPI.DOUBLE], dest=target)

    def recv_quadtree(self, target, tgd, args):
        sizebuf = np.zeros(1, "int64")
        self.comm.Recv(sizebuf, source=target)
        buf = [
            np.empty((sizebuf[0],), "int32"),
            np.empty((sizebuf[0], args[2]), "float64"),
            np.empty((sizebuf[0],), "float64"),
        ]
        self.comm.Recv([buf[0], MPI.INT], source=target)
        self.comm.Recv([buf[1], MPI.DOUBLE], source=target)
        self.comm.Recv([buf[2], MPI.DOUBLE], source=target)
        return buf

    @parallel_passthrough
    def merge_quadtree_buffers(self, qt, merge_style):
        # This is a modified version of pairwise reduction from Lisandro Dalcin,
        # in the reductions demo of mpi4py
        size = self.comm.size
        rank = self.comm.rank

        mask = 1

        buf = qt.tobuffer()
        print("PROC", rank, buf[0].shape, buf[1].shape, buf[2].shape)
        sys.exit()

        args = qt.get_args()  # Will always be the same
        tgd = np.array([args[0], args[1]], dtype="int64")
        sizebuf = np.zeros(1, "int64")

        while mask < size:
            if (mask & rank) != 0:
                target = (rank & ~mask) % size
                # print("SENDING FROM %02i to %02i" % (rank, target))
                buf = qt.tobuffer()
                self.send_quadtree(target, buf, tgd, args)
                # qt = self.recv_quadtree(target, tgd, args)
            else:
                target = rank | mask
                if target < size:
                    # print("RECEIVING FROM %02i on %02i" % (target, rank))
                    buf = self.recv_quadtree(target, tgd, args)
                    qto = QuadTree(tgd, args[2], qt.bounds)
                    qto.frombuffer(buf[0], buf[1], buf[2], merge_style)
                    merge_quadtrees(qt, qto, style=merge_style)
                    del qto
                    # self.send_quadtree(target, qt, tgd, args)
            mask <<= 1

        if rank == 0:
            buf = qt.tobuffer()
            sizebuf[0] = buf[0].size
        self.comm.Bcast([sizebuf, MPI.LONG], root=0)
        if rank != 0:
            buf = [
                np.empty((sizebuf[0],), "int32"),
                np.empty((sizebuf[0], args[2]), "float64"),
                np.empty((sizebuf[0],), "float64"),
            ]
        self.comm.Bcast([buf[0], MPI.INT], root=0)
        self.comm.Bcast([buf[1], MPI.DOUBLE], root=0)
        self.comm.Bcast([buf[2], MPI.DOUBLE], root=0)
        self.refined = buf[0]
        if rank != 0:
            qt = QuadTree(tgd, args[2], qt.bounds)
            qt.frombuffer(buf[0], buf[1], buf[2], merge_style)
        return qt

    def send_array(self, arr, dest, tag=0):
        if not isinstance(arr, np.ndarray):
            self.comm.send((None, None), dest=dest, tag=tag)
            self.comm.send(arr, dest=dest, tag=tag)
            return
        tmp = arr.view(self.__tocast)  # Cast to CHAR
        # communicate type and shape and optionally units
        if isinstance(arr, YTArray):
            unit_metadata = (str(arr.units), arr.units.registry.lut)
            if isinstance(arr, ImageArray):
                unit_metadata += ("ImageArray",)
            else:
                unit_metadata += ("YTArray",)
        else:
            unit_metadata = ()
        self.comm.send((arr.dtype.str, arr.shape) + unit_metadata, dest=dest, tag=tag)
        self.comm.Send([arr, MPI.CHAR], dest=dest, tag=tag)
        del tmp

    def recv_array(self, source, tag=0):
        metadata = self.comm.recv(source=source, tag=tag)
        dt, ne = metadata[:2]
        if ne is None and dt is None:
            return self.comm.recv(source=source, tag=tag)
        arr = np.empty(ne, dtype=dt)
        if len(metadata) == 5:
            registry = UnitRegistry(lut=metadata[3], add_default_symbols=False)
            if metadata[-1] == "ImageArray":
                arr = ImageArray(arr, units=metadata[2], registry=registry)
            else:
                arr = YTArray(arr, metadata[2], registry=registry)
        tmp = arr.view(self.__tocast)
        self.comm.Recv([tmp, MPI.CHAR], source=source, tag=tag)
        return arr

    def alltoallv_array(self, send, total_size, offsets, sizes):
        if len(send.shape) > 1:
            recv = []
            for i in range(send.shape[0]):
                recv.append(
                    self.alltoallv_array(send[i, :].copy(), total_size, offsets, sizes)
                )
            recv = np.array(recv)
            return recv
        offset = offsets[self.comm.rank]
        tmp_send = send.view(self.__tocast)
        recv = np.empty(total_size, dtype=send.dtype)
        if isinstance(send, YTArray):
            # We assume send.units is consistent with the units
            # on the receiving end.
            if isinstance(send, ImageArray):
                recv = ImageArray(recv, units=send.units)
            else:
                recv = YTArray(recv, send.units)
        recv[offset : offset + send.size] = send[:]
        dtr = send.dtype.itemsize / tmp_send.dtype.itemsize  # > 1
        roff = [off * dtr for off in offsets]
        rsize = [siz * dtr for siz in sizes]
        tmp_recv = recv.view(self.__tocast)
        self.comm.Allgatherv(
            (tmp_send, tmp_send.size, MPI.CHAR), (tmp_recv, (rsize, roff), MPI.CHAR)
        )
        return recv

    def probe_loop(self, tag, callback):
        while True:
            st = MPI.Status()
            self.comm.Probe(MPI.ANY_SOURCE, tag=tag, status=st)
            try:
                callback(st)
            except StopIteration:
                mylog.debug("Probe loop ending.")
                break


communication_system = CommunicationSystem()


class ParallelAnalysisInterface:
    comm = None
    _grids = None
    _distributed = None

    def __init__(self, comm=None):
        if comm is None:
            self.comm = communication_system.communicators[-1]
        else:
            self.comm = comm
        self._grids = self.comm._grids
        self._distributed = self.comm._distributed

    def _get_objs(self, attr, *args, **kwargs):
        if self._distributed:
            rr = kwargs.pop("round_robin", False)
            self._initialize_parallel(*args, **kwargs)
            return ParallelObjectIterator(self, attr=attr, round_robin=rr)
        return ObjectIterator(self, attr=attr)

    def _get_grids(self, *args, **kwargs):
        if self._distributed:
            self._initialize_parallel(*args, **kwargs)
            return ParallelObjectIterator(self, attr="_grids")
        return ObjectIterator(self, attr="_grids")

    def _get_grid_objs(self):
        if self._distributed:
            return ParallelObjectIterator(self, True, attr="_grids")
        return ObjectIterator(self, True, attr="_grids")

    def get_dependencies(self, fields):
        deps = []
        fi = self.ds.field_info
        for field in fields:
            if any(getattr(v, "ghost_zones", 0) > 0 for v in fi[field].validators):
                continue
            deps += list(
                always_iterable(fi[field].get_dependencies(ds=self.ds).requested)
            )
        return list(set(deps))

    def _initialize_parallel(self):
        pass

    def _finalize_parallel(self):
        pass

    def partition_index_2d(self, axis):
        if not self._distributed:
            return False, self.index.grid_collection(self.center, self.index.grids)

        xax = self.ds.coordinates.x_axis[axis]
        yax = self.ds.coordinates.y_axis[axis]
        cc = MPI.Compute_dims(self.comm.size, 2)
        mi = self.comm.rank
        cx, cy = np.unravel_index(mi, cc)
        x = np.mgrid[0 : 1 : (cc[0] + 1) * 1j][cx : cx + 2]
        y = np.mgrid[0 : 1 : (cc[1] + 1) * 1j][cy : cy + 2]

        DLE, DRE = self.ds.domain_left_edge.copy(), self.ds.domain_right_edge.copy()
        LE = np.ones(3, dtype="float64") * DLE
        RE = np.ones(3, dtype="float64") * DRE
        LE[xax] = x[0] * (DRE[xax] - DLE[xax]) + DLE[xax]
        RE[xax] = x[1] * (DRE[xax] - DLE[xax]) + DLE[xax]
        LE[yax] = y[0] * (DRE[yax] - DLE[yax]) + DLE[yax]
        RE[yax] = y[1] * (DRE[yax] - DLE[yax]) + DLE[yax]
        mylog.debug("Dimensions: %s %s", LE, RE)

        reg = self.ds.region(self.center, LE, RE)
        return True, reg

    def partition_index_3d(self, ds, padding=0.0, rank_ratio=1):
        LE, RE = np.array(ds.left_edge), np.array(ds.right_edge)
        # We need to establish if we're looking at a subvolume, in which case
        # we *do* want to pad things.
        if (LE == self.ds.domain_left_edge).all() and (
            RE == self.ds.domain_right_edge
        ).all():
            subvol = False
        else:
            subvol = True
        if not self._distributed and not subvol:
            return False, LE, RE, ds
        if not self._distributed and subvol:
            return True, LE, RE, self.ds.region(self.center, LE - padding, RE + padding)
        elif ytcfg.get("yt", "inline"):
            # At this point, we want to identify the root grid tile to which
            # this processor is assigned.
            # The only way I really know how to do this is to get the level-0
            # grid that belongs to this processor.
            grids = self.ds.index.select_grids(0)
            root_grids = [g for g in grids if g.proc_num == self.comm.rank]
            if len(root_grids) != 1:
                raise RuntimeError
            # raise KeyError
            LE = root_grids[0].LeftEdge
            RE = root_grids[0].RightEdge
            return True, LE, RE, self.ds.region(self.center, LE, RE)

        cc = MPI.Compute_dims(self.comm.size / rank_ratio, 3)
        mi = self.comm.rank % (self.comm.size // rank_ratio)
        cx, cy, cz = np.unravel_index(mi, cc)
        x = np.mgrid[LE[0] : RE[0] : (cc[0] + 1) * 1j][cx : cx + 2]
        y = np.mgrid[LE[1] : RE[1] : (cc[1] + 1) * 1j][cy : cy + 2]
        z = np.mgrid[LE[2] : RE[2] : (cc[2] + 1) * 1j][cz : cz + 2]

        LE = np.array([x[0], y[0], z[0]], dtype="float64")
        RE = np.array([x[1], y[1], z[1]], dtype="float64")

        if padding > 0:
            return True, LE, RE, self.ds.region(self.center, LE - padding, RE + padding)

        return False, LE, RE, self.ds.region(self.center, LE, RE)

    def partition_region_3d(self, left_edge, right_edge, padding=0.0, rank_ratio=1):
        """
        Given a region, it subdivides it into smaller regions for parallel
        analysis.
        """
        LE, RE = left_edge[:], right_edge[:]
        if not self._distributed:
            raise NotImplementedError

        cc = MPI.Compute_dims(self.comm.size / rank_ratio, 3)
        mi = self.comm.rank % (self.comm.size // rank_ratio)
        cx, cy, cz = np.unravel_index(mi, cc)
        x = np.mgrid[LE[0] : RE[0] : (cc[0] + 1) * 1j][cx : cx + 2]
        y = np.mgrid[LE[1] : RE[1] : (cc[1] + 1) * 1j][cy : cy + 2]
        z = np.mgrid[LE[2] : RE[2] : (cc[2] + 1) * 1j][cz : cz + 2]

        LE = np.array([x[0], y[0], z[0]], dtype="float64")
        RE = np.array([x[1], y[1], z[1]], dtype="float64")

        if padding > 0:
            return True, LE, RE, self.ds.region(self.center, LE - padding, RE + padding)

        return False, LE, RE, self.ds.region(self.center, LE, RE)

    def partition_index_3d_bisection_list(self):
        """
        Returns an array that is used to drive _partition_index_3d_bisection,
        below.
        """

        def factor(n):
            if n == 1:
                return [1]
            i = 2
            limit = n**0.5
            while i <= limit:
                if n % i == 0:
                    ret = factor(n / i)
                    ret.append(i)
                    return ret
                i += 1
            return [n]

        cc = MPI.Compute_dims(self.comm.size, 3)
        si = self.comm.size

        factors = factor(si)
        xyzfactors = [factor(cc[0]), factor(cc[1]), factor(cc[2])]

        # Each entry of cuts is a two element list, that is:
        # [cut dim, number of cuts]
        cuts = []
        # The higher cuts are in the beginning.
        # We're going to do our best to make the cuts cyclic, i.e. x, then y,
        # then z, etc...
        lastdim = 0
        for f in factors:
            nextdim = (lastdim + 1) % 3
            while True:
                if f in xyzfactors[nextdim]:
                    cuts.append([nextdim, f])
                    topop = xyzfactors[nextdim].index(f)
                    xyzfactors[nextdim].pop(topop)
                    lastdim = nextdim
                    break
                nextdim = (nextdim + 1) % 3
        return cuts


class GroupOwnership(ParallelAnalysisInterface):
    def __init__(self, items):
        ParallelAnalysisInterface.__init__(self)
        self.num_items = len(items)
        self.items = items
        assert self.num_items >= self.comm.size
        self.owned = list(range(self.comm.size))
        self.pointer = 0
        if parallel_capable:
            communication_system.push_with_ids([self.comm.rank])

    def __del__(self):
        if parallel_capable:
            communication_system.pop()

    def inc(self, n=-1):
        old_item = self.item
        if n == -1:
            n = self.comm.size
        for _ in range(n):
            if self.pointer >= self.num_items - self.comm.size:
                break
            self.owned[self.pointer % self.comm.size] += self.comm.size
            self.pointer += 1
        if self.item is not old_item:
            self.switch()

    def dec(self, n=-1):
        old_item = self.item
        if n == -1:
            n = self.comm.size
        for _ in range(n):
            if self.pointer == 0:
                break
            self.owned[(self.pointer - 1) % self.comm.size] -= self.comm.size
            self.pointer -= 1
        if self.item is not old_item:
            self.switch()

    _last = None

    @property
    def item(self):
        own = self.owned[self.comm.rank]
        if self._last != own:
            self._item = self.items[own]
            self._last = own
        return self._item

    def switch(self):
        pass
