import abc
import os
import weakref
from typing import Optional, Tuple

import numpy as np

from yt.config import ytcfg
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.units.yt_array import YTArray
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.io_handler import io_registry
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    ParallelAnalysisInterface,
    parallel_root_only,
)


class Index(ParallelAnalysisInterface, abc.ABC):
    """The base index class"""

    _unsupported_objects: Tuple[str, ...] = ()
    _index_properties: Tuple[str, ...] = ()

    def __init__(self, ds, dataset_type):
        ParallelAnalysisInterface.__init__(self)
        self.dataset = weakref.proxy(ds)
        self.ds = self.dataset

        self._initialize_state_variables()

        mylog.debug("Initializing data storage.")
        self._initialize_data_storage()

        mylog.debug("Setting up domain geometry.")
        self._setup_geometry()

        mylog.debug("Initializing data grid data IO")
        self._setup_data_io()

        # Note that this falls under the "geometry" object since it's
        # potentially quite expensive, and should be done with the indexing.
        mylog.debug("Detecting fields.")
        self._detect_output_fields()

    @abc.abstractmethod
    def _detect_output_fields(self):
        pass

    def _icoords_to_fcoords(
        self,
        icoords: np.ndarray,
        ires: np.ndarray,
        axes: Optional[Tuple[int, ...]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        # What's the use of raising NotImplementedError for this, when it's an
        # abstract base class?  Well, only *some* of the subclasses have it --
        # and for those that *don't*, we should not be calling it -- and since
        # it's a semi-private method, it shouldn't be called outside of yt
        # machinery.  So we shouldn't ever get here!
        raise NotImplementedError

    def _initialize_state_variables(self):
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self.num_grids = None

    def _initialize_data_storage(self):
        if not ytcfg.get("yt", "serialize"):
            return
        fn = self.ds.storage_filename
        if fn is None:
            if os.path.isfile(
                os.path.join(self.directory, f"{self.ds.unique_identifier}.yt")
            ):
                fn = os.path.join(self.directory, f"{self.ds.unique_identifier}.yt")
            else:
                fn = os.path.join(self.directory, f"{self.dataset.basename}.yt")
        dir_to_check = os.path.dirname(fn)
        if dir_to_check == "":
            dir_to_check = "."
        # We have four options:
        #    Writeable, does not exist      : create, open as append
        #    Writeable, does exist          : open as append
        #    Not writeable, does not exist  : do not attempt to open
        #    Not writeable, does exist      : open as read-only
        exists = os.path.isfile(fn)
        if not exists:
            writeable = os.access(dir_to_check, os.W_OK)
        else:
            writeable = os.access(fn, os.W_OK)
        writeable = writeable and not ytcfg.get("yt", "only_deserialize")
        # We now have our conditional stuff
        self.comm.barrier()
        if not writeable and not exists:
            return
        if writeable:
            try:
                if not exists:
                    self.__create_data_file(fn)
                self._data_mode = "a"
            except OSError:
                self._data_mode = None
                return
        else:
            self._data_mode = "r"

        self.__data_filename = fn
        self._data_file = h5py.File(fn, mode=self._data_mode)

    def __create_data_file(self, fn):
        # Note that this used to be parallel_root_only; it no longer is,
        # because we have better logic to decide who owns the file.
        f = h5py.File(fn, mode="a")
        f.close()

    def _setup_data_io(self):
        if getattr(self, "io", None) is not None:
            return
        self.io = io_registry[self.dataset_type](self.dataset)

    @parallel_root_only
    def save_data(
        self, array, node, name, set_attr=None, force=False, passthrough=False
    ):
        """
        Arbitrary numpy data will be saved to the region in the datafile
        described by *node* and *name*.  If data file does not exist, it throws
        no error and simply does not save.
        """

        if self._data_mode != "a":
            return
        try:
            node_loc = self._data_file[node]
            if name in node_loc and force:
                mylog.info("Overwriting node %s/%s", node, name)
                del self._data_file[node][name]
            elif name in node_loc and passthrough:
                return
        except Exception:
            pass
        myGroup = self._data_file["/"]
        for q in node.split("/"):
            if q:
                myGroup = myGroup.require_group(q)
        arr = myGroup.create_dataset(name, data=array)
        if set_attr is not None:
            for i, j in set_attr.items():
                arr.attrs[i] = j
        self._data_file.flush()

    def _reload_data_file(self, *args, **kwargs):
        if self._data_file is None:
            return
        self._data_file.close()
        del self._data_file
        self._data_file = h5py.File(self.__data_filename, mode=self._data_mode)

    def get_data(self, node, name):
        """
        Return the dataset with a given *name* located at *node* in the
        datafile.
        """
        if self._data_file is None:
            return None
        if node[0] != "/":
            node = f"/{node}"

        myGroup = self._data_file["/"]
        for group in node.split("/"):
            if group:
                if group not in myGroup:
                    return None
                myGroup = myGroup[group]
        if name not in myGroup:
            return None

        full_name = f"{node}/{name}"
        try:
            return self._data_file[full_name][:]
        except TypeError:
            return self._data_file[full_name]

    def _get_particle_type_counts(self):
        # this is implemented by subclasses
        raise NotImplementedError

    def _close_data_file(self):
        if self._data_file:
            self._data_file.close()
            del self._data_file
            self._data_file = None

    def _split_fields(self, fields):
        # This will split fields into either generated or read fields
        fields_to_read, fields_to_generate = [], []
        for ftype, fname in fields:
            if fname in self.field_list or (ftype, fname) in self.field_list:
                fields_to_read.append((ftype, fname))
            elif (
                fname in self.ds.derived_field_list
                or (ftype, fname) in self.ds.derived_field_list
            ):
                fields_to_generate.append((ftype, fname))
            else:
                raise YTFieldNotFound((ftype, fname), self.ds)
        return fields_to_read, fields_to_generate

    def _read_particle_fields(self, fields, dobj, chunk=None):
        if len(fields) == 0:
            return {}, []
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        selector = dobj.selector
        if chunk is None:
            self._identify_base_chunk(dobj)
        chunks = self._chunk_io(dobj, cache=False)
        fields_to_return = self.io._read_particle_selection(
            chunks, selector, fields_to_read
        )
        return fields_to_return, fields_to_generate

    def _read_fluid_fields(self, fields, dobj, chunk=None):
        if len(fields) == 0:
            return {}, []
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        selector = dobj.selector
        if chunk is None:
            self._identify_base_chunk(dobj)
            chunk_size = dobj.size
        else:
            chunk_size = chunk.data_size
        fields_to_return = self.io._read_fluid_selection(
            self._chunk_io(dobj), selector, fields_to_read, chunk_size
        )
        return fields_to_return, fields_to_generate

    def _chunk(self, dobj, chunking_style, ngz=0, **kwargs):
        # A chunk is either None or (grids, size)
        if dobj._current_chunk is None:
            self._identify_base_chunk(dobj)
        if ngz != 0 and chunking_style != "spatial":
            raise NotImplementedError
        if chunking_style == "all":
            return self._chunk_all(dobj, **kwargs)
        elif chunking_style == "spatial":
            return self._chunk_spatial(dobj, ngz, **kwargs)
        elif chunking_style == "io":
            return self._chunk_io(dobj, **kwargs)
        else:
            raise NotImplementedError


def cacheable_property(func):
    # not quite equivalent to functools.cached_property
    # this decorator allows cached to be disabled via a self._cache flag attribute
    n = f"_{func.__name__}"

    @property
    def cacheable_func(self):
        if self._cache and getattr(self, n, None) is not None:
            return getattr(self, n)
        if self.data_size is None:
            tr = self._accumulate_values(n[1:])
        else:
            tr = func(self)
        if self._cache:
            setattr(self, n, tr)
        return tr

    return cacheable_func


class YTDataChunk:
    def __init__(
        self,
        dobj,
        chunk_type,
        objs,
        data_size=None,
        field_type=None,
        cache=False,
        fast_index=None,
    ):
        self.dobj = dobj
        self.chunk_type = chunk_type
        self.objs = objs
        self.data_size = data_size
        self._field_type = field_type
        self._cache = cache
        self._fast_index = fast_index

    def _accumulate_values(self, method):
        # We call this generically.  It's somewhat slower, since we're doing
        # costly getattr functions, but this allows us to generalize.
        mname = f"select_{method}"
        arrs = []
        for obj in self._fast_index or self.objs:
            f = getattr(obj, mname)
            arrs.append(f(self.dobj))
        if method == "dtcoords":
            arrs = [arr[0] for arr in arrs]
        elif method == "tcoords":
            arrs = [arr[1] for arr in arrs]
        arrs = uconcatenate(arrs)
        self.data_size = arrs.shape[0]
        return arrs

    @cacheable_property
    def fcoords(self):
        if self._fast_index is not None:
            ci = self._fast_index.select_fcoords(self.dobj.selector, self.data_size)
            ci = YTArray(ci, units="code_length", registry=self.dobj.ds.unit_registry)
            return ci
        ci = np.empty((self.data_size, 3), dtype="float64")
        ci = YTArray(ci, units="code_length", registry=self.dobj.ds.unit_registry)
        if self.data_size == 0:
            return ci
        ind = 0
        for obj in self._fast_index or self.objs:
            c = obj.select_fcoords(self.dobj)
            if c.shape[0] == 0:
                continue
            ci.d[ind : ind + c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cacheable_property
    def icoords(self):
        if self._fast_index is not None:
            ci = self._fast_index.select_icoords(self.dobj.selector, self.data_size)
            return ci
        ci = np.empty((self.data_size, 3), dtype="int64")
        if self.data_size == 0:
            return ci
        ind = 0
        for obj in self._fast_index or self.objs:
            c = obj.select_icoords(self.dobj)
            if c.shape[0] == 0:
                continue
            ci[ind : ind + c.shape[0], :] = c
            ind += c.shape[0]
        return ci

    @cacheable_property
    def fwidth(self):
        if self._fast_index is not None:
            ci = self._fast_index.select_fwidth(self.dobj.selector, self.data_size)
            ci = YTArray(ci, units="code_length", registry=self.dobj.ds.unit_registry)
            return ci
        ci = np.empty((self.data_size, 3), dtype="float64")
        ci = YTArray(ci, units="code_length", registry=self.dobj.ds.unit_registry)
        if self.data_size == 0:
            return ci
        ind = 0
        for obj in self._fast_index or self.objs:
            c = obj.select_fwidth(self.dobj)
            if c.shape[0] == 0:
                continue
            ci.d[ind : ind + c.shape[0], :] = c

            ind += c.shape[0]
        return ci

    @cacheable_property
    def ires(self):
        if self._fast_index is not None:
            ci = self._fast_index.select_ires(self.dobj.selector, self.data_size)
            return ci
        ci = np.empty(self.data_size, dtype="int64")
        if self.data_size == 0:
            return ci
        ind = 0
        for obj in self._fast_index or self.objs:
            c = obj.select_ires(self.dobj)
            if c.shape == 0:
                continue
            ci[ind : ind + c.size] = c
            ind += c.size
        return ci

    @cacheable_property
    def tcoords(self):
        self.dtcoords
        return self._tcoords

    @cacheable_property
    def dtcoords(self):
        ct = np.empty(self.data_size, dtype="float64")
        cdt = np.empty(self.data_size, dtype="float64")
        self._tcoords = ct  # Se this for tcoords
        if self.data_size == 0:
            return cdt
        ind = 0
        for obj in self._fast_index or self.objs:
            gdt, gt = obj.select_tcoords(self.dobj)
            if gt.size == 0:
                continue
            ct[ind : ind + gt.size] = gt
            cdt[ind : ind + gdt.size] = gdt
            ind += gt.size
        return cdt

    @cacheable_property
    def fcoords_vertex(self):
        nodes_per_elem = self.dobj.index.meshes[0].connectivity_indices.shape[1]
        dim = self.dobj.ds.dimensionality
        ci = np.empty((self.data_size, nodes_per_elem, dim), dtype="float64")
        ci = YTArray(ci, units="code_length", registry=self.dobj.ds.unit_registry)
        if self.data_size == 0:
            return ci
        ind = 0
        for obj in self.objs:
            c = obj.select_fcoords_vertex(self.dobj)
            if c.shape[0] == 0:
                continue
            ci.d[ind : ind + c.shape[0], :, :] = c
            ind += c.shape[0]
        return ci


class ChunkDataCache:
    def __init__(self, base_iter, preload_fields, geometry_handler, max_length=256):
        # At some point, max_length should instead become a heuristic function,
        # potentially looking at estimated memory usage.  Note that this never
        # initializes the iterator; it assumes the iterator is already created,
        # and it calls next() on it.
        self.base_iter = base_iter.__iter__()
        self.queue = []
        self.max_length = max_length
        self.preload_fields = preload_fields
        self.geometry_handler = geometry_handler
        self.cache = {}

    def __iter__(self):
        return self

    def __next__(self):
        if len(self.queue) == 0:
            for _ in range(self.max_length):
                try:
                    self.queue.append(next(self.base_iter))
                except StopIteration:
                    break
            # If it's still zero ...
            if len(self.queue) == 0:
                raise StopIteration
            chunk = YTDataChunk(None, "cache", self.queue, cache=False)
            self.cache = (
                self.geometry_handler.io._read_chunk_data(chunk, self.preload_fields)
                or {}
            )
        g = self.queue.pop(0)
        g._initialize_cache(self.cache.pop(g.id, {}))
        return g


def is_curvilinear(geo):
    # tell geometry is curvilinear or not
    if geo in ["polar", "cylindrical", "spherical"]:
        return True
    else:
        return False
