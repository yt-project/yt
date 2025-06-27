import os
import typing
import weakref
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.logger import ytLogger as mylog

from .fields import ChollaFieldInfo

# this is a hacky workaround to get _h5py.File to work in annotations. We can probably
# address this issue more robustly by directly modifying yt.utilities.on_demand_imports
if typing.TYPE_CHECKING:
    import h5py as _h5py
else:
    from yt.utilities.on_demand_imports import _h5py


@dataclass(kw_only=True, slots=True, frozen=True)
class _BlockDiskMapping:
    """Contains info for mapping blockids to locations in hdf5 files"""

    # ``fname_template.format(blockid=...)`` produces the file containing blockid (this
    # can properly handle cases where all blocks are stored in a single file)
    fname_template: str
    # group containing field data (empty string denotes the root group)
    field_group: str
    # maps blockid to an index that select all associated data from a field-dataset
    field_idx_map: Mapping[int, tuple[int | slice, ...]]

    # in the future, we may add particle_group and particle_idx_map


def _infer_blockid_location_arr(fname_template, global_dims, arr_shape):
    # used when hdf5 files don't have an explicit "domain" group
    blockid_location_arr = np.empty(shape=tuple(int(e) for e in arr_shape), dtype="i8")
    if blockid_location_arr.size == 1:
        # primarily intended to handle the result of older concatenation scripts (it
        # also handles the case when only a single block is used, which is okay)
        blockid_location_arr[0, 0, 0] = 0
    else:  # handle distributed cholla datasets
        local_dims, rem = np.divmod(global_dims, blockid_location_arr.shape)
        assert np.all(rem == 0) and np.all(local_dims > 0)
        for blockid in range(0, blockid_location_arr.size):
            with _h5py.File(fname_template.format(blockid=blockid), "r") as f:
                tmp, rem = np.divmod(f.attrs["offset"][:], local_dims)
            assert np.all(rem == 0)  # sanity check
            idx3D = tuple(int(e) for e in tmp)
            blockid_location_arr[idx3D] = blockid
    return blockid_location_arr


def _determine_data_layout(f: _h5py.File) -> tuple[np.ndarray, _BlockDiskMapping]:
    """Determine the data layout of the snapshot

    The premise is that the basic different data formats shouldn't
    matter outside of this function."""
    filename = f.filename

    # STEP 1: infer the template for all Cholla data-files by inspecting filename
    # ===========================================================================
    # There are 2 conventions for the names of Cholla's data-files:
    #  1. "root.h5.{blockid}" is the standard format Cholla uses when writing files
    #     storing a single snapshot. Each MPI-rank will write a separate file and
    #     replace ``{blockid}`` with MPI-rank (Modern Cholla versions without MPI
    #     replace ``{blockid}`` with ``0``)
    #  2. "root.h5": is the standard format used by Cholla's concatenation scripts
    #     (older versions of Cholla without MPI also used this format to name outputs)
    _dir, _base = os.path.split(filename)
    _sep_i = _base.rfind(".")
    no_suffix = (_sep_i == -1) or (_base[_sep_i:] == "") or (_base[:_sep_i] == "")
    if no_suffix or not _base[_sep_i + 1 :].isdecimal():
        inferred_fname_template = filename  # filename doesn't change based on blockid
        cur_filename_suffix = None
    else:
        inferred_fname_template = os.path.join(_dir, _base[:_sep_i]) + ".{blockid}"
        cur_filename_suffix = int(_base[_sep_i + 1 :])

    # STEP 2: Check whether the hdf5 file has a flat structure
    # ========================================================
    # Historically, we would always store datasets directly in the root group of the
    # data file. More recent concatenation scripts store no data in groups.
    flat_structure = any(not isinstance(elem, _h5py.Group) for elem in f.values())

    # STEP 3: Extract basic domain info information from the file(s)
    # ==============================================================
    has_explicit_domain_info = "domain" in f
    if has_explicit_domain_info:
        # this branch primarily handles concatenated files made with newer logic
        blockid_location_arr = f["domain/blockid_location_arr"][...]
        field_idx_map = {
            int(blockid): (i, slice(None), slice(None), slice(None))
            for i, blockid in enumerate(f["domain/stored_blockid_list"][...])
        }
        consolidated_data = len(field_idx_map) == blockid_location_arr.size
        if not consolidated_data:
            # in the near future, we may support one of the 2 cases:
            # > if (flat_structure):
            # >     _common_idx = (slice(None), slice(None), slice(None))
            # > else:
            # >     _common_idx = (0, slice(None), slice(None), slice(None))
            # > field_idx_map = defaultdict(lambda arg=_common_idx: arg)
            raise ValueError(
                "no support for reading Cholla datasets where data is distributed "
                "among files that explicitly encode domain info."
            )
    else:  # (not has_explicit_domain_info)
        # this branch covers distributed datasets (directly written by Cholla) and
        # older concatenated files.
        #
        # historically, when the dataset is concatenated (in post-processing),
        # the "nprocs" hdf5 attribute has been dropped
        blockid_location_arr = _infer_blockid_location_arr(
            fname_template=inferred_fname_template,
            global_dims=f.attrs["dims"].astype("=i8"),
            arr_shape=f.attrs.get("nprocs", np.array([1, 1, 1])).astype("=i8"),
        )
        consolidated_data = blockid_location_arr.size == 1

        def _get_common_idx():
            return (slice(None), slice(None), slice(None))

        field_idx_map = defaultdict(_get_common_idx)

    # STEP 4: Finalize the fname template
    # ===================================
    if consolidated_data:
        fname_template = filename
    elif cur_filename_suffix != 0:
        raise ValueError(  # mostly just a sanity check!
            "filename passed to yt.load for a distributed cholla dataset must "
            "end in '.0'"
        )
    else:
        fname_template = inferred_fname_template

    mapping = _BlockDiskMapping(
        fname_template=fname_template,
        field_group="" if flat_structure else "field",
        field_idx_map=field_idx_map,
    )
    return blockid_location_arr, mapping


def _split_fname_procid_suffix(filename: str):
    """Splits ``filename`` at the '.' separating the beginning part of the
    string from the process-id suffix, and returns both parts in a 2-tuple.

    There are 2 conventions for the names of Cholla's data-files:

    1. "root.h5.{procid}" is the standard format used when Cholla directly
       writes out data-files. When outputting a snaphsot, each MPI-process
       writes a separate file, named using this template (``{procid}`` is
       replaced with MPI-rank). Modern versions of Cholla without MPI replace
       ``{procid}`` with ``0``. When passed a name of this form, the function
       returns ``("root.h5", "{procid}")``
    2. "root.h5": is the standard format used by Cholla's concatenation scirpts
       (older versions of the code compiled without MPI also used this format
       to name outputs). When passed a name of this form, the function returns
       ``("root.h5", "")``

    Notes
    -----
    The following examples illustrate how the behavior differs from that of
    ``os.path.splitext``

    - For "root.h5.3": ``os.path.splitext`` returns ``("root.h5", ".3")``,
      while this function returns ``("root.h5", "3")``
    - For "root.h5": os.path.splitext`` returns ``("root", ".h5")``, while this
      function returns ``("root.h5", "")``
    """

    # at this time, we expect the suffix to be the minimum number of characters
    # that are necessary to represent the process id. For flexibility, we will
    # allow extra zero-padding

    match filename.rpartition("."):
        case ("", "", _) | (_, ".", ""):
            return (filename, "")
        case (prefix, ".", _) if prefix == "" or prefix[-1] == "/":
            raise ValueError(
                f"can't split a process-suffix off of {filename!r} "
                "since the remaining filename would be empty"
            )
        case (prefix, ".", suffix):
            return (prefix, suffix) if suffix.isdecimal() else (filename, "")


class ChollaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dims, filename):
        super().__init__(id, filename=filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims


class ChollaHierarchy(GridIndex):
    grid = ChollaGrid
    _grid_chunksize = 1

    def __init__(self, ds, dataset_type="cholla"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        with _h5py.File(self.index_filename, mode="r") as h5f:
            grp = h5f.get("field", h5f)
            self.field_list = [("cholla", k) for k in grp.keys()]

    def _count_grids(self):
        with _h5py.File(self.index_filename, "r") as f:
            self._blockid_location_arr, self._block_mapping = _determine_data_layout(f)
        self.num_grids = self._blockid_location_arr.size

    def _parse_index(self):
        self.grids = np.empty(self.num_grids, dtype="object")

        shape_arr = np.array(self._blockid_location_arr.shape)
        dims_local = (self.ds.domain_dimensions[:] / shape_arr).astype("=i8")

        for idx3D, blockid in np.ndenumerate(self._blockid_location_arr):
            idx3D_arr = np.array(idx3D)
            left_frac = idx3D_arr / shape_arr
            right_frac = (1 + idx3D_arr) / shape_arr

            level = 0

            self.grids[blockid] = self.grid(
                blockid,
                index=self,
                level=level,
                dims=dims_local,
                filename=self._block_mapping.fname_template.format(blockid=blockid),
            )

            self.grid_left_edge[blockid, :] = left_frac
            self.grid_right_edge[blockid, :] = right_frac
            self.grid_dimensions[blockid, :] = dims_local
            self.grid_levels[blockid, 0] = level
            self.grid_particle_count[blockid, 0] = 0

        slope = self.ds.domain_width / self.ds.arr(np.ones(3), "code_length")
        self.grid_left_edge = self.grid_left_edge * slope + self.ds.domain_left_edge
        self.grid_right_edge = self.grid_right_edge * slope + self.ds.domain_left_edge

        self.max_level = 0

    def _populate_grid_objects(self):
        for i in range(self.num_grids):
            g = self.grids[i]
            g._prepare_grid()
            g._setup_dx()


class ChollaDataset(Dataset):
    _load_requirements = ["h5py"]
    _index_class = ChollaHierarchy
    _field_info_class = ChollaFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="cholla",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.fluid_types += ("cholla",)
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the defaults, but if they are listed
        # in the HDF5 attributes for a file, which is loaded first, then those are
        # used instead.
        #
        if not self.length_unit:
            self.length_unit = self.quan(1.0, "pc")
        if not self.mass_unit:
            self.mass_unit = self.quan(1.0, "Msun")
        if not self.time_unit:
            self.time_unit = self.quan(1000, "yr")
        if not self.velocity_unit:
            self.velocity_unit = self.quan(1.0, "cm/s")
        if not self.magnetic_unit:
            self.magnetic_unit = self.quan(1.0, "gauss")

        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):
        with _h5py.File(self.parameter_filename, mode="r") as h5f:
            attrs = h5f.attrs
            self.parameters = dict(attrs.items())
            self.domain_left_edge = attrs["bounds"][:].astype("=f8")
            self.domain_right_edge = self.domain_left_edge + attrs["domain"][:].astype(
                "=f8"
            )
            self.dimensionality = len(attrs["dims"][:])
            self.domain_dimensions = attrs["dims"][:].astype("=i8")
            self.current_time = attrs["t"][:]
            self._periodicity = tuple(attrs.get("periodicity", (False, False, False)))
            self.gamma = attrs.get("gamma", 5.0 / 3.0)
            if (self.default_species_fields is not None) and "mu" in attrs:
                raise ValueError(
                    'default_species_fields must be None when "mu" is an hdf5 attribute'
                )
            elif "mu" in attrs:
                self.mu = attrs["mu"]
            elif self.default_species_fields is None:
                # other yt-machinery can't handle ds.mu == None, so we simply
                # avoid defining the mu attribute if we don't know its value
                mylog.info(
                    'add the "mu" hdf5 attribute OR use the default_species_fields kwarg '
                    "to compute temperature"
                )
            self.refine_by = 1

            # If header specifies code units, default to those (in CGS)
            length_unit = attrs.get("length_unit", None)
            mass_unit = attrs.get("mass_unit", None)
            time_unit = attrs.get("time_unit", None)
            velocity_unit = attrs.get("velocity_unit", None)
            magnetic_unit = attrs.get("magnetic_unit", None)
            if length_unit:
                self.length_unit = self.quan(length_unit[0], "cm")
            if mass_unit:
                self.mass_unit = self.quan(mass_unit[0], "g")
            if time_unit:
                self.time_unit = self.quan(time_unit[0], "s")
            if velocity_unit:
                self.velocity_unit = self.quan(velocity_unit[0], "cm/s")
            if magnetic_unit:
                self.magnetic_unit = self.quan(magnetic_unit[0], "gauss")

            # this minimalistic implementation fills the requirements for
            # this frontend to run, change it to make it run _correctly_ !
            for key, unit in self.__class__.default_units.items():
                setdefaultattr(self, key, self.quan(1, unit))

        # CHOLLA cannot yet be run as a cosmological simulation
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

        # CHOLLA datasets are always unigrid cartesian
        self.geometry = Geometry.CARTESIAN

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        if cls._missing_load_requirements():
            return False

        try:
            fileh = _h5py.File(filename, mode="r")
        except OSError:
            return False

        try:
            attrs = fileh.attrs
        except AttributeError:
            return False
        else:
            return (
                "bounds" in attrs
                and "domain" in attrs
                and attrs.get("data_type") != "yt_light_ray"
            )
        finally:
            fileh.close()
