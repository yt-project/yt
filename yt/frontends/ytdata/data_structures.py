import os
import sys
import weakref
from collections import defaultdict
from numbers import Number as numeric_type
from typing import Tuple, Type

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.particle_unions import ParticleUnion
from yt.data_objects.profiles import (
    Profile1DFromDataset,
    Profile2DFromDataset,
    Profile3DFromDataset,
)
from yt.data_objects.static_output import Dataset, ParticleFile, validate_index_order
from yt.fields.field_exceptions import NeedsGridType
from yt.fields.field_info_container import FieldInfoContainer
from yt.funcs import is_root, parse_h5_attr
from yt.geometry.geometry_handler import Index
from yt.geometry.grid_geometry_handler import GridIndex
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.units import dimensions
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.units.unit_registry import UnitRegistry  # type: ignore
from yt.units.yt_array import YTQuantity
from yt.utilities.exceptions import GenerationInProgress, YTFieldTypeNotFound
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_root_only
from yt.utilities.tree_container import TreeContainer

from .fields import YTDataContainerFieldInfo, YTGridFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property
_grid_data_containers = ["arbitrary_grid", "covering_grid", "smoothed_covering_grid"]
_set_attrs = {"periodicity": "_periodicity"}


class SavedDataset(Dataset):
    """
    Base dataset class for products of calling save_as_dataset.
    """

    _con_attrs: Tuple[str, ...] = ()

    def _parse_parameter_file(self):
        self.refine_by = 2
        with h5py.File(self.parameter_filename, mode="r") as f:
            for key in f.attrs.keys():
                v = parse_h5_attr(f, key)
                if key == "con_args":
                    try:
                        v = eval(v)
                    except ValueError:
                        # support older ytdata outputs
                        v = v.astype("str")
                    except NameError:
                        # This is the most common error we expect, and it
                        # results from having the eval do a concatenated decoded
                        # set of the values.
                        v = [_.decode("utf8") for _ in v]
                self.parameters[key] = v
            self._with_parameter_file_open(f)

        # if saved, restore unit registry from the json string
        if "unit_registry_json" in self.parameters:
            self.unit_registry = UnitRegistry.from_json(
                self.parameters["unit_registry_json"]
            )
            # reset self.arr and self.quan to use new unit_registry
            self._arr = None
            self._quan = None
            for dim in [
                "length",
                "mass",
                "pressure",
                "temperature",
                "time",
                "velocity",
            ]:
                cu = "code_" + dim
                if cu not in self.unit_registry:
                    self.unit_registry.add(cu, 1.0, getattr(dimensions, dim))
            if "code_magnetic" not in self.unit_registry:
                self.unit_registry.add(
                    "code_magnetic", 0.1**0.5, dimensions.magnetic_field_cgs
                )

        # if saved, set unit system
        if "unit_system_name" in self.parameters:
            unit_system = self.parameters["unit_system_name"]
            del self.parameters["unit_system_name"]
        else:
            unit_system = "cgs"
        # reset unit system since we may have a new unit registry
        self._assign_unit_system(unit_system)

        # assign units to parameters that have associated unit string
        del_pars = []
        for par in self.parameters:
            ustr = f"{par}_units"
            if ustr in self.parameters:
                if isinstance(self.parameters[par], np.ndarray):
                    to_u = self.arr
                else:
                    to_u = self.quan
                self.parameters[par] = to_u(self.parameters[par], self.parameters[ustr])
                del_pars.append(ustr)
        for par in del_pars:
            del self.parameters[par]

        for attr in self._con_attrs:
            try:
                sattr = _set_attrs.get(attr, attr)
                setattr(self, sattr, self.parameters.get(attr))
            except TypeError:
                # some Dataset attributes are properties with setters
                # which may not accept None as an input
                pass

        if self.geometry is None:
            self.geometry = "cartesian"

    def _with_parameter_file_open(self, f):
        # This allows subclasses to access the parameter file
        # while it's still open to get additional information.
        pass

    def set_units(self):
        if "unit_registry_json" in self.parameters:
            self._set_code_unit_attributes()
            del self.parameters["unit_registry_json"]
        else:
            super().set_units()

    def _set_code_unit_attributes(self):
        attrs = (
            "length_unit",
            "mass_unit",
            "time_unit",
            "velocity_unit",
            "magnetic_unit",
        )
        cgs_units = ("cm", "g", "s", "cm/s", "gauss")
        base_units = np.ones(len(attrs))
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if attr in self.parameters and isinstance(
                self.parameters[attr], YTQuantity
            ):
                uq = self.parameters[attr]
            elif attr in self.parameters and f"{attr}_units" in self.parameters:
                uq = self.quan(self.parameters[attr], self.parameters[f"{attr}_units"])
                del self.parameters[attr]
                del self.parameters[f"{attr}_units"]
            elif isinstance(unit, str):
                uq = self.quan(1.0, unit)
            elif isinstance(unit, numeric_type):
                uq = self.quan(unit, cgs_unit)
            elif isinstance(unit, YTQuantity):
                uq = unit
            elif isinstance(unit, tuple):
                uq = self.quan(unit[0], unit[1])
            else:
                raise RuntimeError(f"{attr} ({unit}) is invalid.")
            setattr(self, attr, uq)


class YTDataset(SavedDataset):
    """Base dataset class for all ytdata datasets."""

    _con_attrs = (
        "cosmological_simulation",
        "current_time",
        "current_redshift",
        "hubble_constant",
        "omega_matter",
        "omega_lambda",
        "dimensionality",
        "domain_dimensions",
        "geometry",
        "periodicity",
        "domain_left_edge",
        "domain_right_edge",
        "container_type",
        "data_type",
    )

    def _with_parameter_file_open(self, f):
        self.num_particles = {
            group: parse_h5_attr(f[group], "num_elements")
            for group in f
            if group != self.default_fluid_type
        }

    def create_field_info(self):
        self.field_dependencies = {}
        self.derived_field_list = []
        self.filtered_particle_types = []
        self.field_info = self._field_info_class(self, self.field_list)
        self.coordinates.setup_fields(self.field_info)
        self.field_info.setup_fluid_fields()
        for ptype in self.particle_types:
            self.field_info.setup_particle_fields(ptype)

        self._setup_gas_alias()
        self.field_info.setup_fluid_index_fields()

        if "all" not in self.particle_types:
            mylog.debug("Creating Particle Union 'all'")
            pu = ParticleUnion("all", list(self.particle_types_raw))
            self.add_particle_union(pu)
        self.field_info.setup_extra_union_fields()
        mylog.debug("Loading field plugins.")
        self.field_info.load_all_plugins()
        deps, unloaded = self.field_info.check_derived_fields()
        self.field_dependencies.update(deps)

    def _setup_gas_alias(self):
        pass

    def _setup_override_fields(self):
        pass


class YTDataHDF5File(ParticleFile):
    def __init__(self, ds, io, filename, file_id, range):
        with h5py.File(filename, mode="r") as f:
            self.header = {field: parse_h5_attr(f, field) for field in f.attrs.keys()}

        super().__init__(ds, io, filename, file_id, range)


class YTDataContainerDataset(YTDataset):
    """Dataset for saved geometric data containers."""

    _index_class = ParticleIndex
    _file_class = YTDataHDF5File
    _field_info_class: Type[FieldInfoContainer] = YTDataContainerFieldInfo
    _suffix = ".h5"
    fluid_types = ("grid", "gas", "deposit", "index")

    def __init__(
        self,
        filename,
        dataset_type="ytdatacontainer_hdf5",
        index_order=None,
        index_filename=None,
        units_override=None,
        unit_system="cgs",
    ):
        self.index_order = validate_index_order(index_order)
        self.index_filename = index_filename
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        self.particle_types_raw = tuple(self.num_particles.keys())
        self.particle_types = self.particle_types_raw
        self.filename_template = self.parameter_filename
        self.file_count = 1
        self.domain_dimensions = np.ones(3, "int32")

    def _setup_gas_alias(self):
        "Alias the grid type to gas by making a particle union."

        if "grid" in self.particle_types and "gas" not in self.particle_types:
            pu = ParticleUnion("gas", ["grid"])
            self.add_particle_union(pu)
        # We have to alias this because particle unions only
        # cover the field_list.
        self.field_info.alias(("gas", "cell_volume"), ("grid", "cell_volume"))

    @cached_property
    def data(self):
        """
        Return a data container configured like the original used to
        create this dataset.
        """

        # Some data containers can't be reconstructed in the same way
        # since this is now particle-like data.
        data_type = self.parameters.get("data_type")
        container_type = self.parameters.get("container_type")
        ex_container_type = ["cutting", "quad_proj", "ray", "slice", "cut_region"]
        if data_type == "yt_light_ray" or container_type in ex_container_type:
            mylog.info("Returning an all_data data container.")
            return self.all_data()

        my_obj = getattr(self, self.parameters["container_type"])
        my_args = [self.parameters[con_arg] for con_arg in self.parameters["con_args"]]
        return my_obj(*my_args)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            cont_type = parse_h5_attr(f, "container_type")
            if data_type is None:
                return False
            if (
                data_type == "yt_data_container"
                and cont_type not in _grid_data_containers
            ):
                return True
        return False


class YTDataLightRayDataset(YTDataContainerDataset):
    """Dataset for saved LightRay objects."""

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        self._restore_light_ray_solution()

    def _restore_light_ray_solution(self):
        """
        Restore all information associate with the light ray solution
        to its original form.
        """
        key = "light_ray_solution"
        self.light_ray_solution = []
        lrs_fields = [
            par for par in self.parameters if key in par and not par.endswith("_units")
        ]
        if len(lrs_fields) == 0:
            return
        self.light_ray_solution = [{} for val in self.parameters[lrs_fields[0]]]
        for sp3 in ["unique_identifier", "filename"]:
            ksp3 = f"{key}_{sp3}"
            if ksp3 not in lrs_fields:
                continue
            self.parameters[ksp3] = self.parameters[ksp3].astype(str)
        for field in lrs_fields:
            field_name = field[len(key) + 1 :]
            for i in range(self.parameters[field].shape[0]):
                self.light_ray_solution[i][field_name] = self.parameters[field][i]

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            if data_type in ["yt_light_ray"]:
                return True
        return False


class YTSpatialPlotDataset(YTDataContainerDataset):
    """Dataset for saved slices and projections."""

    _field_info_class = YTGridFieldInfo

    def __init__(self, *args, **kwargs):
        super().__init__(*args, dataset_type="ytspatialplot_hdf5", **kwargs)

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        if self.parameters["container_type"] == "proj":
            if (
                isinstance(self.parameters["weight_field"], str)
                and self.parameters["weight_field"] == "None"
            ):
                self.parameters["weight_field"] = None
            elif isinstance(self.parameters["weight_field"], np.ndarray):
                self.parameters["weight_field"] = tuple(self.parameters["weight_field"])

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            cont_type = parse_h5_attr(f, "container_type")
            if data_type == "yt_data_container" and cont_type in [
                "cutting",
                "proj",
                "slice",
                "quad_proj",
            ]:
                return True
        return False


class YTGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, gid, index, filename=None):
        AMRGridPatch.__init__(self, gid, filename=filename, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = 0
        self.LeftEdge = self.index.ds.domain_left_edge
        self.RightEdge = self.index.ds.domain_right_edge

    def __getitem__(self, key):
        tr = super(AMRGridPatch, self).__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        finfo = self.ds._get_field_info(*fields[0])
        if not finfo.sampling_type == "particle":
            return tr.reshape(self.ActiveDimensions[: self.ds.dimensionality])
        return tr

    @property
    def Parent(self):
        return None

    @property
    def Children(self):
        return []


class YTDataHierarchy(GridIndex):
    def __init__(self, ds, dataset_type=None):
        self.dataset_type = dataset_type
        self.float_type = "float64"
        self.dataset = weakref.proxy(ds)
        self.directory = os.getcwd()
        super().__init__(ds, dataset_type)

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_dimensions[:] = self.ds.domain_dimensions
        self.grid_left_edge[:] = self.ds.domain_left_edge
        self.grid_right_edge[:] = self.ds.domain_right_edge
        self.grid_levels[:] = np.zeros(self.num_grids)
        self.grid_procs = np.zeros(self.num_grids)
        self.grid_particle_count[:] = sum(self.ds.num_particles.values())
        self.grids = []
        for gid in range(self.num_grids):
            self.grids.append(self.grid(gid, self))
            self.grids[gid].Level = self.grid_levels[gid, 0]
        self.max_level = self.grid_levels.max()
        temp_grids = np.empty(self.num_grids, dtype="object")
        for i, grid in enumerate(self.grids):
            grid.filename = self.ds.parameter_filename
            grid._prepare_grid()
            grid.proc_num = self.grid_procs[i]
            temp_grids[i] = grid
        self.grids = temp_grids

    def _detect_output_fields(self):
        self.field_list = []
        self.ds.field_units = self.ds.field_units or {}
        with h5py.File(self.ds.parameter_filename, mode="r") as f:
            for group in f:
                for field in f[group]:
                    field_name = (str(group), str(field))
                    self.field_list.append(field_name)
                    self.ds.field_units[field_name] = parse_h5_attr(
                        f[group][field], "units"
                    )


class YTGridHierarchy(YTDataHierarchy):
    grid = YTGrid

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
        self.max_level = self.grid_levels.max()


class YTGridDataset(YTDataset):
    """Dataset for saved covering grids, arbitrary grids, and FRBs."""

    _index_class: Type[Index] = YTGridHierarchy
    _field_info_class = YTGridFieldInfo
    _dataset_type = "ytgridhdf5"
    geometry = "cartesian"
    default_fluid_type = "grid"
    fluid_types: Tuple[str, ...] = ("grid", "gas", "deposit", "index")

    def __init__(self, filename, unit_system="cgs"):
        super().__init__(filename, self._dataset_type, unit_system=unit_system)
        self.data = self.index.grids[0]

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        self.num_particles.pop(self.default_fluid_type, None)
        self.particle_types_raw = tuple(self.num_particles.keys())
        self.particle_types = self.particle_types_raw

        # correct domain dimensions for the covering grid dimension
        self.base_domain_left_edge = self.domain_left_edge
        self.base_domain_right_edge = self.domain_right_edge
        self.base_domain_dimensions = self.domain_dimensions

        if self.container_type in _grid_data_containers:
            self.domain_left_edge = self.parameters["left_edge"]

            if "level" in self.parameters["con_args"]:
                dx = (self.base_domain_right_edge - self.base_domain_left_edge) / (
                    self.domain_dimensions * self.refine_by ** self.parameters["level"]
                )
                self.domain_right_edge = (
                    self.domain_left_edge + self.parameters["ActiveDimensions"] * dx
                )
                self.domain_dimensions = (
                    (self.domain_right_edge - self.domain_left_edge) / dx
                ).astype(int)
            else:
                self.domain_right_edge = self.parameters["right_edge"]
                self.domain_dimensions = self.parameters["ActiveDimensions"]
                dx = (
                    self.domain_right_edge - self.domain_left_edge
                ) / self.domain_dimensions

            periodicity = (
                np.abs(self.domain_left_edge - self.base_domain_left_edge) < 0.5 * dx
            )
            periodicity &= (
                np.abs(self.domain_right_edge - self.base_domain_right_edge) < 0.5 * dx
            )
            self._periodicity = periodicity

        elif self.data_type == "yt_frb":
            dle = self.domain_left_edge
            self.domain_left_edge = uconcatenate(
                [self.parameters["left_edge"].to(dle.units), [0] * dle.uq]
            )
            dre = self.domain_right_edge
            self.domain_right_edge = uconcatenate(
                [self.parameters["right_edge"].to(dre.units), [1] * dre.uq]
            )
            self.domain_dimensions = np.concatenate(
                [self.parameters["ActiveDimensions"], [1]]
            )

    def _setup_gas_alias(self):
        "Alias the grid type to gas with a field alias."

        for ftype, field in self.field_list:
            if ftype == "grid":
                self.field_info.alias(("gas", field), ("grid", field))

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            cont_type = parse_h5_attr(f, "container_type")
            if data_type == "yt_frb":
                return True
            if data_type == "yt_data_container" and cont_type in _grid_data_containers:
                return True
        return False


class YTNonspatialGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, gid, index, filename=None):
        super().__init__(gid, filename=filename, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = 0
        self.LeftEdge = self.index.ds.domain_left_edge
        self.RightEdge = self.index.ds.domain_right_edge

    def __getitem__(self, key):
        tr = super(AMRGridPatch, self).__getitem__(key)
        try:
            fields = self._determine_fields(key)
        except YTFieldTypeNotFound:
            return tr
        self.ds._get_field_info(*fields[0])
        return tr

    def get_data(self, fields=None):
        if fields is None:
            return
        nfields = []
        apply_fields = defaultdict(list)
        for field in self._determine_fields(fields):
            if field[0] in self.ds.filtered_particle_types:
                f = self.ds.known_filters[field[0]]
                apply_fields[field[0]].append((f.filtered_type, field[1]))
            else:
                nfields.append(field)
        for filter_type in apply_fields:
            f = self.ds.known_filters[filter_type]
            with f.apply(self):
                self.get_data(apply_fields[filter_type])
        fields = nfields
        if len(fields) == 0:
            return
        # Now we collect all our fields
        # Here is where we need to perform a validation step, so that if we
        # have a field requested that we actually *can't* yet get, we put it
        # off until the end.  This prevents double-reading fields that will
        # need to be used in spatial fields later on.
        fields_to_get = []
        # This will be pre-populated with spatial fields
        fields_to_generate = []
        for field in self._determine_fields(fields):
            if field in self.field_data:
                continue
            finfo = self.ds._get_field_info(*field)
            try:
                finfo.check_available(self)
            except NeedsGridType:
                fields_to_generate.append(field)
                continue
            fields_to_get.append(field)
        if len(fields_to_get) == 0 and len(fields_to_generate) == 0:
            return
        elif self._locked:
            raise GenerationInProgress(fields)
        # Track which ones we want in the end
        ofields = set(list(self.field_data.keys()) + fields_to_get + fields_to_generate)
        # At this point, we want to figure out *all* our dependencies.
        fields_to_get = self._identify_dependencies(fields_to_get, self._spatial)
        # We now split up into readers for the types of fields
        fluids, particles = [], []
        finfos = {}
        for ftype, fname in fields_to_get:
            finfo = self.ds._get_field_info(ftype, fname)
            finfos[ftype, fname] = finfo
            if finfo.sampling_type == "particle":
                particles.append((ftype, fname))
            elif (ftype, fname) not in fluids:
                fluids.append((ftype, fname))

        # The _read method will figure out which fields it needs to get from
        # disk, and return a dict of those fields along with the fields that
        # need to be generated.
        read_fluids, gen_fluids = self.index._read_fluid_fields(
            fluids, self, self._current_chunk
        )
        for f, v in read_fluids.items():
            convert = True
            if v.dtype != np.float64:
                if finfos[f].units == "":
                    self.field_data[f] = v
                    convert = False
                else:
                    v = v.astype(np.float64)
            if convert:
                self.field_data[f] = self.ds.arr(v, units=finfos[f].units)
                self.field_data[f].convert_to_units(finfos[f].output_units)

        read_particles, gen_particles = self.index._read_fluid_fields(
            particles, self, self._current_chunk
        )
        for f, v in read_particles.items():
            convert = True
            if v.dtype != np.float64:
                if finfos[f].units == "":
                    self.field_data[f] = v
                    convert = False
                else:
                    v = v.astype(np.float64)
            if convert:
                self.field_data[f] = self.ds.arr(v, units=finfos[f].units)
                self.field_data[f].convert_to_units(finfos[f].output_units)

        fields_to_generate += gen_fluids + gen_particles
        self._generate_fields(fields_to_generate)
        for field in list(self.field_data.keys()):
            if field not in ofields:
                self.field_data.pop(field)

    @property
    def Parent(self):
        return None

    @property
    def Children(self):
        return []


class YTNonspatialHierarchy(YTDataHierarchy):
    grid = YTNonspatialGrid

    def _populate_grid_objects(self):
        for g in self.grids:
            g._setup_dx()
            # this is non-spatial, so remove the code_length units
            g.dds = self.ds.arr(g.dds.d, "")
            g.ActiveDimensions = self.ds.domain_dimensions
        self.max_level = self.grid_levels.max()

    def _read_fluid_fields(self, fields, dobj, chunk=None):
        if len(fields) == 0:
            return {}, []
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        selector = dobj.selector
        fields_to_return = self.io._read_fluid_selection(dobj, selector, fields_to_read)
        return fields_to_return, fields_to_generate


class YTNonspatialDataset(YTGridDataset):
    """Dataset for general array data."""

    _index_class = YTNonspatialHierarchy
    _field_info_class = YTGridFieldInfo
    _dataset_type = "ytnonspatialhdf5"
    geometry = "cartesian"
    default_fluid_type = "data"
    fluid_types: Tuple[str, ...] = ("data", "gas")

    def _parse_parameter_file(self):
        super(YTGridDataset, self)._parse_parameter_file()
        self.num_particles.pop(self.default_fluid_type, None)
        self.particle_types_raw = tuple(self.num_particles.keys())
        self.particle_types = self.particle_types_raw

    def _set_derived_attrs(self):
        # set some defaults just to make things go
        default_attrs = {
            "dimensionality": 3,
            "domain_dimensions": np.ones(3, dtype="int64"),
            "domain_left_edge": np.zeros(3),
            "domain_right_edge": np.ones(3),
            "_periodicity": np.ones(3, dtype="bool"),
        }
        for att, val in default_attrs.items():
            if getattr(self, att, None) is None:
                setattr(self, att, val)

    def _setup_classes(self):
        # We don't allow geometric selection for non-spatial datasets
        self.objects = []

    @parallel_root_only
    def print_key_parameters(self):
        for a in [
            "current_time",
            "domain_dimensions",
            "domain_left_edge",
            "domain_right_edge",
            "cosmological_simulation",
        ]:
            v = getattr(self, a)
            if v is not None:
                mylog.info("Parameters: %-25s = %s", a, v)
        if hasattr(self, "cosmological_simulation") and self.cosmological_simulation:
            for a in [
                "current_redshift",
                "omega_lambda",
                "omega_matter",
                "hubble_constant",
            ]:
                v = getattr(self, a)
                if v is not None:
                    mylog.info("Parameters: %-25s = %s", a, v)
        mylog.warning("Geometric data selection not available for this dataset type.")

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            if data_type == "yt_array_data":
                return True
        return False


class YTProfileDataset(YTNonspatialDataset):
    """Dataset for saved profile objects."""

    fluid_types = ("data", "gas", "standard_deviation")

    def __init__(self, filename, unit_system="cgs"):
        super().__init__(filename, unit_system=unit_system)

    _profile = None

    @property
    def profile(self):
        if self._profile is not None:
            return self._profile
        if self.dimensionality == 1:
            self._profile = Profile1DFromDataset(self)
        elif self.dimensionality == 2:
            self._profile = Profile2DFromDataset(self)
        elif self.dimensionality == 3:
            self._profile = Profile3DFromDataset(self)
        else:
            self._profile = None
        return self._profile

    def _parse_parameter_file(self):
        super(YTGridDataset, self)._parse_parameter_file()

        if (
            isinstance(self.parameters["weight_field"], str)
            and self.parameters["weight_field"] == "None"
        ):
            self.parameters["weight_field"] = None
        elif isinstance(self.parameters["weight_field"], np.ndarray):
            self.parameters["weight_field"] = tuple(
                self.parameters["weight_field"].astype(str)
            )

        for a in ["profile_dimensions"] + [
            f"{ax}_{attr}" for ax in "xyz"[: self.dimensionality] for attr in ["log"]
        ]:
            setattr(self, a, self.parameters[a])

        self.base_domain_left_edge = self.domain_left_edge
        self.base_domain_right_edge = self.domain_right_edge
        self.base_domain_dimensions = self.domain_dimensions

        domain_dimensions = np.ones(3, dtype="int64")
        domain_dimensions[: self.dimensionality] = self.profile_dimensions
        self.domain_dimensions = domain_dimensions
        domain_left_edge = np.zeros(3)
        domain_right_edge = np.ones(3)
        for i, ax in enumerate("xyz"[: self.dimensionality]):
            range_name = f"{ax}_range"
            my_range = self.parameters[range_name]
            if getattr(self, f"{ax}_log", False):
                my_range = np.log10(my_range)
            domain_left_edge[i] = my_range[0]
            domain_right_edge[i] = my_range[1]
            setattr(self, range_name, self.parameters[range_name])

            bin_field = f"{ax}_field"
            if (
                isinstance(self.parameters[bin_field], str)
                and self.parameters[bin_field] == "None"
            ):
                self.parameters[bin_field] = None
            elif isinstance(self.parameters[bin_field], np.ndarray):
                self.parameters[bin_field] = tuple(
                    ["data", self.parameters[bin_field].astype(str)[1]]
                )
            setattr(self, bin_field, self.parameters[bin_field])
        self.domain_left_edge = domain_left_edge
        self.domain_right_edge = domain_right_edge

    def _setup_gas_alias(self):
        "Alias the grid type to gas with a field alias."
        for ftype, field in self.field_list:
            if ftype == "data":
                self.field_info.alias(("gas", field), (ftype, field))

    def create_field_info(self):
        super().create_field_info()
        if self.parameters["weight_field"] is not None:
            self.field_info.alias(
                self.parameters["weight_field"], (self.default_fluid_type, "weight")
            )

    def _set_derived_attrs(self):
        self.domain_center = 0.5 * (self.domain_right_edge + self.domain_left_edge)
        self.domain_width = self.domain_right_edge - self.domain_left_edge

    def print_key_parameters(self):
        if is_root():
            mylog.info("YTProfileDataset")
            for a in ["dimensionality", "profile_dimensions"] + [
                f"{ax}_{attr}"
                for ax in "xyz"[: self.dimensionality]
                for attr in ["field", "range", "log"]
            ]:
                v = getattr(self, a)
                mylog.info("Parameters: %-25s = %s", a, v)
        super().print_key_parameters()

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            if data_type == "yt_profile":
                return True
        return False


class YTClumpContainer(TreeContainer):
    def __init__(
        self, clump_id, global_id, parent_id, contour_key, contour_id, ds=None
    ):
        self.clump_id = clump_id
        self.global_id = global_id
        self.parent_id = parent_id
        self.contour_key = contour_key
        self.contour_id = contour_id
        self.parent = None
        self.ds = ds
        TreeContainer.__init__(self)

    def add_child(self, child):
        if self.children is None:
            self.children = []
        self.children.append(child)
        child.parent = self

    def __repr__(self):
        return "Clump[%d]" % self.clump_id

    def __getitem__(self, field):
        g = self.ds.data
        f = g._determine_fields(field)[0]
        if f[0] == "clump":
            return g[f][self.global_id]
        if self.contour_id == -1:
            return g[f]
        cfield = (f[0], f"contours_{self.contour_key.decode('utf-8')}")
        if f[0] == "grid":
            return g[f][g[cfield] == self.contour_id]
        return self.parent[f][g[cfield] == self.contour_id]


class YTClumpTreeDataset(YTNonspatialDataset):
    """Dataset for saved clump-finder data."""

    def __init__(self, filename, unit_system="cgs"):
        super().__init__(filename, unit_system=unit_system)
        self._load_tree()

    def _load_tree(self):
        my_tree = {}
        for i, clump_id in enumerate(self.data[("clump", "clump_id")]):
            my_tree[clump_id] = YTClumpContainer(
                clump_id,
                i,
                self.data["clump", "parent_id"][i],
                self.data["clump", "contour_key"][i],
                self.data["clump", "contour_id"][i],
                self,
            )
        for clump in my_tree.values():
            if clump.parent_id == -1:
                self.tree = clump
            else:
                parent = my_tree[clump.parent_id]
                parent.add_child(clump)

    @cached_property
    def leaves(self):
        return [clump for clump in self.tree if clump.children is None]

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if not filename.endswith(".h5"):
            return False
        with h5py.File(filename, mode="r") as f:
            data_type = parse_h5_attr(f, "data_type")
            if data_type is None:
                return False
            if data_type == "yt_clump_tree":
                return True
        return False
