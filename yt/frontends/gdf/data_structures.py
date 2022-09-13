import os
import sys
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import just_one, setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.units.dimensions import dimensionless as sympy_one  # type: ignore
from yt.units.unit_object import Unit  # type: ignore
from yt.units.unit_systems import unit_system_registry  # type: ignore
from yt.utilities.exceptions import YTGDFUnknownGeometry
from yt.utilities.lib.misc_utilities import get_box_grids_level
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import GDFFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property
GEOMETRY_TRANS = {
    0: "cartesian",
    1: "polar",
    2: "cylindrical",
    3: "spherical",
}


class GDFGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, start, dimensions):
        AMRGridPatch.__init__(self, id, filename=index.index_filename, index=index)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.ds.refine_by
        else:
            LE, RE = self.index.grid_left_edge[id, :], self.index.grid_right_edge[id, :]
            self.dds = np.array((RE - LE) / self.ActiveDimensions)
        self.field_data["dx"], self.field_data["dy"], self.field_data["dz"] = self.dds
        self.dds = self.ds.arr(self.dds, "code_length")


class GDFHierarchy(GridIndex):

    grid = GDFGrid

    def __init__(self, ds, dataset_type="grid_data_format"):
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        h5f = h5py.File(self.index_filename, mode="r")
        self.dataset_type = dataset_type
        GridIndex.__init__(self, ds, dataset_type)
        self.directory = os.path.dirname(self.index_filename)
        h5f.close()

    def _detect_output_fields(self):
        h5f = h5py.File(self.index_filename, mode="r")
        self.field_list = [("gdf", str(f)) for f in h5f["field_types"].keys()]
        h5f.close()

    def _count_grids(self):
        h5f = h5py.File(self.index_filename, mode="r")
        self.num_grids = h5f["/grid_parent_id"].shape[0]
        h5f.close()

    def _parse_index(self):
        h5f = h5py.File(self.index_filename, mode="r")
        dxs = []
        self.grids = np.empty(self.num_grids, dtype="object")
        levels = (h5f["grid_level"][:]).copy()
        glis = (h5f["grid_left_index"][:]).copy()
        gdims = (h5f["grid_dimensions"][:]).copy()
        active_dims = ~(
            (np.max(gdims, axis=0) == 1) & (self.dataset.domain_dimensions == 1)
        )

        for i in range(levels.shape[0]):
            self.grids[i] = self.grid(i, self, levels[i], glis[i], gdims[i])
            self.grids[i]._level_id = levels[i]

            dx = (
                self.dataset.domain_right_edge - self.dataset.domain_left_edge
            ) / self.dataset.domain_dimensions
            dx[active_dims] /= self.dataset.refine_by ** levels[i]
            dxs.append(dx.in_units("code_length"))
        dx = self.dataset.arr(dxs, units="code_length")
        self.grid_left_edge = self.dataset.domain_left_edge + dx * glis
        self.grid_dimensions = gdims.astype("int32")
        self.grid_right_edge = self.grid_left_edge + dx * self.grid_dimensions
        self.grid_particle_count = h5f["grid_particle_count"][:]
        del levels, glis, gdims
        h5f.close()

    def _populate_grid_objects(self):
        mask = np.empty(self.grids.size, dtype="int32")
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

        for gi, g in enumerate(self.grids):
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
            get_box_grids_level(
                self.grid_left_edge[gi, :],
                self.grid_right_edge[gi, :],
                self.grid_levels[gi],
                self.grid_left_edge,
                self.grid_right_edge,
                self.grid_levels,
                mask,
            )
            m = mask.astype("bool")
            m[gi] = False
            siblings = self.grids[gi:][m[gi:]]
            if len(siblings) > 0:
                g.OverlappingSiblings = siblings.tolist()
        self.max_level = self.grid_levels.max()

    def _get_box_grids(self, left_edge, right_edge):
        """
        Gets back all the grids between a left edge and right edge
        """
        eps = np.finfo(np.float64).eps
        grid_i = np.where(
            np.all((self.grid_right_edge - left_edge) > eps, axis=1)
            & np.all((right_edge - self.grid_left_edge) > eps, axis=1)
        )

        return self.grids[grid_i], grid_i

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype="bool")
        grids, grid_ind = self._get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]


class GDFDataset(Dataset):
    _index_class = GDFHierarchy
    _field_info_class = GDFFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="grid_data_format",
        storage_filename=None,
        geometry=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.geometry = geometry
        self.fluid_types += ("gdf",)
        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units
        based on the parameter file
        """

        # This should be improved.
        h5f = h5py.File(self.parameter_filename, mode="r")
        for field_name in h5f["/field_types"]:
            current_field = h5f[f"/field_types/{field_name}"]
            if "field_to_cgs" in current_field.attrs:
                field_conv = current_field.attrs["field_to_cgs"]
                self.field_units[field_name] = just_one(field_conv)
            elif "field_units" in current_field.attrs:
                field_units = current_field.attrs["field_units"]
                if isinstance(field_units, str):
                    current_field_units = current_field.attrs["field_units"]
                else:
                    current_field_units = just_one(current_field.attrs["field_units"])
                self.field_units[field_name] = current_field_units.decode("utf8")
            else:
                self.field_units[field_name] = ""

        if "dataset_units" in h5f:
            for unit_name in h5f["/dataset_units"]:
                current_unit = h5f[f"/dataset_units/{unit_name}"]
                value = current_unit[()]
                unit = current_unit.attrs["unit"]
                # need to convert to a Unit object and check dimensions
                # because unit can be things like
                # 'dimensionless/dimensionless**3' so naive string
                # comparisons are insufficient
                unit = Unit(unit, registry=self.unit_registry)
                if unit_name.endswith("_unit") and unit.dimensions is sympy_one:
                    # Catch code units and if they are dimensionless,
                    # assign CGS units. setdefaultattr will catch code units
                    # which have already been set via units_override.
                    un = unit_name[:-5]
                    un = un.replace("magnetic", "magnetic_field_cgs", 1)
                    unit = unit_system_registry["cgs"][un]
                    setdefaultattr(self, unit_name, self.quan(value, unit))
                setdefaultattr(self, unit_name, self.quan(value, unit))
                if unit_name in h5f["/field_types"]:
                    if unit_name in self.field_units:
                        mylog.warning(
                            "'field_units' was overridden by 'dataset_units/%s'",
                            unit_name,
                        )
                    self.field_units[unit_name] = str(unit)
        else:
            setdefaultattr(self, "length_unit", self.quan(1.0, "cm"))
            setdefaultattr(self, "mass_unit", self.quan(1.0, "g"))
            setdefaultattr(self, "time_unit", self.quan(1.0, "s"))

        h5f.close()

    @cached_property
    def unique_identifier(self) -> str:
        with h5py.File(self.parameter_filename, mode="r") as handle:
            return str(handle["/simulation_parameters"].attrs["unique_identifier"])

    def _parse_parameter_file(self):
        self._handle = h5py.File(self.parameter_filename, mode="r")
        if "data_software" in self._handle["gridded_data_format"].attrs:
            self.data_software = self._handle["gridded_data_format"].attrs[
                "data_software"
            ]
        else:
            self.data_software = "unknown"
        sp = self._handle["/simulation_parameters"].attrs
        if self.geometry is None:
            geometry = just_one(sp.get("geometry", 0))
            try:
                self.geometry = GEOMETRY_TRANS[geometry]
            except KeyError as e:
                raise YTGDFUnknownGeometry(geometry) from e
        self.parameters.update(sp)
        self.domain_left_edge = sp["domain_left_edge"][:]
        self.domain_right_edge = sp["domain_right_edge"][:]
        self.domain_dimensions = sp["domain_dimensions"][:]
        refine_by = sp["refine_by"]
        if refine_by is None:
            refine_by = 2
        self.refine_by = refine_by
        self.dimensionality = sp["dimensionality"]
        self.current_time = sp["current_time"]
        self.cosmological_simulation = sp["cosmological_simulation"]
        if sp["num_ghost_zones"] != 0:
            raise RuntimeError
        self.num_ghost_zones = sp["num_ghost_zones"]
        self.field_ordering = sp["field_ordering"]
        self.boundary_conditions = sp["boundary_conditions"][:]
        self._periodicity = tuple(bnd == 0 for bnd in self.boundary_conditions[::2])
        if self.cosmological_simulation:
            self.current_redshift = sp["current_redshift"]
            self.omega_lambda = sp["omega_lambda"]
            self.omega_matter = sp["omega_matter"]
            self.hubble_constant = sp["hubble_constant"]
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0
        self.parameters["Time"] = 1.0  # Hardcode time conversion for now.
        # Hardcode for now until field staggering is supported.
        self.parameters["HydroMethod"] = 0
        self._handle.close()
        del self._handle

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            fileh = h5py.File(filename, mode="r")
            if "gridded_data_format" in fileh:
                fileh.close()
                return True
            fileh.close()
        except Exception:
            pass
        return False

    def __str__(self):
        return self.basename.rsplit(".", 1)[0]
