"""
CF Radial data structures



"""

import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.on_demand_imports import _xarray as xr

from .fields import CFRadialFieldInfo


class CFRadialGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dimensions):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dimensions

    def __repr__(self):
        return "CFRadialGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class CFRadialHierarchy(GridIndex):
    grid = CFRadialGrid

    def __init__(self, ds, dataset_type="cf_radial"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # our index file is the dataset itself:
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        # This sets self.field_list, containing all the available on-disk fields and
        # records the units for each field.
        self.field_list = []
        units = {}
        for key in self.ds._handle.variables.keys():
            if all(x in self.ds._handle[key].dims for x in ["time", "z", "y", "x"]):
                fld = ("cf_radial", key)
                self.field_list.append(fld)
                units[fld] = self.ds._handle[key].units

        self.ds.field_units.update(units)

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 1
        self.max_level = 1

    def _populate_grid_objects(self):
        # only a single grid, no need to loop
        g = self.grid(0, self, self.grid_levels.flat[0], self.grid_dimensions[0])
        g._prepare_grid()
        g._setup_dx()
        self.grids = np.array([g], dtype="object")


class CFRadialDataset(Dataset):
    _index_class = CFRadialHierarchy
    _field_info_class = CFRadialFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="cf_radial",
        grid_shape=(40, 200, 200),
        grid_limits=((0.0, 2000.0), (-1e5, 1e5), (-1e5, 1e5)),
        storage_filename=None,
        units_override=None,
    ):
        self.fluid_types += ("cf_radial",)
        self._handle = xr.open_dataset(filename)

        if "x" not in self._handle.coords:
            if storage_filename is None:
                f_base, f_ext = os.path.splitext(filename)
                storage_filename = f_base + "_grid" + f_ext

            if not os.path.isfile(storage_filename):
                from yt.utilities.on_demand_imports import _pyart as pyart

                radar = pyart.io.read_cfradial(filename)
                self.grid_shape = grid_shape
                self.grid_limits = grid_limits
                grid = pyart.map.grid_from_radars(
                    (radar,), grid_shape=self.grid_shape, grid_limits=self.grid_limits
                )
                mylog.warning(
                    "Saving a cartesian grid for file %s at %s. Data will be loaded "
                    "from the cartesian grid.",
                    filename,
                    storage_filename,
                )
                grid.write(storage_filename)
            self._handle = xr.open_dataset(storage_filename)
            filename = storage_filename
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename
        self.refine_by = 2  # refinement factor between a grid and its subgrid

    def _set_code_unit_attributes(self):
        length_unit = self._handle.variables["x"].attrs["units"]
        self.length_unit = self.quan(1.0, length_unit)
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")

    def _parse_parameter_file(self):
        self.parameters = {}

        x, y, z = (self._handle.coords[d] for d in "xyz")

        self.origin_latitude = self._handle.origin_latitude[0]
        self.origin_longitude = self._handle.origin_longitude[0]

        self.domain_left_edge = np.array([x.min(), y.min(), z.min()])
        self.domain_right_edge = np.array([x.max(), y.max(), z.max()])
        self.dimensionality = 3
        dims = [self._handle.dims[d] for d in "xyz"]
        self.domain_dimensions = np.array(dims, dtype="int64")
        self._periodicity = (False, False, False)
        self.current_time = float(self._handle.time.values)

        # Cosmological information set to zero (not in space).
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.

        try:
            ds = xr.open_dataset(filename)
        except (ImportError, OSError, AttributeError, TypeError):
            # catch all these to avoid errors when xarray cant handle a file
            return False

        if hasattr(ds, "attrs") and isinstance(ds.attrs, dict):
            con = "Conventions"
            return "CF/Radial" in ds.attrs.get(con, "") + ds.attrs.get(con.lower(), "")
        return False
