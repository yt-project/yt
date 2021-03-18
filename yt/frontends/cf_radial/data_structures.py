"""
CF Radial data structures



"""

import os
import stat
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
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _initialize_state_variables(self):
        super()._initialize_state_variables()

    def _detect_output_fields(self):
        # This sets self.field_list, containing all the available on-disk fields and
        # records the units for each field.
        self.field_list = []
        units = {}
        for key in self.ds._handle.variables.keys():
            if (
                all(x in self.ds._handle[key].dims for x in ["time", "z", "y", "x"])
                is True
            ):
                fld = ("cf_radial", key)
                self.field_list.append(fld)
                units[fld] = self.ds._handle[key].units

        self.ds.field_units.update(units)

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        # This needs to fill the following arrays, where N is self.num_grids:
        #   self.grid_left_edge         (N, 3) <= float64
        #   self.grid_right_edge        (N, 3) <= float64
        #   self.grid_dimensions        (N, 3) <= int
        #   self.grid_particle_count    (N, 1) <= int
        #   self.grid_levels            (N, 1) <= int
        #   self.grids                  (N, 1) <= grid objects
        #   self.max_level = self.grid_levels.max()
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            g = self.grid(i, self, self.grid_levels.flat[i], self.grid_dimensions[i])
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g
        self.max_level = 1

    def _populate_grid_objects(self):
        # For each grid, this must call:
        #   grid._prepare_grid()
        #   grid._setup_dx()
        # This must also set:
        #   grid.Children <= list of child grids
        #   grid.Parent   <= parent grid
        # This is handled by the frontend because often the children must be
        # identified.
        pass


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
        self.refine_by = 2

        if "x" not in self._handle.coords:
            if storage_filename is None:
                f_base, f_ext = os.path.splitext(filename)
                storage_filename = f_base + "_grid" + f_ext

            if not os.path.isfile(storage_filename):
                from yt.utilities.on_demand_imports import _pyart

                pyart = _pyart.pyart

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
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        length_unit = self._handle.variables["x"].attrs["units"]
        self.length_unit = self.quan(1.0, length_unit)
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.

        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        #   self.parameters             <= full of code-specific items of use
        self.parameters = {}

        x, y, z = [self._handle.coords[d] for d in "xyz"]

        self.origin_latitude = self._handle.origin_latitude[0]
        self.origin_longitude = self._handle.origin_longitude[0]

        #   self.domain_left_edge       <= array of float64
        dle = [x.min(), y.min(), z.min()]
        self.domain_left_edge = np.array(dle)
        #   self.domain_right_edge      <= array of float64
        dre = [x.max(), y.max(), z.max()]
        self.domain_right_edge = np.array(dre)
        #   self.dimensionality         <= int
        self.dimensionality = 3
        #   self.domain_dimensions      <= array of int64
        dims = [self._handle.dims[d] for d in "xyz"]
        self.domain_dimensions = np.array(dims, dtype="int64")
        #   self.periodicity            <= three-element tuple of booleans
        self._periodicity = (False, False, False)
        #   self.current_time           <= simulation time in code units
        self.current_time = float(self._handle.time.values)
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1
        #   self.current_redshift           <= float
        #   self.omega_lambda               <= float
        #   self.omega_matter               <= float
        #   self.hubble_constant            <= float
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
            import xarray
        except ImportError:
            return False

        try:
            ds = xarray.open_dataset(filename)
        except OSError:
            return False

        conventions = ""
        if hasattr(ds, "attrs") and isinstance(ds.attrs, dict):
            conventions += ds.attrs.get("Conventions", "")
            conventions += ds.attrs.get("conventions", "")
        return "CF/Radial" in conventions
