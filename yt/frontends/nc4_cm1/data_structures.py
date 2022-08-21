import os
import weakref
from collections import OrderedDict

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.file_handler import NetCDF4FileHandler, warn_netcdf
from yt.utilities.logger import ytLogger as mylog

from .fields import CM1FieldInfo


class CM1Grid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dimensions):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dimensions


class CM1Hierarchy(GridIndex):
    grid = CM1Grid

    def __init__(self, ds, dataset_type="cm1"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        # build list of on-disk fields for dataset_type 'cm1'
        vnames = self.dataset.parameters["variable_names"]
        self.field_list = [("cm1", vname) for vname in vnames]

    def _count_grids(self):
        # This needs to set self.num_grids
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 0
        self.max_level = 0

    def _populate_grid_objects(self):
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            g = self.grid(i, self, self.grid_levels.flat[i], self.grid_dimensions[i])
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g


class CM1Dataset(Dataset):
    _index_class = CM1Hierarchy
    _field_info_class = CM1FieldInfo

    def __init__(
        self,
        filename,
        dataset_type="cm1",
        storage_filename=None,
        units_override=None,
        unit_system="mks",
    ):
        self.fluid_types += ("cm1",)
        self._handle = NetCDF4FileHandler(filename)
        # refinement factor between a grid and its subgrid.
        self.refine_by = 1
        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.storage_filename = storage_filename

    def _setup_coordinate_handler(self):
        # ensure correct ordering of axes so plots aren't rotated (z should always be
        # on the vertical axis).
        super()._setup_coordinate_handler()
        self.coordinates._x_pairs = (("x", "y"), ("y", "x"), ("z", "x"))
        self.coordinates._y_pairs = (("x", "z"), ("y", "z"), ("z", "y"))

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        with self._handle.open_ds() as _handle:
            length_unit = _handle.variables["xh"].units
        self.length_unit = self.quan(1.0, length_unit)
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "m/s")
        self.time_unit = self.quan(1.0, "s")

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.

        self.parameters = {}  # code-specific items
        with self._handle.open_ds() as _handle:
            # _handle here is a netcdf Dataset object, we need to parse some metadata
            # for constructing our yt ds.

            # TO DO: generalize this to be coordinate variable name agnostic in order to
            # make useful for WRF or climate data. For now, we're hard coding for CM1
            # specifically and have named the classes appropriately. Additionally, we
            # are only handling the cell-centered grid ("xh","yh","zh") at present.
            # The cell-centered grid contains scalar fields and interpolated velocities.
            dims = [_handle.dimensions[i].size for i in ["xh", "yh", "zh"]]
            xh, yh, zh = (_handle.variables[i][:] for i in ["xh", "yh", "zh"])
            self.domain_left_edge = np.array(
                [xh.min(), yh.min(), zh.min()], dtype="float64"
            )
            self.domain_right_edge = np.array(
                [xh.max(), yh.max(), zh.max()], dtype="float64"
            )

            # loop over the variable names in the netCDF file, record only those on the
            # "zh","yh","xh" grid.
            varnames = []
            for key, var in _handle.variables.items():
                if all(x in var.dimensions for x in ["time", "zh", "yh", "xh"]):
                    varnames.append(key)
            self.parameters["variable_names"] = varnames
            self.parameters["lofs_version"] = _handle.cm1_lofs_version
            self.parameters["is_uniform"] = _handle.uniform_mesh
            self.current_time = _handle.variables["time"][:][0]

            # record the dimension metadata: __handle.dimensions contains netcdf
            # objects so we need to manually copy over attributes.
            dim_info = OrderedDict()
            for dim, meta in _handle.dimensions.items():
                dim_info[dim] = meta.size
            self.parameters["dimensions"] = dim_info

        self.dimensionality = 3
        self.domain_dimensions = np.array(dims, dtype="int64")
        self._periodicity = (False, False, False)

        # Set cosmological information to zero for non-cosmological.
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.

        warn_netcdf(filename)
        try:
            nc4_file = NetCDF4FileHandler(filename)
            with nc4_file.open_ds(keepweakref=True) as _handle:
                is_cm1_lofs = hasattr(_handle, "cm1_lofs_version")
                is_cm1 = hasattr(_handle, "cm1 version")  # not a typo, it is a space...

                # ensure coordinates of each variable array exists in the dataset
                coords = _handle.dimensions  # get the dataset wide coordinates
                failed_vars = []  # list of failed variables
                for var in _handle.variables:  # iterate over the variables
                    vcoords = _handle[var].dimensions  # get the dims for the variable
                    ncoords = len(vcoords)  # number of coordinates in variable
                    # number of coordinates that pass for a variable
                    coordspassed = sum(vc in coords for vc in vcoords)
                    if coordspassed != ncoords:
                        failed_vars.append(var)

                if failed_vars:
                    mylog.warning(
                        "Trying to load a cm1_lofs netcdf file but the "
                        "coordinates of the following fields do not match the "
                        "coordinates of the dataset: %s",
                        failed_vars,
                    )
                    return False

            if not is_cm1_lofs:
                if is_cm1:
                    mylog.warning(
                        "It looks like you are trying to load a cm1 netcdf file, "
                        "but at present yt only supports cm1_lofs output. Until"
                        " support is added, you can likely use"
                        " yt.load_uniform_grid() to load your cm1 file manually."
                    )
                return False
        except (OSError, AttributeError, ImportError):
            return False

        return True
