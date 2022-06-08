"""
CF Radial data structures



"""
import contextlib
import os
import weakref
from typing import Optional, Tuple

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
        with self.ds._handle() as xr_ds_handle:
            for key in xr_ds_handle.variables.keys():
                if all(x in xr_ds_handle[key].dims for x in ["time", "z", "y", "x"]):
                    fld = ("cf_radial", key)
                    self.field_list.append(fld)
                    units[fld] = xr_ds_handle[key].units

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
        storage_filename: Optional[str] = None,
        storage_overwrite: Optional[bool] = False,
        grid_shape: Optional[Tuple[int, int, int]] = None,
        grid_limit_x: Optional[Tuple[float, float]] = None,
        grid_limit_y: Optional[Tuple[float, float]] = None,
        grid_limit_z: Optional[Tuple[float, float]] = None,
        units_override=None,
    ):
        """

        Parameters
        ----------
        filename
        dataset_type
        storage_filename: Optional[str]
            the filename to store gridded file to if necessary. If not provided,
            the string "_yt_grid" will be appended to the dataset filename.
        storage_overwrite: Optional[bool]
            if True and if any gridding parameters are set, then the
            storage_filename will be over-written if it exists. Default is False.
        grid_shape : Optional[Tuple[int, int, int]]
            when gridding to cartesian, grid_shape is the number of cells in the
            z, y, x coordinatess. If not provided, yt attempts to calculate a
            reasonable shape based on the resolution of the original cfradial grid
        grid_limit_x : Optional[Tuple[float, float]]
            The x range of the cartesian-gridded data in the form (xmin, xmax) with
            x in the native radar range units
        grid_limit_y
            The y range of the cartesian-gridded data in the form (ymin, ymax) with
            y in the native radar range units
        grid_limit_z
            The z range of the cartesian-gridded data in the form (zmin, zmax) with
            z in the native radar range units
        units_override
        """
        self.fluid_types += ("cf_radial",)

        with self._handle(filename=filename) as xr_ds_handle:
            if "x" not in xr_ds_handle.coords:
                if storage_filename is None:
                    f_base, f_ext = os.path.splitext(filename)
                    storage_filename = f_base + "_yt_grid" + f_ext

                regrid = True
                if os.path.exists(storage_filename):
                    # pyart grid.write will error if the filename exists, so this logic
                    # forces some explicit behavior to minimize confusion and avoid
                    # overwriting or deleting without explicit user consent.
                    if storage_overwrite:
                        os.remove(storage_filename)
                    elif any([grid_shape, grid_limit_x, grid_limit_y, grid_limit_z]):
                        mylog.warning(
                            "Ignoring provided grid parameters because %s exists.",
                            storage_filename,
                        )
                        mylog.warning(
                            "To re-grid, either provide a unique storage_filename or set "
                            "storage_overwrite to True to overwrite %s.",
                            storage_filename,
                        )
                        regrid = False
                    else:
                        mylog.info(
                            "loading existing re-gridded file: %s", storage_filename
                        )
                        regrid = False

                if regrid:
                    mylog.info("Building cfradial grid")
                    from yt.utilities.on_demand_imports import _pyart as pyart

                    radar = pyart.io.read_cfradial(filename)

                    grid_shape, grid_limits = self._validate_grid(
                        radar, grid_shape, grid_limit_x, grid_limit_y, grid_limit_z
                    )

                    # grid_shape and grid_limits are in (z, y, x) order.
                    self.grid_shape = grid_shape
                    self.grid_limits = grid_limits
                    mylog.info("Calling pyart.map.grid_from_radars ... ")
                    # this is fairly slow
                    grid = pyart.map.grid_from_radars(
                        (radar,),
                        grid_shape=self.grid_shape,
                        grid_limits=self.grid_limits,
                    )
                    mylog.info(
                        "Successfully built cfradial grid, writing to %s",
                        storage_filename,
                    )
                    mylog.info(
                        "Subsequent loads of %s will load the gridded file by default",
                        filename,
                    )
                    grid.write(storage_filename)

                filename = storage_filename

        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename
        self.refine_by = 2  # refinement factor between a grid and its subgrid

    @contextlib.contextmanager
    def _handle(self, filename: Optional[str] = None):
        if filename is None:
            if hasattr(self, "filename"):
                filename = self.filename
            else:
                raise RuntimeError("Dataset has no filename yet.")

        with xr.open_dataset(filename) as xrds:
            yield xrds

    def _validate_grid(
        self,
        radar,
        grid_shape: Optional[Tuple[int, int, int]] = None,
        grid_limit_x: Optional[Tuple[float, float]] = None,
        grid_limit_y: Optional[Tuple[float, float]] = None,
        grid_limit_z: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Tuple, Tuple]:
        # this function sets all the re-gridding parameters if they are not set

        # radar is a pyart radar handle, from pyart.io.read_cfradial()
        max_range_meters = radar.range["data"].max()  # along-ray distance
        elevation = radar.elevation["data"]  # angle above horizontal
        if radar.elevation["units"] == "degrees":
            elevation = elevation * np.pi / 180.0
        max_height = max_range_meters * np.sin(elevation.max())
        max_distance = max_range_meters * np.cos(elevation.max())

        show_help = False
        if grid_limit_z is None:
            # radar.altitude['data'].max() does not seem to be consistently set,
            # so we calculate a max altitude in meters from elevation and range parameters
            grid_limit_z = (0.0, max_height)
            mylog.info(
                "grid_limit_z not provided, using max height range in data: (%f, %f)",
                *grid_limit_z,
            )
            show_help = True

        if grid_limit_x is None:
            grid_limit_x = (-max_distance, max_distance)
            mylog.info(
                "grid_limit_x not provided, using max horizontal range in data: (%f, %f)",
                *grid_limit_x,
            )
            show_help = True

        if grid_limit_y is None:
            grid_limit_y = (-max_distance, max_distance)
            mylog.info(
                "grid_limit_y not provided, using max horizontal range in data: (%f, %f)",
                *grid_limit_y,
            )
            show_help = True

        grid_limits = (grid_limit_z, grid_limit_y, grid_limit_x)

        if grid_shape is None:
            grid_shape = (100, 100, 100)
            mylog.info(
                "grid_shape not provided, using (nz, ny, nx) = (%i, %i, %i)",
                *grid_shape,
            )
            show_help = True

        if show_help:
            mylog.info(
                "to override default grid values, use any of the following "
                "parameters: grid_limit_x, grid_limit_y, grid_limit_z, grid_shape"
                " with yt.load(cf_radial_file, ...)"
            )

        return grid_shape, grid_limits

    def _set_code_unit_attributes(self):
        with self._handle() as xr_ds_handle:
            length_unit = xr_ds_handle.variables["x"].attrs["units"]
        self.length_unit = self.quan(1.0, length_unit)
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")

    def _parse_parameter_file(self):
        self.parameters = {}

        with self._handle() as xr_ds_handle:
            x, y, z = (xr_ds_handle.coords[d] for d in "xyz")

            self.origin_latitude = xr_ds_handle.origin_latitude[0]
            self.origin_longitude = xr_ds_handle.origin_longitude[0]

            self.domain_left_edge = np.array([x.min(), y.min(), z.min()])
            self.domain_right_edge = np.array([x.max(), y.max(), z.max()])
            self.dimensionality = 3
            dims = [xr_ds_handle.dims[d] for d in "xyz"]
            self.domain_dimensions = np.array(dims, dtype="int64")
            self._periodicity = (False, False, False)
            self.current_time = float(xr_ds_handle.time.values)

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
