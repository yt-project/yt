"""
CF Radial data structures



"""

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
            x, y, z coordinates. If not provided, yt attempts to calculate a
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
        self._handle = xr.open_dataset(filename)

        if "x" not in self._handle.coords:
            if storage_filename is None:
                f_base, f_ext = os.path.splitext(filename)
                storage_filename = f_base + "_yt_grid" + f_ext

            regrid = True
            if os.path.exists(storage_filename):
                # pyart grid.write will error if the filename exists, so this logic
                # forces some explicit behavior to minimize confusion and avoid
                # overwriting or deleting without explicit user input.
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
                    mylog.info("loading existing re-gridded file: %s", storage_filename)
                    regrid = False

            if regrid:
                from yt.utilities.on_demand_imports import _pyart as pyart

                radar = pyart.io.read_cfradial(filename)
                grid_shape, grid_limits = self._validate_grid(
                    radar, grid_shape, grid_limit_x, grid_limit_y, grid_limit_z
                )

                self.grid_shape = grid_shape
                self.grid_limits = grid_limits

                mylog.info("Building cfradial grid and writing to %s", storage_filename)
                grid = pyart.map.grid_from_radars(
                    (radar,), grid_shape=self.grid_shape, grid_limits=self.grid_limits
                )
                grid.write(storage_filename)

            self._handle = xr.open_dataset(storage_filename)
            filename = storage_filename

        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename
        self.refine_by = 2  # refinement factor between a grid and its subgrid

    def _validate_grid(
        self,
        radar,
        grid_shape: Optional[Tuple[int, int, int]] = None,
        grid_limit_x: Optional[Tuple[float, float]] = None,
        grid_limit_y: Optional[Tuple[float, float]] = None,
        grid_limit_z: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Tuple, Tuple]:

        # radar: a pyart radar handle, from pyart.io.read_cfradial
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
            mylog.warning(
                "grid_limit_z not provided, using max height range in data: %f",
                grid_limit_z,
            )
            show_help = True

        if grid_limit_x is None:
            grid_limit_x = (-max_distance, max_distance)
            mylog.warning(
                "grid_limit_x not provided, using max horizontal range in data: %f",
                grid_limit_x,
            )
            show_help = True

        if grid_limit_y is None:
            grid_limit_y = (-max_distance, max_distance)
            mylog.warning(
                "grid_limit_y not provided, using max horizontal range in data: %f",
                grid_limit_y,
            )
            show_help = True

        grid_limits = (grid_limit_x, grid_limit_y, grid_limit_z)

        if grid_shape is None:
            # choose a resolution that respects the original grid to some degree
            # the following uses only the radar.elevation and radar.range coordinates.
            # It ignores radar.azimuth and so may result in some empty data when
            # using pyart.map.grid_from_radars if the azimuth does not span a full 360
            # degrees. In that case, the user will have to adjust parameters to
            # improve the gridding...

            # characteristic change in along-ray distance
            rdata = radar.range["data"]
            x_vals = rdata * np.cos(elevation.min())
            d_dist = (x_vals[1:] - x_vals[:-1]).min()
            dx = d_dist
            dy = dx

            # representative element close to available resolution in z
            z_vals = max_range_meters * np.sin(elevation)
            dz = np.nanmean(z_vals[1:] - z_vals[:-1])

            # calculate number of elements in x, y, z for that grid resolution
            grid_l = np.array(grid_limits)
            grid_w = grid_l[:, 1] - grid_l[:, 0]  # width in x, y, z
            grid_d = np.array([dx, dy, dz])
            grid_shape = grid_w / grid_d

            # we do not want the initial gridding to take a long time, so
            # limit the max resolution
            max_allowed = 400
            if grid_shape.max() > max_allowed:
                # reduce, but maintain aspect ratio if we can
                aspect = grid_shape / grid_shape.max()
                grid_shape = aspect * max_allowed

            # limit the min resolution
            grid_shape[grid_shape < 100] = 100
            grid_shape = tuple(grid_shape.astype(int).tolist())

            mylog.warning("estimating a good grid_shape: %s", str(grid_shape))
            show_help = True

        if show_help:
            mylog.warning(
                "to override default grid values, use any of the following "
                "parameters: grid_limit_x, grid_limit_y, grid_limit_z, grid_shape"
            )

        return grid_shape, grid_limits

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
