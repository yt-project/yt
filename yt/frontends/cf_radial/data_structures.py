"""
CF Radial data structures



"""
import contextlib
import os
import weakref
from typing import Optional, Tuple

import numpy as np
from unyt import unyt_array

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.file_handler import NetCDF4FileHandler, warn_netcdf
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
        self.grid_levels[0][0] = 0
        self.max_level = 0

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
        storage_filename=None,
        storage_overwrite: bool = False,
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
        storage_overwrite: bool
            if True and if any gridding parameters are set, then the
            storage_filename will be over-written if it exists. Default is False.
        grid_shape : Optional[Tuple[int, int, int]]
            when gridding to cartesian, grid_shape is the number of cells in the
            z, y, x coordinatess. If not provided, yt attempts to calculate a
            reasonable shape based on the resolution of the original cfradial grid
        grid_limit_x : Optional[Tuple[float, float]]
            The x range of the cartesian-gridded data in the form (xmin, xmax) with
            x in the native radar range units
        grid_limit_y : Optional[Tuple[float, float]]
            The y range of the cartesian-gridded data in the form (ymin, ymax) with
            y in the native radar range units
        grid_limit_z : Optional[Tuple[float, float]]
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

                    grid_limit_z = self._validate_grid_dim(radar, "z", grid_limit_z)
                    grid_limit_x = self._validate_grid_dim(radar, "x", grid_limit_x)
                    grid_limit_y = self._validate_grid_dim(radar, "y", grid_limit_y)
                    grid_limits = (grid_limit_z, grid_limit_y, grid_limit_x)
                    grid_shape = self._validate_grid_shape(grid_shape)
                    # note: grid_shape must be in (z, y, x) order.

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

    def _validate_grid_dim(
        self, radar, dim: str, grid_limit: Optional[Tuple[float, float]] = None
    ) -> Tuple[float, float]:
        if grid_limit is None:
            if dim.lower() == "z":
                gate_alt = radar.gate_altitude["data"]
                gate_alt_units = radar.gate_altitude["units"]
                grid_limit = (gate_alt.min(), gate_alt.max())
                grid_limit = self._round_grid_guess(grid_limit, gate_alt_units)
                mylog.info(
                    "grid_limit_z not provided, using max height range in data: (%f, %f)",
                    *grid_limit,
                )
            else:
                max_range = radar.range["data"].max()
                grid_limit = self._round_grid_guess(
                    (-max_range, max_range), radar.range["units"]
                )
                mylog.info(
                    "grid_limit_%s not provided, using max horizontal range in data: (%f, %f)",
                    dim,
                    *grid_limit,
                )
        if len(grid_limit) != 2:
            raise ValueError(
                f"grid_limit_{dim} must have 2 dimensions, but it has {len(grid_limit)}"
            )

        return grid_limit

    def _validate_grid_shape(
        self, grid_shape: Optional[Tuple[int, int, int]] = None
    ) -> Tuple[int, int, int]:
        if grid_shape is None:
            grid_shape = (100, 100, 100)
            mylog.info(
                "grid_shape not provided, using (nz, ny, nx) = (%i, %i, %i)",
                *grid_shape,
            )
        if len(grid_shape) != 3:
            raise ValueError(
                f"grid_shape must have 3 dimensions, but it has {len(grid_shape)}"
            )
        return grid_shape

    def _round_grid_guess(self, bounds: Tuple[float, float], unit_str: str):
        # rounds the bounds to the closest 10 km increment that still contains
        # the grid_limit
        for findstr, repstr in self._field_info_class.unit_subs:
            unit_str = unit_str.replace(findstr, repstr)
        limits = unyt_array(bounds, unit_str).to("km")
        limits[0] = np.floor(limits[0] / 10.0) * 10.0
        limits[1] = np.ceil(limits[1] / 10.0) * 10.0
        return tuple(limits.to(unit_str).tolist())

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

        if not xr.__is_available__:
            return False

        warn_netcdf(filename)
        is_cfrad = False
        try:
            # note that we use the NetCDF4FileHandler here to avoid some
            # issues with xarray opening datasets it cannot handle. Once
            # a dataset is as identified as a CFRadialDataset, xarray is used
            # for opening. See https://github.com/yt-project/yt/issues/3987
            nc4_file = NetCDF4FileHandler(filename)
            with nc4_file.open_ds(keepweakref=True) as ds:
                con = "Conventions"  # the attribute to check for file conventions
                cons = ""  # the value of the Conventions attribute
                for c in [con, con.lower()]:
                    if hasattr(ds, c):
                        cons += getattr(ds, c)
                is_cfrad = "CF/Radial" in cons
        except (OSError, AttributeError, ImportError):
            return False

        return is_cfrad
