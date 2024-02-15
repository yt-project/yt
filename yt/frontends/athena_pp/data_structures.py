import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.index_subobjects.stretched_grid import StretchedGrid
from yt.data_objects.static_output import Dataset
from yt.fields.magnetic_field import get_magnetic_normalization
from yt.funcs import mylog
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.chemical_formulas import compute_mu
from yt.utilities.file_handler import HDF5FileHandler

from .fields import AthenaPPFieldInfo

geom_map = {
    "cartesian": "cartesian",
    "cylindrical": "cylindrical",
    "spherical_polar": "spherical",
    "minkowski": "cartesian",
    "tilted": "cartesian",
    "sinusoidal": "cartesian",
    "schwarzschild": "spherical",
    "kerr-schild": "spherical",
}


class AthenaPPGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        LE, RE = self.index.grid_left_edge[id, :], self.index.grid_right_edge[id, :]
        self.dds = self.ds.arr((RE - LE) / self.ActiveDimensions, "code_length")
        if self.ds.dimensionality < 2:
            self.dds[1] = 1.0
        if self.ds.dimensionality < 3:
            self.dds[2] = 1.0
        self.field_data["dx"], self.field_data["dy"], self.field_data["dz"] = self.dds


class AthenaPPStretchedGrid(StretchedGrid):
    _id_offset = 0

    def __init__(self, id, cell_widths, index, level):
        super().__init__(id, cell_widths, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level


class AthenaPPHierarchy(GridIndex):
    _dataset_type = "athena_pp"
    _data_file = None

    def __init__(self, ds, dataset_type="athena_pp"):
        self.dataset = weakref.proxy(ds)
        self.grid = AthenaPPStretchedGrid if self.dataset._nonuniform else AthenaPPGrid
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.filename
        self._handle = ds._handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [("athena_pp", k) for k in self.dataset._field_map]

    def _count_grids(self):
        self.num_grids = self._handle.attrs["NumMeshBlocks"]

    def _parse_index(self):
        num_grids = self._handle.attrs["NumMeshBlocks"]

        self.grid_left_edge = np.zeros((num_grids, 3), dtype="float64")
        self.grid_right_edge = np.zeros((num_grids, 3), dtype="float64")
        self.grid_dimensions = np.zeros((num_grids, 3), dtype="int32")

        # TODO: In an unlikely case this would use too much memory, implement
        #       chunked read along 1 dim
        x = self._handle["x1f"][:, :].astype("float64")
        y = self._handle["x2f"][:, :].astype("float64")
        z = self._handle["x3f"][:, :].astype("float64")
        dx = np.diff(x, axis=1)
        dy = np.diff(y, axis=1)
        dz = np.diff(z, axis=1)
        mesh_block_size = self._handle.attrs["MeshBlockSize"]

        for i in range(num_grids):
            self.grid_left_edge[i] = np.array([x[i, 0], y[i, 0], z[i, 0]])
            self.grid_right_edge[i] = np.array([x[i, -1], y[i, -1], z[i, -1]])
            self.grid_dimensions[i] = mesh_block_size
        levels = self._handle["Levels"][:]

        self.grid_left_edge = self.ds.arr(self.grid_left_edge, "code_length")
        self.grid_right_edge = self.ds.arr(self.grid_right_edge, "code_length")

        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(num_grids):
            if self.dataset._nonuniform:
                self.grids[i] = self.grid(i, [dx[i], dy[i], dz[i]], self, levels[i])
            else:
                self.grids[i] = self.grid(i, self, levels[i])

        if self.dataset.dimensionality <= 2:
            self.grid_right_edge[:, 2] = self.dataset.domain_right_edge[2]
        if self.dataset.dimensionality == 1:
            self.grid_right_edge[:, 1:] = self.dataset.domain_right_edge[1:]
        self.grid_particle_count = np.zeros([self.num_grids, 1], dtype="int64")

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self._handle.attrs["MaxLevel"]


class AthenaPPDataset(Dataset):
    _field_info_class = AthenaPPFieldInfo
    _dataset_type = "athena_pp"
    _index_class = AthenaPPHierarchy

    def __init__(
        self,
        filename,
        dataset_type="athena_pp",
        storage_filename=None,
        parameters=None,
        units_override=None,
        unit_system="code",
        default_species_fields=None,
        magnetic_normalization="gaussian",
    ):
        self.fluid_types += ("athena_pp",)
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters
        if units_override is None:
            units_override = {}
        self._handle = HDF5FileHandler(filename)
        xrat = self._handle.attrs["RootGridX1"][2]
        yrat = self._handle.attrs["RootGridX2"][2]
        zrat = self._handle.attrs["RootGridX3"][2]
        self._nonuniform = xrat != 1.0 or yrat != 1.0 or zrat != 1.0
        self._magnetic_factor = get_magnetic_normalization(magnetic_normalization)

        geom = self._handle.attrs["Coordinates"].decode("utf-8")
        self.geometry = Geometry(geom_map[geom])
        if self.geometry == "cylindrical":
            axis_order = ("r", "theta", "z")
        else:
            axis_order = None

        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
            axis_order=axis_order,
        )
        if storage_filename is None:
            storage_filename = self.basename + ".yt"
        self.storage_filename = storage_filename

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file
        """
        if "length_unit" not in self.units_override:
            self.no_cgs_equiv_length = True
        for unit, cgs in [
            ("length", "cm"),
            ("time", "s"),
            ("mass", "g"),
            ("temperature", "K"),
        ]:
            # We set these to cgs for now, but they may have been overridden
            if getattr(self, unit + "_unit", None) is not None:
                continue
            mylog.warning("Assuming 1.0 = 1.0 %s", cgs)
            setattr(self, f"{unit}_unit", self.quan(1.0, cgs))

        self.magnetic_unit = np.sqrt(
            self._magnetic_factor
            * self.mass_unit
            / (self.time_unit**2 * self.length_unit)
        )
        self.magnetic_unit.convert_to_units("gauss")
        self.velocity_unit = self.length_unit / self.time_unit

    def _parse_parameter_file(self):
        xmin, xmax = self._handle.attrs["RootGridX1"][:2]
        ymin, ymax = self._handle.attrs["RootGridX2"][:2]
        zmin, zmax = self._handle.attrs["RootGridX3"][:2]

        self.domain_left_edge = np.array([xmin, ymin, zmin], dtype="float64")
        self.domain_right_edge = np.array([xmax, ymax, zmax], dtype="float64")

        self.domain_width = self.domain_right_edge - self.domain_left_edge
        self.domain_dimensions = self._handle.attrs["RootGridSize"]

        self._field_map = {}
        k = 0
        for dname, num_var in zip(
            self._handle.attrs["DatasetNames"], self._handle.attrs["NumVariables"]
        ):
            for j in range(num_var):
                fname = self._handle.attrs["VariableNames"][k].decode("ascii", "ignore")
                self._field_map[fname] = (dname.decode("ascii", "ignore"), j)
                k += 1

        self.refine_by = 2
        dimensionality = 3
        if self.domain_dimensions[2] == 1:
            dimensionality = 2
        if self.domain_dimensions[1] == 1:
            dimensionality = 1
        self.dimensionality = dimensionality
        self.current_time = self._handle.attrs["Time"]
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.field_ordering = "fortran"
        self.boundary_conditions = [1] * 6
        self._periodicity = tuple(
            self.specified_parameters.get("periodicity", (True, True, True))
        )
        if "gamma" in self.specified_parameters:
            self.gamma = float(self.specified_parameters["gamma"])
        else:
            self.gamma = 5.0 / 3.0

        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0
        self.cosmological_simulation = 0
        self.parameters["Time"] = self.current_time  # Hardcode time conversion for now.
        self.parameters[
            "HydroMethod"
        ] = 0  # Hardcode for now until field staggering is supported.
        if "gamma" in self.specified_parameters:
            self.parameters["Gamma"] = self.specified_parameters["gamma"]
        else:
            self.parameters["Gamma"] = 5.0 / 3.0
        self.mu = self.specified_parameters.get(
            "mu", compute_mu(self.default_species_fields)
        )

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        return filename.endswith(".athdf")

    @property
    def _skip_cache(self):
        return True

    def __str__(self):
        return self.basename.rsplit(".", 1)[0]
