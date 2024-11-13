import os
import warnings
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.fields.magnetic_field import get_magnetic_normalization
from yt.funcs import mylog, setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.chemical_formulas import compute_mu
from yt.utilities.file_handler import HDF5FileHandler

from .fields import ParthenonFieldInfo

_geom_map = {
    "UniformCartesian": Geometry.CARTESIAN,
    "UniformCylindrical": Geometry.CYLINDRICAL,
    "UniformSpherical": Geometry.SPHERICAL,
}

# fmt: off
_cis = np.array(
    [
        [0, 0, 0],
        [0, 0, 1],
        [0, 1, 0],
        [0, 1, 1],
        [1, 0, 0],
        [1, 0, 1],
        [1, 1, 0],
        [1, 1, 1],
    ],
    dtype="int64",
)
# fmt: on


class ParthenonGrid(AMRGridPatch):
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

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False, smoothed=False):
        if smoothed:
            warnings.warn(
                "ghost-zones interpolation/smoothing is not "
                "currently supported for Parthenon data.",
                category=RuntimeWarning,
                stacklevel=2,
            )
            smoothed = False
        return super().retrieve_ghost_zones(
            n_zones, fields, all_levels=all_levels, smoothed=smoothed
        )


class ParthenonHierarchy(GridIndex):
    grid = ParthenonGrid
    _dataset_type = "parthenon"
    _data_file = None

    def __init__(self, ds, dataset_type="parthenon"):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.filename
        self._handle = ds._handle
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [("parthenon", k) for k in self.dataset._field_map]

    def _count_grids(self):
        self.num_grids = self._handle["Info"].attrs["NumMeshBlocks"]

    def _parse_index(self):
        num_grids = self._handle["Info"].attrs["NumMeshBlocks"]

        # TODO: In an unlikely case this would use too much memory, implement
        #       chunked read along 1 dim
        x = self._handle["Locations"]["x"][:, :]
        y = self._handle["Locations"]["y"][:, :]
        z = self._handle["Locations"]["z"][:, :]
        mesh_block_size = self._handle["Info"].attrs["MeshBlockSize"]

        self.grids = np.empty(self.num_grids, dtype="object")
        levels = self._handle["Levels"][:]
        for i in range(num_grids):
            self.grid_left_edge[i] = np.array(
                [x[i, 0], y[i, 0], z[i, 0]], dtype="float64"
            )
            self.grid_right_edge[i] = np.array(
                [x[i, -1], y[i, -1], z[i, -1]], dtype="float64"
            )
            self.grid_dimensions[i] = mesh_block_size
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
        self.max_level = self._handle["Info"].attrs["MaxLevel"]


class ParthenonDataset(Dataset):
    _load_requirements = ["h5py"]
    _field_info_class = ParthenonFieldInfo
    _dataset_type = "parthenon"
    _index_class = ParthenonHierarchy

    def __init__(
        self,
        filename,
        dataset_type="parthenon",
        storage_filename=None,
        parameters=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
        magnetic_normalization="gaussian",
    ):
        self.fluid_types += ("parthenon",)
        if parameters is None:
            parameters = {}
        self.specified_parameters = parameters
        if units_override is None:
            units_override = {}
        self._handle = HDF5FileHandler(filename)
        xrat = self._handle["Info"].attrs["RootGridDomain"][2]
        yrat = self._handle["Info"].attrs["RootGridDomain"][5]
        zrat = self._handle["Info"].attrs["RootGridDomain"][8]
        if xrat != 1.0 or yrat != 1.0 or zrat != 1.0:
            raise NotImplementedError(
                "Logarithmic grids not yet supported/tested in Parthenon frontend."
            )

        self._magnetic_factor = get_magnetic_normalization(magnetic_normalization)

        self.geometry = _geom_map[self._handle["Info"].attrs["Coordinates"]]

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
        for unit, cgs in [
            ("length", "cm"),
            ("time", "s"),
            ("mass", "g"),
        ]:
            unit_param = f"Hydro/code_{unit}_cgs"
            # Use units, if provided in output
            if unit_param in self.parameters:
                setdefaultattr(
                    self, f"{unit}_unit", self.quan(self.parameters[unit_param], cgs)
                )
            # otherwise use code = cgs
            else:
                mylog.warning(f"Assuming 1.0 code_{unit} = 1.0 {cgs}")
                setdefaultattr(self, f"{unit}_unit", self.quan(1.0, cgs))

        self.magnetic_unit = np.sqrt(
            self._magnetic_factor
            * self.mass_unit
            / (self.time_unit**2 * self.length_unit)
        )
        self.magnetic_unit.convert_to_units("gauss")
        self.velocity_unit = self.length_unit / self.time_unit

    def _parse_parameter_file(self):
        self.parameters.update(self.specified_parameters)
        for key, val in self._handle["Params"].attrs.items():
            if key in self.parameters.keys():
                mylog.warning(
                    f"Overriding existing {key!r} key in ds.parameters from data 'Params'"
                )
            self.parameters[key] = val

        xmin, xmax = self._handle["Info"].attrs["RootGridDomain"][0:2]
        ymin, ymax = self._handle["Info"].attrs["RootGridDomain"][3:5]
        zmin, zmax = self._handle["Info"].attrs["RootGridDomain"][6:8]

        self.domain_left_edge = np.array([xmin, ymin, zmin], dtype="float64")
        self.domain_right_edge = np.array([xmax, ymax, zmax], dtype="float64")

        self.domain_width = self.domain_right_edge - self.domain_left_edge
        self.domain_dimensions = self._handle["Info"].attrs["RootGridSize"]

        self._field_map = {}
        k = 0

        dnames = self._handle["Info"].attrs["OutputDatasetNames"]
        num_components = self._handle["Info"].attrs["NumComponents"]

        if "OutputFormatVersion" in self._handle["Info"].attrs.keys():
            self.output_format_version = self._handle["Info"].attrs[
                "OutputFormatVersion"
            ]
        else:
            raise NotImplementedError("Could not determine OutputFormatVersion.")

        # For a single variable, we need to convert it to a list for the following
        # zip to work.
        if isinstance(num_components, np.uint64):
            dnames = (dnames,)
            num_components = (num_components,)

        component_name_offset = 0
        for dname, num_component in zip(dnames, num_components, strict=False):
            for j in range(num_component):
                fname = self._handle["Info"].attrs["ComponentNames"][
                    j + component_name_offset
                ]
                self._field_map[fname] = (dname, j)
                k += 1
            component_name_offset = int(component_name_offset + num_component)

        self.refine_by = 2
        dimensionality = 3
        if self.domain_dimensions[2] == 1:
            dimensionality = 2
        if self.domain_dimensions[1] == 1:
            dimensionality = 1
        self.dimensionality = dimensionality
        self.current_time = self._handle["Info"].attrs["Time"]
        self.num_ghost_zones = 0
        self.field_ordering = "fortran"
        self.boundary_conditions = [1] * 6
        self.cosmological_simulation = False

        if "periodicity" in self.parameters:
            self._periodicity = tuple(self.parameters["periodicity"])
        else:
            boundary_conditions = self._handle["Info"].attrs["BoundaryConditions"]

            inner_bcs = boundary_conditions[::2]
            # outer_bcs = boundary_conditions[1::2]
            ##Check self consistency
            # for inner_bc,outer_bc in zip(inner_bcs,outer_bcs):
            #    if( (inner_bc == "periodicity" or outer_bc == "periodic") and inner_bc != outer_bc ):
            #        raise Exception("Inconsistent periodicity in boundary conditions")

            self._periodicity = tuple(bc == "periodic" for bc in inner_bcs)

        if "gamma" in self.parameters:
            self.gamma = float(self.parameters["gamma"])
        elif "Hydro/AdiabaticIndex" in self.parameters:
            self.gamma = self.parameters["Hydro/AdiabaticIndex"]
        else:
            mylog.warning(
                "Adiabatic index gamma could not be determined. Falling back to 5/3."
            )
            self.gamma = 5.0 / 3.0

        if "mu" in self.parameters:
            self.mu = self.parameters["mu"]
        elif "Hydro/mu" in self.parameters:
            self.mu = self.parameters["Hydro/mu"]
        # Legacy He_mass_fraction parameter implemented in AthenaPK
        elif "Hydro/He_mass_fraction" in self.parameters:
            He_mass_fraction = self.parameters["Hydro/He_mass_fraction"]
            self.mu = 1 / (He_mass_fraction * 3.0 / 4.0 + (1 - He_mass_fraction) * 2)
        # Fallback to primorial gas composition (and show warning)
        else:
            mylog.warning(
                "Plasma composition could not be determined in data file. Falling back to fully ionized primodial composition."
            )
            self.mu = self.parameters.get("mu", compute_mu(self.default_species_fields))

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if cls._missing_load_requirements():
            return False
        return filename.endswith((".phdf", ".rhdf"))

    @property
    def _skip_cache(self):
        return True

    def __str__(self):
        return self.basename.rsplit(".", 1)[0]
