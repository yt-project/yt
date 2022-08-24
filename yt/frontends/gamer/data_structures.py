import os
import weakref

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.cosmology import Cosmology
from yt.utilities.file_handler import HDF5FileHandler

from .definitions import geometry_parameters
from .fields import GAMERFieldInfo


class GAMERGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=index.index_filename, index=index)
        self.Parent = None  # do NOT initialize Parent as []
        self.Children = []
        self.Level = level


class GAMERHierarchy(GridIndex):
    grid = GAMERGrid
    _preload_implemented = True  # since gamer defines "_read_chunk_data" in io.py

    def __init__(self, ds, dataset_type="gamer"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._handle = ds._handle
        self._group_grid = ds._group_grid
        self._group_particle = ds._group_particle
        self.float_type = "float64"  # fixed even when FLOAT8 is off
        self._particle_handle = ds._particle_handle
        self.refine_by = ds.refine_by
        self.pgroup = self.refine_by**3  # number of patches in a patch group
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # find all field names in the current dataset
        # grid fields
        self.field_list = [("gamer", v) for v in self._group_grid.keys()]

        # particle fields
        if self._group_particle is not None:
            self.field_list += [("io", v) for v in self._group_particle.keys()]

    def _count_grids(self):
        # count the total number of patches at all levels
        self.num_grids = self.dataset.parameters["NPatch"].sum() // self.pgroup

    def _parse_index(self):
        parameters = self.dataset.parameters
        gid0 = 0
        grid_corner = self._handle["Tree/Corner"][()][:: self.pgroup]
        convert2physical = self._handle["Tree/Corner"].attrs["Cvt2Phy"]

        self.grid_dimensions[:] = parameters["PatchSize"] * self.refine_by

        for lv in range(0, parameters["NLevel"]):
            num_grids_level = parameters["NPatch"][lv] // self.pgroup
            if num_grids_level == 0:
                break

            patch_scale = (
                parameters["PatchSize"] * parameters["CellScale"][lv] * self.refine_by
            )

            # set the level and edge of each grid
            # (left/right_edge are YT arrays in code units)
            self.grid_levels.flat[gid0 : gid0 + num_grids_level] = lv
            self.grid_left_edge[gid0 : gid0 + num_grids_level] = (
                grid_corner[gid0 : gid0 + num_grids_level] * convert2physical
            )
            self.grid_right_edge[gid0 : gid0 + num_grids_level] = (
                grid_corner[gid0 : gid0 + num_grids_level] + patch_scale
            ) * convert2physical

            gid0 += num_grids_level
        self.grid_left_edge += self.dataset.domain_left_edge
        self.grid_right_edge += self.dataset.domain_left_edge

        # allocate all grid objects
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels.flat[i])

        # maximum level with patches (which can be lower than MAX_LEVEL)
        self.max_level = self.grid_levels.max()

        # number of particles in each grid
        try:
            self.grid_particle_count[:] = np.sum(
                self._handle["Tree/NPar"][()].reshape(-1, self.pgroup), axis=1
            )[:, None]
        except KeyError:
            self.grid_particle_count[:] = 0.0

        # calculate the starting particle indices for each grid (starting from 0)
        # --> note that the last element must store the total number of particles
        #    (see _read_particle_coords and _read_particle_fields in io.py)
        self._particle_indices = np.zeros(self.num_grids + 1, dtype="int64")
        np.add.accumulate(
            self.grid_particle_count.squeeze(), out=self._particle_indices[1:]
        )

    def _populate_grid_objects(self):
        son_list = self._handle["Tree/Son"][()]

        for gid in range(self.num_grids):
            grid = self.grids[gid]
            son_gid0 = (
                son_list[gid * self.pgroup : (gid + 1) * self.pgroup] // self.pgroup
            )

            # set up the parent-children relationship
            grid.Children = [self.grids[t] for t in son_gid0[son_gid0 >= 0]]

            for son_grid in grid.Children:
                son_grid.Parent = grid

            # set up other grid attributes
            grid._prepare_grid()
            grid._setup_dx()

        # validate the parent-children relationship in the debug mode
        if self.dataset._debug:
            self._validate_parent_children_relationship()

    # for _debug mode only
    def _validate_parent_children_relationship(self):
        mylog.info("Validating the parent-children relationship ...")

        father_list = self._handle["Tree/Father"][()]

        for grid in self.grids:
            # parent->children == itself
            if grid.Parent is not None:
                assert (
                    grid in grid.Parent.Children
                ), "Grid %d, Parent %d, Parent->Children[0] %d" % (
                    grid.id,
                    grid.Parent.id,
                    grid.Parent.Children[0].id,
                )

            # children->parent == itself
            for c in grid.Children:
                assert c.Parent is grid, "Grid %d, Children %d, Children->Parent %d" % (
                    grid.id,
                    c.id,
                    c.Parent.id,
                )

            # all refinement grids should have parent
            if grid.Level > 0:
                assert (
                    grid.Parent is not None and grid.Parent.id >= 0
                ), "Grid %d, Level %d, Parent %d" % (
                    grid.id,
                    grid.Level,
                    grid.Parent.id if grid.Parent is not None else -999,
                )

            # parent index is consistent with the loaded dataset
            if grid.Level > 0:
                father_gid = father_list[grid.id * self.pgroup] // self.pgroup
                assert (
                    father_gid == grid.Parent.id
                ), "Grid %d, Level %d, Parent_Found %d, Parent_Expect %d" % (
                    grid.id,
                    grid.Level,
                    grid.Parent.id,
                    father_gid,
                )

            # edges between children and parent
            for c in grid.Children:
                for d in range(0, 3):
                    msgL = (
                        "Grid %d, Child %d, Grid->EdgeL %14.7e, Children->EdgeL %14.7e"
                        % (grid.id, c.id, grid.LeftEdge[d], c.LeftEdge[d])
                    )
                    msgR = (
                        "Grid %d, Child %d, Grid->EdgeR %14.7e, Children->EdgeR %14.7e"
                        % (grid.id, c.id, grid.RightEdge[d], c.RightEdge[d])
                    )
                    if not grid.LeftEdge[d] <= c.LeftEdge[d]:
                        raise ValueError(msgL)

                    if not grid.RightEdge[d] >= c.RightEdge[d]:
                        raise ValueError(msgR)

        mylog.info("Check passed")


class GAMERDataset(Dataset):
    _index_class = GAMERHierarchy
    _field_info_class = GAMERFieldInfo
    _handle = None
    _group_grid = None
    _group_particle = None
    _debug = False  # debug mode for the GAMER frontend

    def __init__(
        self,
        filename,
        dataset_type="gamer",
        storage_filename=None,
        particle_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        if self._handle is not None:
            return

        self.fluid_types += ("gamer",)
        self._handle = HDF5FileHandler(filename)
        self.particle_filename = particle_filename

        # to catch both the new and old data formats for the grid data
        try:
            self._group_grid = self._handle["GridData"]
        except KeyError:
            self._group_grid = self._handle["Data"]

        if "Particle" in self._handle:
            self._group_particle = self._handle["Particle"]

        if self.particle_filename is None:
            self._particle_handle = self._handle
        else:
            self._particle_handle = HDF5FileHandler(self.particle_filename)

        # currently GAMER only supports refinement by a factor of 2
        self.refine_by = 2

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
        if self.parameters["Opt__Unit"]:
            # GAMER units are always in CGS
            setdefaultattr(
                self, "length_unit", self.quan(self.parameters["Unit_L"], "cm")
            )
            setdefaultattr(self, "mass_unit", self.quan(self.parameters["Unit_M"], "g"))
            setdefaultattr(self, "time_unit", self.quan(self.parameters["Unit_T"], "s"))

            if self.mhd:
                setdefaultattr(
                    self, "magnetic_unit", self.quan(self.parameters["Unit_B"], "gauss")
                )

        else:
            if len(self.units_override) == 0:
                mylog.warning(
                    "Cannot determine code units ==> "
                    "Use units_override to specify the units"
                )

            for unit, value, cgs in [
                ("length", 1.0, "cm"),
                ("time", 1.0, "s"),
                ("mass", 1.0, "g"),
                ("magnetic", np.sqrt(4.0 * np.pi), "gauss"),
            ]:
                setdefaultattr(self, f"{unit}_unit", self.quan(value, cgs))

                if len(self.units_override) == 0:
                    mylog.warning("Assuming %8s unit = %f %s", unit, value, cgs)

    def _parse_parameter_file(self):

        # code-specific parameters
        for t in self._handle["Info"]:
            info_category = self._handle["Info"][t]
            for v in info_category.dtype.names:
                self.parameters[v] = info_category[v]

        # shortcut for self.parameters
        parameters = self.parameters

        # reset 'Model' to be more readable
        # (no longer regard MHD as a separate model)
        if parameters["Model"] == 1:
            parameters["Model"] = "Hydro"
        elif parameters["Model"] == 3:
            parameters["Model"] = "ELBDM"
        else:
            parameters["Model"] = "Unknown"

        # simulation time and domain
        self.dimensionality = 3  # always 3D
        self.domain_left_edge = parameters.get(
            "BoxEdgeL", np.array([0.0, 0.0, 0.0])
        ).astype("f8")
        self.domain_right_edge = parameters.get(
            "BoxEdgeR", parameters["BoxSize"]
        ).astype("f8")
        self.domain_dimensions = parameters["NX0"].astype("int64")

        # periodicity
        if parameters["FormatVersion"] >= 2106:
            periodic_bc = 1
        else:
            periodic_bc = 0
        self._periodicity = (
            bool(parameters["Opt__BC_Flu"][0] == periodic_bc),
            bool(parameters["Opt__BC_Flu"][2] == periodic_bc),
            bool(parameters["Opt__BC_Flu"][4] == periodic_bc),
        )

        # cosmological parameters
        if parameters["Comoving"]:
            self.cosmological_simulation = 1
            # here parameters["Time"][0] is the scale factor a at certain redshift
            self.current_redshift = 1.0 / parameters["Time"][0] - 1.0
            self.omega_matter = parameters["OmegaM0"]
            self.omega_lambda = 1.0 - self.omega_matter
            # default to 0.7 for old data format
            self.hubble_constant = parameters.get("Hubble0", 0.7)

            # use the cosmological age computed by the given cosmological parameters as the current time when COMOVING is on; cosmological age is computed by subtracting the lookback time at self.current_redshift from that at z = 1e6 (i.e., very early universe)
            cosmo = Cosmology(
                hubble_constant=self.hubble_constant,
                omega_matter=self.omega_matter,
                omega_lambda=self.omega_lambda,
            )
            self.current_time = cosmo.lookback_time(self.current_redshift, 1e6)
        else:
            self.cosmological_simulation = 0
            self.current_redshift = 0.0
            self.omega_matter = 0.0
            self.omega_lambda = 0.0
            self.hubble_constant = 0.0

            # use parameters["Time"][0] as current time when COMOVING is off
            self.current_time = parameters["Time"][0]

        # make aliases to some frequently used variables
        if parameters["Model"] == "Hydro":
            self.gamma = parameters["Gamma"]
            self.eos = parameters.get("EoS", 1)  # Assume gamma-law by default
            # default to 0.6 for old data format
            self.mu = parameters.get(
                "MolecularWeight", 0.6
            )  # Assume ionized primordial by default
            self.mhd = parameters.get("Magnetohydrodynamics", 0)
            self.srhd = parameters.get("SRHydrodynamics", 0)
        else:
            self.mhd = 0
            self.srhd = 0

        # old data format (version < 2210) did not contain any information of code units
        self.parameters.setdefault("Opt__Unit", 0)

        self.geometry = geometry_parameters[parameters.get("Coordinate", 1)]

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        try:
            # define a unique way to identify GAMER datasets
            f = HDF5FileHandler(filename)
            if "Info" in f["/"].keys() and "KeyInfo" in f["/Info"].keys():
                return True
        except Exception:
            pass
        return False
