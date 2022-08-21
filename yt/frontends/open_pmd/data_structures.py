from functools import reduce
from operator import mul
from os import listdir, path
from re import match
from typing import List, Optional

import numpy as np
from packaging.version import Version

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.data_objects.time_series import DatasetSeries
from yt.frontends.open_pmd.fields import OpenPMDFieldInfo
from yt.frontends.open_pmd.misc import get_component, is_const_component
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.file_handler import HDF5FileHandler, warn_h5py
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

ompd_known_versions = [Version(_) for _ in ("1.0.0", "1.0.1", "1.1.0")]
opmd_required_attributes = ["openPMD", "basePath"]


class OpenPMDGrid(AMRGridPatch):
    """Represents chunk of data on-disk.

    This defines the index and offset for every mesh and particle type.
    It also defines parents and children grids. Since openPMD does not have multiple
    levels of refinement there are no parents or children for any grid.
    """

    _id_offset = 0
    __slots__ = ["_level_id"]
    # Every particle species and mesh might have different hdf5-indices and offsets

    ftypes: Optional[List[str]] = []
    ptypes: Optional[List[str]] = []
    findex = 0
    foffset = 0
    pindex = 0
    poffset = 0

    def __init__(self, gid, index, level=-1, fi=0, fo=0, pi=0, po=0, ft=None, pt=None):
        AMRGridPatch.__init__(self, gid, filename=index.index_filename, index=index)
        if ft is None:
            ft = []
        if pt is None:
            pt = []
        self.findex = fi
        self.foffset = fo
        self.pindex = pi
        self.poffset = po
        self.ftypes = ft
        self.ptypes = pt
        self.Parent = None
        self.Children = []
        self.Level = level

    def __str__(self):
        return "OpenPMDGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class OpenPMDHierarchy(GridIndex):
    """Defines which fields and particles are created and read from disk.

    Furthermore it defines the characteristics of the grids.
    """

    grid = OpenPMDGrid

    def __init__(self, ds, dataset_type="openPMD"):
        self.dataset_type = dataset_type
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = path.dirname(self.index_filename)
        GridIndex.__init__(self, ds, dataset_type)

    def _get_particle_type_counts(self):
        """Reads the active number of particles for every species.

        Returns
        -------
        dict
            keys are ptypes
            values are integer counts of the ptype
        """
        result = {}
        f = self.dataset._handle
        bp = self.dataset.base_path
        pp = self.dataset.particles_path

        try:
            for ptype in self.ds.particle_types_raw:
                if str(ptype) == "io":
                    spec = list(f[bp + pp].keys())[0]
                else:
                    spec = ptype
                axis = list(f[bp + pp + "/" + spec + "/position"].keys())[0]
                pos = f[bp + pp + "/" + spec + "/position/" + axis]
                if is_const_component(pos):
                    result[ptype] = pos.attrs["shape"]
                else:
                    result[ptype] = pos.len()
        except (KeyError):
            result["io"] = 0

        return result

    def _detect_output_fields(self):
        """Populates ``self.field_list`` with native fields (mesh and particle) on disk.

        Each entry is a tuple of two strings. The first element is the on-disk fluid
        type or particle type. The second element is the name of the field in yt.
        This string is later used for accessing the data.
        Convention suggests that the on-disk fluid type should be "openPMD",
        the on-disk particle type (for a single species of particles) is "io"
        or (for multiple species of particles) the particle name on-disk.
        """
        f = self.dataset._handle
        bp = self.dataset.base_path
        mp = self.dataset.meshes_path
        pp = self.dataset.particles_path

        mesh_fields = []
        try:
            meshes = f[bp + mp]
            for mname in meshes.keys():
                try:
                    mesh = meshes[mname]
                    for axis in mesh.keys():
                        mesh_fields.append(mname.replace("_", "-") + "_" + axis)
                except AttributeError:
                    # This is a h5py.Dataset (i.e. no axes)
                    mesh_fields.append(mname.replace("_", "-"))
        except (KeyError, TypeError, AttributeError):
            pass
        self.field_list = [("openPMD", str(field)) for field in mesh_fields]

        particle_fields = []
        try:
            particles = f[bp + pp]
            for pname in particles.keys():
                species = particles[pname]
                for recname in species.keys():
                    record = species[recname]
                    if is_const_component(record):
                        # Record itself (e.g. particle_mass) is constant
                        particle_fields.append(
                            pname.replace("_", "-") + "_" + recname.replace("_", "-")
                        )
                    elif "particlePatches" not in recname:
                        try:
                            # Create a field for every axis (x,y,z) of every
                            # property (position) of every species (electrons)
                            axes = list(record.keys())
                            if str(recname) == "position":
                                recname = "positionCoarse"
                            for axis in axes:
                                particle_fields.append(
                                    pname.replace("_", "-")
                                    + "_"
                                    + recname.replace("_", "-")
                                    + "_"
                                    + axis
                                )
                        except AttributeError:
                            # Record is a dataset, does not have axes (e.g. weighting)
                            particle_fields.append(
                                pname.replace("_", "-")
                                + "_"
                                + recname.replace("_", "-")
                            )
                            pass
                    else:
                        pass
            if len(list(particles.keys())) > 1:
                # There is more than one particle species,
                # use the specific names as field types
                self.field_list.extend(
                    [
                        (
                            str(field).split("_")[0],
                            ("particle_" + "_".join(str(field).split("_")[1:])),
                        )
                        for field in particle_fields
                    ]
                )
            else:
                # Only one particle species, fall back to "io"
                self.field_list.extend(
                    [
                        ("io", ("particle_" + "_".join(str(field).split("_")[1:])))
                        for field in particle_fields
                    ]
                )
        except (KeyError, TypeError, AttributeError):
            pass

    def _count_grids(self):
        """Sets ``self.num_grids`` to be the total number of grids in the simulation.

        The number of grids is determined by their respective memory footprint.
        """
        f = self.dataset._handle
        bp = self.dataset.base_path
        mp = self.dataset.meshes_path
        pp = self.dataset.particles_path

        self.meshshapes = {}
        self.numparts = {}

        self.num_grids = 0

        try:
            meshes = f[bp + mp]
            for mname in meshes.keys():
                mesh = meshes[mname]
                if isinstance(mesh, h5py.Group):
                    shape = mesh[list(mesh.keys())[0]].shape
                else:
                    shape = mesh.shape
                spacing = tuple(mesh.attrs["gridSpacing"])
                offset = tuple(mesh.attrs["gridGlobalOffset"])
                unit_si = mesh.attrs["gridUnitSI"]
                self.meshshapes[mname] = (shape, spacing, offset, unit_si)
        except (KeyError, TypeError, AttributeError):
            pass
        try:
            particles = f[bp + pp]
            for pname in particles.keys():
                species = particles[pname]
                if "particlePatches" in species.keys():
                    for (patch, size) in enumerate(
                        species["/particlePatches/numParticles"]
                    ):
                        self.numparts[f"{pname}#{patch}"] = size
                else:
                    axis = list(species["/position"].keys())[0]
                    if is_const_component(species["/position/" + axis]):
                        self.numparts[pname] = species["/position/" + axis].attrs[
                            "shape"
                        ]
                    else:
                        self.numparts[pname] = species["/position/" + axis].len()
        except (KeyError, TypeError, AttributeError):
            pass

        # Limit values per grid by resulting memory footprint
        self.vpg = int(self.dataset.gridsize / 4)  # 4Byte per value (f32)

        # Meshes of the same size do not need separate chunks
        for shape, *_ in set(self.meshshapes.values()):
            self.num_grids += min(
                shape[0], int(np.ceil(reduce(mul, shape) * self.vpg**-1))
            )

        # Same goes for particle chunks if they are not inside particlePatches
        patches = {}
        no_patches = {}
        for (k, v) in self.numparts.items():
            if "#" in k:
                patches[k] = v
            else:
                no_patches[k] = v
        for size in set(no_patches.values()):
            self.num_grids += int(np.ceil(size * self.vpg**-1))
        for size in patches.values():
            self.num_grids += int(np.ceil(size * self.vpg**-1))

    def _parse_index(self):
        """Fills each grid with appropriate properties (extent, dimensions, ...)

        This calculates the properties of every OpenPMDGrid based on the total number of
        grids in the simulation. The domain is divided into ``self.num_grids`` (roughly)
        equally sized chunks along the x-axis. ``grid_levels`` is always equal to 0
        since we only have one level of refinement in openPMD.

        Notes
        -----
        ``self.grid_dimensions`` is rounded to the nearest integer. Grid edges are
        calculated from this dimension. Grids with dimensions [0, 0, 0] are particle
        only. The others do not have any particles affiliated with them.
        """
        f = self.dataset._handle
        bp = self.dataset.base_path
        pp = self.dataset.particles_path

        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype="object")

        grid_index_total = 0

        # Mesh grids
        for mesh in set(self.meshshapes.values()):
            (shape, spacing, offset, unit_si) = mesh
            shape = np.asarray(shape)
            spacing = np.asarray(spacing)
            offset = np.asarray(offset)
            # Total dimension of this grid
            domain_dimension = np.asarray(shape, dtype=np.int32)
            domain_dimension = np.append(
                domain_dimension, np.ones(3 - len(domain_dimension))
            )
            # Number of grids of this shape
            num_grids = min(shape[0], int(np.ceil(reduce(mul, shape) * self.vpg**-1)))
            gle = offset * unit_si  # self.dataset.domain_left_edge
            gre = (
                domain_dimension[: spacing.size] * unit_si * spacing + gle
            )  # self.dataset.domain_right_edge
            gle = np.append(gle, np.zeros(3 - len(gle)))
            gre = np.append(gre, np.ones(3 - len(gre)))
            grid_dim_offset = np.linspace(
                0, domain_dimension[0], num_grids + 1, dtype=np.int32
            )
            grid_edge_offset = (
                grid_dim_offset * float(domain_dimension[0]) ** -1 * (gre[0] - gle[0])
                + gle[0]
            )
            mesh_names = []
            for (mname, mdata) in self.meshshapes.items():
                if mesh == mdata:
                    mesh_names.append(str(mname))
            prev = 0
            for grid in np.arange(num_grids):
                self.grid_dimensions[grid_index_total] = domain_dimension
                self.grid_dimensions[grid_index_total][0] = (
                    grid_dim_offset[grid + 1] - grid_dim_offset[grid]
                )
                self.grid_left_edge[grid_index_total] = gle
                self.grid_left_edge[grid_index_total][0] = grid_edge_offset[grid]
                self.grid_right_edge[grid_index_total] = gre
                self.grid_right_edge[grid_index_total][0] = grid_edge_offset[grid + 1]
                self.grid_particle_count[grid_index_total] = 0
                self.grids[grid_index_total] = self.grid(
                    grid_index_total,
                    self,
                    0,
                    fi=prev,
                    fo=self.grid_dimensions[grid_index_total][0],
                    ft=mesh_names,
                )
                prev += self.grid_dimensions[grid_index_total][0]
                grid_index_total += 1

        handled_ptypes = []

        # Particle grids
        for (species, count) in self.numparts.items():
            if "#" in species:
                # This is a particlePatch
                spec = species.split("#")
                patch = f[bp + pp + "/" + spec[0] + "/particlePatches"]
                domain_dimension = np.ones(3, dtype=np.int32)
                for (ind, axis) in enumerate(list(patch["extent"].keys())):
                    domain_dimension[ind] = patch["extent/" + axis][()][int(spec[1])]
                num_grids = int(np.ceil(count * self.vpg**-1))
                gle = []
                for axis in patch["offset"].keys():
                    gle.append(
                        get_component(patch, "offset/" + axis, int(spec[1]), 1)[0]
                    )
                gle = np.asarray(gle)
                gle = np.append(gle, np.zeros(3 - len(gle)))
                gre = []
                for axis in patch["extent"].keys():
                    gre.append(
                        get_component(patch, "extent/" + axis, int(spec[1]), 1)[0]
                    )
                gre = np.asarray(gre)
                gre = np.append(gre, np.ones(3 - len(gre)))
                np.add(gle, gre, gre)
                npo = patch["numParticlesOffset"][()].item(int(spec[1]))
                particle_count = np.linspace(
                    npo, npo + count, num_grids + 1, dtype=np.int32
                )
                particle_names = [str(spec[0])]
            elif str(species) not in handled_ptypes:
                domain_dimension = self.dataset.domain_dimensions
                num_grids = int(np.ceil(count * self.vpg**-1))
                gle = self.dataset.domain_left_edge
                gre = self.dataset.domain_right_edge
                particle_count = np.linspace(0, count, num_grids + 1, dtype=np.int32)
                particle_names = []
                for (pname, size) in self.numparts.items():
                    if size == count:
                        # Since this is not part of a particlePatch,
                        # we can include multiple same-sized ptypes
                        particle_names.append(str(pname))
                        handled_ptypes.append(str(pname))
            else:
                # A grid with this exact particle count has already been created
                continue
            for grid in np.arange(num_grids):
                self.grid_dimensions[grid_index_total] = domain_dimension
                self.grid_left_edge[grid_index_total] = gle
                self.grid_right_edge[grid_index_total] = gre
                self.grid_particle_count[grid_index_total] = (
                    particle_count[grid + 1] - particle_count[grid]
                ) * len(particle_names)
                self.grids[grid_index_total] = self.grid(
                    grid_index_total,
                    self,
                    0,
                    pi=particle_count[grid],
                    po=particle_count[grid + 1] - particle_count[grid],
                    pt=particle_names,
                )
                grid_index_total += 1

    def _populate_grid_objects(self):
        """This initializes all grids.

        Additionally, it should set up Children and Parent lists on each grid object.
        openPMD is not adaptive and thus there are no Children and Parents for any grid.
        """
        for i in np.arange(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0


class OpenPMDDataset(Dataset):
    """Contains all the required information of a single iteration of the simulation.

    Notes
    -----
    It is assumed that
    - all meshes cover the same region. Their resolution can be different.
    - all particles reside in this same region exclusively.
    - particle and mesh positions are *absolute* with respect to the simulation origin.
    """

    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="openPMD",
        storage_filename=None,
        units_override=None,
        unit_system="mks",
        **kwargs,
    ):
        self._handle = HDF5FileHandler(filename)
        self.gridsize = kwargs.pop("open_pmd_virtual_gridsize", 10**9)
        self.standard_version = Version(self._handle.attrs["openPMD"].decode())
        self.iteration = kwargs.pop("iteration", None)
        self._set_paths(self._handle, path.dirname(filename), self.iteration)
        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.storage_filename = storage_filename
        self.fluid_types += ("openPMD",)
        try:
            particles = tuple(
                str(c)
                for c in self._handle[self.base_path + self.particles_path].keys()
            )
            if len(particles) > 1:
                # Only use on-disk particle names if there is more than one species
                self.particle_types = particles
            mylog.debug("self.particle_types: %s", self.particle_types)
            self.particle_types_raw = self.particle_types
            self.particle_types = tuple(self.particle_types)
        except (KeyError, TypeError, AttributeError):
            pass

    def _set_paths(self, handle, path, iteration):
        """Parses relevant hdf5-paths out of ``handle``.

        Parameters
        ----------
        handle : h5py.File
        path : str
            (absolute) filepath for current hdf5 container
        """
        iterations = []
        if iteration is None:
            iteration = list(handle["/data"].keys())[0]
        encoding = handle.attrs["iterationEncoding"].decode()
        if "groupBased" in encoding:
            iterations = list(handle["/data"].keys())
            mylog.info("Found %s iterations in file", len(iterations))
        elif "fileBased" in encoding:
            itformat = handle.attrs["iterationFormat"].decode().split("/")[-1]
            regex = "^" + itformat.replace("%T", "[0-9]+") + "$"
            if path == "":
                mylog.warning(
                    "For file based iterations, please use absolute file paths!"
                )
                pass
            for filename in listdir(path):
                if match(regex, filename):
                    iterations.append(filename)
            mylog.info("Found %s iterations in directory", len(iterations))

        if len(iterations) == 0:
            mylog.warning("No iterations found!")
        if "groupBased" in encoding and len(iterations) > 1:
            mylog.warning("Only chose to load one iteration (%s)", iteration)

        self.base_path = f"/data/{iteration}/"
        try:
            self.meshes_path = self._handle["/"].attrs["meshesPath"].decode()
            handle[self.base_path + self.meshes_path]
        except (KeyError):
            if self.standard_version <= Version("1.1.0"):
                mylog.info(
                    "meshesPath not present in file. "
                    "Assuming file contains no meshes and has a domain extent of 1m^3!"
                )
                self.meshes_path = None
            else:
                raise
        try:
            self.particles_path = self._handle["/"].attrs["particlesPath"].decode()
            handle[self.base_path + self.particles_path]
        except (KeyError):
            if self.standard_version <= Version("1.1.0"):
                mylog.info(
                    "particlesPath not present in file."
                    " Assuming file contains no particles!"
                )
                self.particles_path = None
            else:
                raise

    def _set_code_unit_attributes(self):
        """Handle conversion between different physical units and the code units.

        Every dataset in openPMD can have different code <-> physical scaling.
        The individual factor is obtained by multiplying with "unitSI" reading getting
        data from disk.
        """
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "m/s"))
        setdefaultattr(self, "magnetic_unit", self.quan(1.0, "T"))

    def _parse_parameter_file(self):
        """Read in metadata describing the overall data on-disk."""
        f = self._handle
        bp = self.base_path
        mp = self.meshes_path

        self.parameters = 0
        self._periodicity = np.zeros(3, dtype="bool")
        self.refine_by = 1
        self.cosmological_simulation = 0

        try:
            shapes = {}
            left_edges = {}
            right_edges = {}
            meshes = f[bp + mp]
            for mname in meshes.keys():
                mesh = meshes[mname]
                if isinstance(mesh, h5py.Group):
                    shape = np.asarray(mesh[list(mesh.keys())[0]].shape)
                else:
                    shape = np.asarray(mesh.shape)
                spacing = np.asarray(mesh.attrs["gridSpacing"])
                offset = np.asarray(mesh.attrs["gridGlobalOffset"])
                unit_si = np.asarray(mesh.attrs["gridUnitSI"])
                le = offset * unit_si
                re = le + shape * unit_si * spacing
                shapes[mname] = shape
                left_edges[mname] = le
                right_edges[mname] = re
            lowest_dim = np.min([len(i) for i in shapes.values()])
            shapes = np.asarray([i[:lowest_dim] for i in shapes.values()])
            left_edges = np.asarray([i[:lowest_dim] for i in left_edges.values()])
            right_edges = np.asarray([i[:lowest_dim] for i in right_edges.values()])
            fs = []
            dle = []
            dre = []
            for i in np.arange(lowest_dim):
                fs.append(np.max(shapes.transpose()[i]))
                dle.append(np.min(left_edges.transpose()[i]))
                dre.append(np.min(right_edges.transpose()[i]))
            self.dimensionality = len(fs)
            self.domain_dimensions = np.append(fs, np.ones(3 - self.dimensionality))
            self.domain_left_edge = np.append(dle, np.zeros(3 - len(dle)))
            self.domain_right_edge = np.append(dre, np.ones(3 - len(dre)))
        except (KeyError, TypeError, AttributeError):
            if self.standard_version <= Version("1.1.0"):
                self.dimensionality = 3
                self.domain_dimensions = np.ones(3, dtype=np.float64)
                self.domain_left_edge = np.zeros(3, dtype=np.float64)
                self.domain_right_edge = np.ones(3, dtype=np.float64)
            else:
                raise

        self.current_time = f[bp].attrs["time"] * f[bp].attrs["timeUnitSI"]

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        """Checks whether the supplied file can be read by this frontend."""
        warn_h5py(filename)
        try:
            with h5py.File(filename, mode="r") as f:
                attrs = list(f["/"].attrs.keys())
                for i in opmd_required_attributes:
                    if i not in attrs:
                        return False

                if Version(f.attrs["openPMD"].decode()) not in ompd_known_versions:
                    return False

                if f.attrs["iterationEncoding"].decode() == "fileBased":
                    return True

                return False
        except (OSError, ImportError):
            return False


class OpenPMDDatasetSeries(DatasetSeries):
    _pre_outputs = ()
    _dataset_cls = OpenPMDDataset
    parallel = True
    setup_function = None
    mixed_dataset_types = False

    def __init__(self, filename):
        super().__init__([])
        self.handle = h5py.File(filename, mode="r")
        self.filename = filename
        self._pre_outputs = sorted(
            np.asarray(list(self.handle["/data"].keys()), dtype="int64")
        )

    def __iter__(self):
        for it in self._pre_outputs:
            ds = self._load(it, **self.kwargs)
            self._setup_function(ds)
            yield ds

    def __getitem__(self, key):
        if isinstance(key, int):
            o = self._load(key)
            self._setup_function(o)
            return o
        else:
            raise KeyError(f"Unknown iteration {key}")

    def _load(self, it, **kwargs):
        return OpenPMDDataset(self.filename, iteration=it)


class OpenPMDGroupBasedDataset(Dataset):
    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo

    def __new__(cls, filename, *args, **kwargs):
        ret = object.__new__(OpenPMDDatasetSeries)
        ret.__init__(filename)
        return ret

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        warn_h5py(filename)
        try:
            with h5py.File(filename, mode="r") as f:
                attrs = list(f["/"].attrs.keys())
                for i in opmd_required_attributes:
                    if i not in attrs:
                        return False

                if Version(f.attrs["openPMD"].decode()) not in ompd_known_versions:
                    return False

                if f.attrs["iterationEncoding"].decode() == "groupBased":
                    return True

                return False
        except (OSError, ImportError):
            return False
