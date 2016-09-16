"""
openPMD data structures


"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from functools import reduce
from operator import mul
from os import \
    path, \
    listdir
from re import match

import numpy as np

from yt.data_objects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.frontends.open_pmd.fields import OpenPMDFieldInfo
from yt.frontends.open_pmd.misc import is_const_component
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.file_handler import HDF5FileHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5


class OpenPMDGrid(AMRGridPatch):
    """Represents chunk of data on-disk.

    This defines the index and offset for every mesh and particle type.
    It also defines parents and children grids. Since openPMD does not have multiple levels of refinement,
    there are no parents or children for any grid.
    """
    _id_offset = 0
    __slots__ = ["_level_id"]
    # Every particle species and mesh might have different hdf5-indices and offsets
    ftypes=[]
    ptypes=[]
    findex = 0
    foffset = 0
    pindex = 0
    poffset = 0

    def __init__(self, gid, index, level=-1, fi=0, fo=0, pi=0, po=0, ft=[], pt=[]):
        AMRGridPatch.__init__(self, gid, filename=index.index_filename,
                              index=index)
        self.findex = fi
        self.foffset = fo
        self.pindex = pi
        self.poffset = po
        self.ftypes = ft
        self.ptypes = pt
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
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
        """Sets the active number of particles for every species.

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
        for ptype in self.ds.particle_types_raw:
            if ptype in "io":
                spec = list(f[bp + pp].keys())[0]
            else:
                spec = ptype
            axis = list(f[bp + pp + "/" + spec + "/position"].keys())[0]
            pos = f[bp + pp + "/" + spec + "/position/" + axis]
            if is_const_component(pos):
                result[ptype] = pos.attrs["shape"]
            else:
                result[ptype] = pos.len()
        return result

    def _detect_output_fields(self):
        """Populates ``self.field_list`` with native fields (mesh and particle) on disk.

        Each entry is a tuple of two strings. The first element is the on-disk fluid type or particle type.
        The second element is the name of the field in yt. This string is later used for accessing the data.
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
            for field in list(f[bp + mp].keys()):
                try:
                    for axis in list(f[bp + mp + field].keys()):
                        mesh_fields.append(field + "_" + axis)
                except AttributeError:
                    # This is a h5.Dataset (i.e. no axes)
                    mesh_fields.append(field.replace("_", "-"))
        except KeyError:
            # There are no mesh fields
            pass
        self.field_list = [("openPMD", str(field)) for field in mesh_fields]

        particle_fields = []
        try:
            for species in list(f[bp + pp].keys()):
                for record in list(f[bp + pp + species].keys()):
                    if is_const_component(f[bp + pp + species + "/" + record]):
                        # Record itself (e.g. particle_mass) is constant
                        particle_fields.append(species + "_" + record)
                    elif "particlePatches" not in record:
                        try:
                            # Create a field for every axis (x,y,z) of every property (position)
                            # of every species (electrons)
                            axes = list(f[bp + pp + species + "/" + record].keys())
                            if record in "position":
                                record = "positionCoarse"
                            for axis in axes:
                                particle_fields.append(species + "_" + record + "_" + axis)
                        except AttributeError:
                            # Record is a dataset, does not have axes (e.g. weighting)
                            particle_fields.append(species + "_" + record)
                            pass
                    else:
                        pass
            if len(list(f[bp + pp].keys())) > 1:
                # There is more than one particle species, use the specific names as field types
                self.field_list.extend(
                    [(str(field).split("_")[0],
                      ("particle_" + "_".join(str(field).split("_")[1:]))) for field in particle_fields])
            else:
                # Only one particle species, fall back to "io"
                self.field_list.extend(
                    [("io",
                      ("particle_" + "_".join(str(field).split("_")[1:]))) for field in particle_fields])
        except KeyError:
            # There are no particle fields
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

        for mesh in list(f[bp + mp].keys()):
            if type(f[bp + mp + mesh]) is h5.Group:
                self.meshshapes[mesh] = f[bp + mp + mesh + "/" + list(f[bp + mp + mesh].keys())[0]].shape
            else:
                self.meshshapes[mesh] = f[bp + mp + mesh].shape
        for species in list(f[bp + pp].keys()):
            if "particlePatches" in list(f[bp + pp + "/" + species].keys()):
                for patch, size in enumerate(f[bp + pp + "/" + species + "/particlePatches/numParticles"]):
                    self.numparts[species + "#" + str(patch)] = size
            else:
                axis = list(f[bp + pp + species + "/position"].keys())[0]
                if is_const_component(f[bp + pp + species + "/position/" + axis]):
                    self.numparts[species] = f[bp + pp + species + "/position/" + axis].attrs["shape"]
                else:
                    self.numparts[species] = f[bp + pp + species + "/position/" + axis].len()

        # Limit values per grid by resulting memory footprint
        self.vpg = int(self.dataset.gridsize / 4)  # 4Byte per value (f32)

        # Meshes of the same size do not need separate chunks
        for shape in set(self.meshshapes.values()):
            self.num_grids += min(shape[0], int(np.ceil(reduce(mul, shape) * self.vpg**-1)))

        # Same goes for particle chunks if they are not inside particlePatches
        patches = {}
        no_patches = {}
        for k, v in list(self.numparts.items()):
            if "#" in k:
                patches[k] = v
            else:
                no_patches[k] = v
        for size in set(no_patches.values()):
            self.num_grids += int(np.ceil(size * self.vpg ** -1))
        for size in list(patches.values()):
            self.num_grids += int(np.ceil(size * self.vpg ** -1))

    def _parse_index(self):
        """Fills each grid with appropriate properties (extent, dimensions, ...)

        This calculates the properties of every OpenPMDGrid based on the total number of grids in the simulation.
        The domain is divided into ``self.num_grids`` (roughly) equally sized chunks along the x-axis.
        ``grid_levels`` is always equal to 0 since we only have one level of refinement in openPMD.

        Notes
        -----
        ``self.grid_dimensions`` is rounded to the nearest integer. Grid edges are calculated from this dimension.
        Grids with dimensions [0, 0, 0] are particle only. The others do not have any particles affiliated with them.
        """
        f = self.dataset._handle
        bp = self.dataset.base_path
        pp = self.dataset.particles_path

        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype="object")

        def get_component(group, component_name, index, offset):
            record_component = group[component_name]
            unit_si = record_component.attrs["unitSI"]
            return np.multiply(record_component[index:index+offset], unit_si)

        grid_index_total = 0

        # Mesh grids
        for shape in set(self.meshshapes.values()):
            # Total dimension of this grid
            domain_dimension = np.asarray(shape, dtype=np.int32)
            domain_dimension = np.append(domain_dimension, np.ones(3 - len(domain_dimension)))
            # Number of grids of this shape
            num_grids = min(shape[0], int(np.ceil(reduce(mul, shape) * self.vpg ** -1)))
            gle = self.dataset.domain_left_edge
            gre = self.dataset.domain_right_edge
            grid_dim_offset = np.linspace(0, domain_dimension[0], num_grids + 1, dtype=np.int32)
            grid_edge_offset = grid_dim_offset * np.float(domain_dimension[0]) ** -1 * (gre[0] - gle[0]) + gle[0]
            mesh_names = []
            for (mname, mshape) in list(self.meshshapes.items()):
                if shape == mshape:
                    mesh_names.append(str(mname))
            prev = 0
            for grid in range(num_grids):
                self.grid_dimensions[grid_index_total] = domain_dimension
                self.grid_dimensions[grid_index_total][0] = grid_dim_offset[grid + 1] - grid_dim_offset[grid]
                self.grid_left_edge[grid_index_total] = gle
                self.grid_left_edge[grid_index_total][0] = grid_edge_offset[grid]
                self.grid_right_edge[grid_index_total] = gre
                self.grid_right_edge[grid_index_total][0] = grid_edge_offset[grid + 1]
                self.grid_particle_count[grid_index_total] = 0
                self.grids[grid_index_total] = self.grid(grid_index_total, self, 0,
                                                         fi=prev,
                                                         fo=self.grid_dimensions[grid_index_total][0],
                                                         ft=mesh_names)
                prev += self.grid_dimensions[grid_index_total][0]
                grid_index_total += 1

        handled_ptypes = []

        # Particle grids
        for (species, count) in list(self.numparts.items()):
            if "#" in species:
                # This is a particlePatch
                spec = species.split("#")
                patch = f[bp + pp + "/" + spec[0] + "/particlePatches"]
                num_grids = int(np.ceil(count * self.vpg ** -1))
                gle = []
                for axis in list(patch["offset"].keys()):
                    gle.append(get_component(patch, "offset/" + axis, int(spec[1]), 1)[0])
                gle = np.asarray(gle)
                gle = np.append(gle, np.zeros(3 - len(gle)))
                gre = []
                for axis in list(patch["extent"].keys()):
                    gre.append(get_component(patch, "extent/" + axis, int(spec[1]), 1)[0])
                gre = np.asarray(gre)
                gre = np.append(gre, np.ones(3 - len(gre)))
                np.add(gle, gre, gre)
                npo = patch["numParticlesOffset"].value.item(int(spec[1]))
                particle_count = np.linspace(npo, npo + count, num_grids + 1,
                                             dtype=np.int32)
                particle_names = [str(spec[0])]
            elif str(species) not in handled_ptypes:
                num_grids = int(np.ceil(count * self.vpg ** -1))
                gle = self.dataset.domain_left_edge
                gre = self.dataset.domain_right_edge
                particle_count = np.linspace(0, count, num_grids + 1, dtype=np.int32)
                particle_names = []
                for (pname, size) in list(self.numparts.items()):
                    if size == count:
                        # Since this is not part of a particlePatch, we can include multiple same-sized ptypes
                        particle_names.append(str(pname))
                        handled_ptypes.append(str(pname))
            else:
                # A grid with this exact particle count has already been created
                continue
            for grid in range(num_grids):
                self.grid_dimensions[grid_index_total] = [0, 0, 0]  # Counted as mesh-size, thus no dimensional extent
                self.grid_left_edge[grid_index_total] = gle
                self.grid_right_edge[grid_index_total] = gre
                self.grid_particle_count[grid_index_total] = (particle_count[grid + 1] - particle_count[grid]) * len(
                    particle_names)
                self.grids[grid_index_total] = self.grid(grid_index_total, self, 0,
                                                         pi=particle_count[grid],
                                                         po=particle_count[grid + 1] - particle_count[grid],
                                                         pt=particle_names)
                grid_index_total += 1

    def _populate_grid_objects(self):
        """This initializes all grids.

        Additionally, it should set up Children and Parent lists on each grid object.
        openPMD is not adaptive and thus there are no Children and Parents for any grid.
        """
        for i in range(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0


class OpenPMDDataset(Dataset):
    """Contains all the required information of a single iteration of the simulation.

    Notes
    -----
    It is assumed that all meshes cover the same region. Their resolution can be different.
    It is assumed that all particles reside in this same region exclusively.
    It is assumed that the particle and mesh positions are *absolute* with respect to the simulation origin.
    """
    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo

    def __init__(self,
                 filename,
                 dataset_type="openPMD",
                 storage_filename=None,
                 units_override=None,
                 unit_system="mks",
                 **kwargs):
        self._handle = HDF5FileHandler(filename)
        self._set_paths(self._handle, path.dirname(filename))
        Dataset.__init__(self,
                         filename,
                         dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.fluid_types += ("openPMD",)
        particles = tuple(str(c) for c in list(self._handle[self.base_path + self.particles_path].keys()))
        if len(particles) > 1:
            # Only use on-disk particle names if there is more than one species
            self.particle_types = particles
        mylog.debug("open_pmd - self.particle_types: {}".format(self.particle_types))
        self.gridsize = kwargs.pop("open_pmd_virtual_gridsize", 10**9)
        self.particle_types_raw = self.particle_types
        self.particle_types = tuple(self.particle_types)

    def _set_paths(self, handle, path):
        """Parses relevant hdf5-paths out of ``handle``.

        Parameters
        ----------
        handle : h5py.File
        path : str
            (absolute) filepath for current hdf5 container
        """
        iterations = []
        encoding = handle.attrs["iterationEncoding"].decode()
        if "groupBased" in encoding:
            iterations = list(handle["/data"].keys())
            mylog.info("open_pmd - found {} iterations in file".format(len(iterations)))
        elif "fileBased" in encoding:
            regex = "^" + handle.attrs["iterationFormat"].decode().replace("%T", "[0-9]+") + "$"
            if path is "":
                mylog.warning("open_pmd - For file based iterations, please use absolute file paths!")
                pass
            for filename in listdir(path):
                if match(regex, filename):
                    iterations.append(filename)
            mylog.info("open_pmd - found {} iterations in directory".format(len(iterations)))

        if len(iterations) == 0:
            mylog.warning("open_pmd - no iterations found!")
        if "groupBased" in encoding and len(iterations) > 1:
            mylog.warning("open_pmd - only choose to load one iteration ({})".format(list(handle["/data"].keys())[0]))

        self.base_path = "/data/{}/".format(list(handle["/data"].keys())[0])
        self.meshes_path = self._handle["/"].attrs["meshesPath"].decode()
        self.particles_path = self._handle["/"].attrs["particlesPath"].decode()

    def _set_code_unit_attributes(self):
        """Handle conversion between different physical units and the code units.

        Every dataset in openPMD can have different code <-> physical scaling.
        The individual factor is obtained by multiplying with "unitSI" reading getting data from disk.
        """
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "m/s"))
        setdefaultattr(self, "magnetic_unit", self.quan(1.0, "T"))

    def _parse_parameter_file(self):
        """Read in metadata describing the overall data on-disk.
        """
        f = self._handle
        bp = self.base_path
        mp = self.meshes_path

        self.unique_identifier = 0
        self.parameters = 0
        self.periodicity = np.zeros(3, dtype=np.bool)
        self.refine_by = 1
        self.cosmological_simulation = 0

        try:
            mesh = list(f[bp + mp].keys())[0]
            axis = list(f[bp + mp + "/" + mesh].keys())[0]
            fshape = np.asarray(f[bp + mp + "/" + mesh + "/" + axis].shape, dtype=np.int64)
        except:
            fshape = np.array([1, 1, 1], dtype=np.int64)
            mylog.warning("open_pmd - Could not detect shape of simulated field! "
                          "Assuming a single cell and thus setting fshape to [1, 1, 1]!")
        self.dimensionality = len(fshape)

        fshape = np.append(fshape, np.ones(3 - self.dimensionality))
        self.domain_dimensions = fshape

        try:
            mesh = list(f[bp + mp].keys())[0]
            spacing = np.asarray(f[bp + mp + "/" + mesh].attrs["gridSpacing"])
            offset = np.asarray(f[bp + mp + "/" + mesh].attrs["gridGlobalOffset"])
            unit_si = np.asarray(f[bp + mp + "/" + mesh].attrs["gridUnitSI"])
            self.domain_left_edge = offset * unit_si
            # self.domain_left_edge = np.zeros(3, dtype=np.float64)
            self.domain_right_edge = self.domain_dimensions[:spacing.size] * unit_si * spacing
            self.domain_right_edge += self.domain_left_edge
            self.domain_left_edge = np.append(self.domain_left_edge, np.zeros(3 - self.domain_left_edge.size))
            self.domain_right_edge = np.append(self.domain_right_edge, np.ones(3 - self.domain_right_edge.size))
        except:
            mylog.warning(
                "open_pmd - The domain extent could not be calculated! Setting the field extent to 1m**3!")
            self.domain_left_edge = np.zeros(3, dtype=np.float64)
            self.domain_right_edge = np.ones(3, dtype=np.float64)

        self.current_time = f[bp].attrs["time"] * f[bp].attrs["timeUnitSI"]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """Checks whether the supplied file can be read by this frontend.
        """
        try:
            f = h5.File(args[0])
        except:
            return False

        requirements = ["openPMD", "basePath", "meshesPath", "particlesPath"]
        attrs = list(f["/"].attrs.keys())
        for i in requirements:
            if i not in attrs:
                return False

        versions = ["1.0.0", "1.0.1"]
        for i in versions:
            if i in f.attrs["openPMD"].decode():
                return True

        return False
