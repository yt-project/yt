"""
openPMD data structures


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .fields import OpenPMDFieldInfo

from yt.utilities.file_handler import \
    HDF5FileHandler

import yt.frontends.open_pmd.misc as validator

import numpy as np
import os
import re
from math import ceil, floor
from yt.utilities.logger import ytLogger as mylog


def is_const_component(record_component):
    """Determines whether a group or dataset in the HDF5 file is constant.

    Parameters
    ----------
    record_component : HDF5 group or dataset

    Returns
    -------
    bool
        True if constant, False otherwise

    References
    ----------
    .. https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md,
       section 'Constant Record Components'
    """
    return "value" in record_component.attrs.keys()


class OpenPMDGrid(AMRGridPatch):
    """Represents a disjoined chunk of data on-disk.

    The chunks are sliced off the original field along the x-axis.
    This defines the index and offset for every mesh and particle type.
    It also defines further required characteristics for every grid.
    """
    _id_offset = 0
    __slots__ = ["_level_id"]
    particle_index = []   # Every particle species might have distinct hdf5-indices
    particle_offset = []  # and offsets. Contain tuples (ptype, index) and (ptype, offset)
    # TODO (maybe) consider these for every ftype
    mesh_index = 0
    mesh_offset = 0

    def __init__(self, id, index, level=-1, pi=None, po=None, mi=0, mo=0):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        if pi:
            self.particle_index = pi
        if po:
            self.particle_offset = po
        self.mesh_index = mi
        self.mesh_offset = mo
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "OpenPMDGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class OpenPMDHierarchy(GridIndex):
    """
    Defines which fields and particles are created and read from the hard disk
    Furthermore it defines the characteristics of the grids
    """
    grid = OpenPMDGrid

    def __init__(self, ds, dataset_type='openPMD'):
        self.dataset_type = dataset_type
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        GridIndex.__init__(self, ds, dataset_type)

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
            for field in f[bp + mp].keys():
                try:
                    for ax in f[bp + mp + field].keys():
                        mesh_fields.append(field + "_" + ax)
                except AttributeError:
                    # This is for dataSets, they do not have keys
                    mesh_fields.append(field.replace("_","-"))
        except KeyError:
            # There are no mesh fields
            pass
        self.field_list = [("openPMD", str(c)) for c in mesh_fields]

        particle_fields = []
        try:
            for species in f[bp + pp].keys():
                for record in f[bp + pp + species].keys():
                    if is_const_component(f[bp + pp + species + "/" + record]):
                        # Record itself (e.g. particle_mass) is constant
                        particle_fields.append(species + "_" + record)
                    elif "particlePatches" not in record:
                        try:
                            # Create a field for every axis (x,y,z) of every property (position)
                            # of every species (electrons)
                            axes = f[bp + pp + species + "/" + record].keys()
                            if record in "position":
                                record = "positionCoarse"
                            for ax in axes:
                                particle_fields.append(species + "_" + record + "_" + ax)
                        except AttributeError:
                            # Record is a dataset, does not have axes (e.g. weighting)
                            particle_fields.append(species + "_" + record)
                            pass
                    else:
                        # We probably do not want particlePatches as accessible field lists
                        pass
            if len(f[bp + pp].keys()) > 1:
                # There is more than one particle species, use the specific names as field types
                self.field_list.extend(
                    [(str(c).split("_")[0], ("particle_" + "_".join(str(c).split("_")[1:]))) for c in particle_fields])
            else:
                # Only one particle species, fall back to "io"
                self.field_list.extend(
                    [("io", ("particle_" + "_".join(str(c).split("_")[1:]))) for c in particle_fields])
        except KeyError:
            # There are no particle fields
            pass

    def _count_grids(self):
        """Sets ``self.num_grids`` to be the total number of grids in the simulation.

        The number of grids is determined by their respective memory footprint.
        You may change ``gridsize`` accordingly. Increase if you have few cores or run single-threaded.
        """
        # TODO Calculate the ppg not solely on particle count, also on meshsize
        f = self.dataset._handle
        bp = self.dataset.base_path
        mp = self.dataset.meshes_path
        pp = self.dataset.particles_path

        gridsize = 10 * 10**6  # Byte
        self.numparts = {}

        for species in f[bp + pp].keys():
            axis = f[bp + pp + species + "/position"].keys()[0]
            if is_const_component(f[bp + pp + species + "/position/" + axis]):
                self.numparts[species] = f[bp + pp + species + "/position/" + axis].attrs["shape"]
            else:
                self.numparts[species] = f[bp + pp + species + "/position/" + axis].len()
        # Limit particles per grid by resulting memory footprint
        ppg = int(gridsize/(self.dataset.dimensionality*4))  # 4 Byte per value per dimension (f32)
        # Use an upper bound of equally sized grids, last one might be smaller
        self.num_grids = int(ceil(np.max(self.numparts.values()) * ppg**-1))

    def _parse_index(self):
        """
            Parses dimensions from self._handle into self.

            From yt doc:
            this must fill in
                grid_left_edge,
                grid_right_edge,
                grid_particle_count,
                grid_dimensions and
                grid_levels
            with the appropriate information.
            Each of these variables is an array with an entry for each of the self.num_grids grids.
            Additionally, grids must be an array of AMRGridPatch objects that already know their IDs.
        """
        # There is only one refinement level in openPMD
        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype='object')

        nrp = self.numparts.copy()  # Number of remaining particles from the dataset
        pci = {}  # Index for particle chunk
        for spec in nrp:
            pci[spec] = 0
        remaining = self.dataset.domain_dimensions[0]
        meshindex = 0
        meshedge = self.dataset.domain_left_edge.copy()[0]
        for i in range(self.num_grids):
            self.grid_dimensions[i] = self.dataset.domain_dimensions
            prev = remaining
            remaining -= self.grid_dimensions[i][0] * self.num_grids**-1
            self.grid_dimensions[i][0] = int(round(prev, 0) - round(remaining, 0))
            self.grid_left_edge[i] = self.dataset.domain_left_edge.copy()
            self.grid_left_edge[i][0] = meshedge
            self.grid_right_edge[i] = self.dataset.domain_right_edge.copy()
            self.grid_right_edge[i][0] = self.grid_left_edge[i][0] \
                                         + self.grid_dimensions[i][0]\
                                          * self.dataset.domain_dimensions[0]**-1\
                                          * self.dataset.domain_right_edge[0]
            meshedge = self.grid_right_edge[i][0]
            particlecount = []
            particleindex = []
            for spec in self.numparts:
                particleindex += [(spec, pci[spec])]
                if i is (self.num_grids - 1):
                    # The last grid need not be the same size as the previous ones
                    num = nrp[spec]
                else:
                    num = int(floor(self.numparts[spec] * self.num_grids ** -1))
                particlecount += [(spec, num)]
                nrp[spec] -= num
                self.grid_particle_count[i] += num
            self.grids[i] = self.grid(
                i, self, self.grid_levels[i, 0],
                pi=particleindex,
                po=particlecount,
                mi=meshindex,
                mo=self.grid_dimensions[i][0])
            for spec, val in particlecount:
                pci[spec] += val
            meshindex += self.grid_dimensions[i][0]
            remaining -= self.grid_dimensions[i][0]

    def _populate_grid_objects(self):
        """
            This function initializes the grids

            From yt doc:
            this initializes the grids by calling _prepare_grid() and _setup_dx() on all of them.
            Additionally, it should set up Children and Parent lists on each grid object.
        """
        for i in range(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        self.max_level = 0


class OpenPMDDataset(Dataset):
    """
    A dataset object contains all the information of the simulation and
    is intialized with yt.load()
    
    TODO Ideally, a data set object should only contain a single data set.
         afaik, yt.load() can load multiple data sets and also supports
         multiple iteration-loading if done that way, e.g., from a prefix
         of files.
    """
    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo
    _nonstandard = False

    def __init__(self, filename, dataset_type='openPMD',
                 storage_filename=None,
                 units_override=None,
                 unit_system="mks"):
        self._handle = HDF5FileHandler(filename)
        self._filepath = os.path.dirname(filename)
        if self._nonstandard:
            self._set_nonstandard_paths(self._handle)
        else:
            self._set_paths(self._handle, self._filepath)
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.fluid_types += ('openPMD',)
        parts = tuple(str(c) for c in self._handle[self.base_path + self.particles_path].keys())
        if len(parts) > 1:
            # Only use infile particle names if there is more than one species
            self.particle_types = parts
        mylog.debug("open_pmd - self.particle_types: {}".format(self.particle_types))
        self.particle_types_raw = self.particle_types
        self.particle_types = tuple(self.particle_types)

    def _set_nonstandard_paths(self, handle):
        iteration = handle["/data"].keys()[0]
        self.base_path = "/data/{}/".format(iteration)
        self.meshes_path = "fields/"
        self.particles_path = "particles/"

    def _set_paths(self, handle, filepath):
        list_iterations = []
        if u"groupBased" in handle.attrs["iterationEncoding"]:
            for i in list(handle["/data"].keys()):
                list_iterations.append(i)
            mylog.info("open_pmd: found {} iterations in file".format(len(list_iterations)))
        elif u"fileBased" in handle.attrs["iterationEncoding"]:
            regex = u"^" + handle.attrs["iterationFormat"].replace('%T', '[0-9]+') + u"$"
            if filepath is '':
                mylog.warning("open_pmd: For file based iterations, please use absolute file paths!")
                pass
            for filename in os.listdir(filepath):
                if re.match(regex, filename):
                    list_iterations.append(filename)
            mylog.info("open_pmd: found {} iterations in directory".format(len(list_iterations)))
        else:
            mylog.warning(
                "openOMD: File does not have valid iteration encoding: {}".format(handle.attrs["iterationEncoding"]))

        if len(list_iterations) == 0:
            mylog.warning("openOMD: No iterations found!")
        if u"groupBased" in handle.attrs["iterationEncoding"] and len(list_iterations) > 1:
            mylog.warning("open_pmd: only choose to load one iteration ({})".format(handle["/data"].keys()[0]))

        self.base_path = "{}/{}/".format("/data", handle["/data"].keys()[0])
        self.meshes_path = self._handle["/"].attrs["meshesPath"]
        self.particles_path = self._handle["/"].attrs["particlesPath"]

    def _set_code_unit_attributes(self):
        """
            From yt doc:
            handle conversion between the different physical units and the code units
        """
        # We hardcode these to 1.0 since every dataset can have different code <-> physical scaling
        # We get the actual unit by multiplying with "unitSI" when getting our data from disk
        self.length_unit = self.quan(1.0, "m")
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "m/s")
        self.magnetic_unit = self.quan(1.0, "T")

    def _parse_parameter_file(self):
        """
            From yt doc:
            read in metadata describing the overall data on disk
        """
        f = self._handle
        bp = self.base_path
        mp = self.meshes_path

        self.unique_identifier = 0
        self.parameters = 0

        # We assume all fields to have the same shape
        try:
            mesh = f[bp + mp].keys()[0]
            axis = f[bp + mp + "/" + mesh].keys()[0]
            fshape = np.asarray(f[bp + mp + "/" + mesh + "/" + axis].shape, dtype=np.int64)
        except:
            fshape = np.array([1, 1, 1], dtype=np.int64)
            mylog.warning("Could not detect shape of simulated field! "
                          "Assuming a single cell and thus setting fshape to [1, 1, 1]!")
        self.dimensionality = len(fshape)

        fshape = np.append(fshape, np.ones(3 - self.dimensionality))
        self.domain_dimensions = fshape

        self.domain_left_edge = np.zeros(3, dtype=np.float64)
        try:
            mesh = f[bp + mp].keys()[0]
            if self._nonstandard:
                width = f[bp].attrs['cell_width']
                height = f[bp].attrs['cell_height']
                depth = f[bp].attrs['cell_depth']
                spacing = np.asarray([width, height, depth])
                unitSI = f[bp].attrs['unit_length']
            else:
                spacing = np.asarray(f[bp + mp + "/" + mesh].attrs["gridSpacing"])
                unitSI = f[bp + mp + "/" + mesh].attrs["gridUnitSI"]
            self.domain_right_edge = self.domain_dimensions[:spacing.size] * unitSI * spacing
            self.domain_right_edge = np.append(self.domain_right_edge, np.ones(3 - self.domain_right_edge.size))
        except Exception as e:
            mylog.warning("The domain extent could not be calculated! ({}) Setting the field extent to 1m**3! "
                          "This WILL break particle-overplotting!".format(e))
            self.domain_right_edge = np.ones(3, dtype=np.float64)

        if self._nonstandard:
            self.current_time = 0
        else:
            self.current_time = f[bp].attrs["time"]

        self.periodicity = np.zeros(3, dtype=np.bool)
        self.refine_by = 1
        self.cosmological_simulation = 0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """
            Checks whether the supplied file adheres to the required openPMD standards
            and thus can be read by this frontend
        """
        try:
            f = validator.open_file(args[0])
        except:
            return False
        verbose = False
        extension_pic = False
        # root attributes at "/"
        result_array = np.array([0, 0])
        result_array += validator.check_root_attr(f, verbose, extension_pic)

        # Go through all the iterations, checking both the particles
        # and the meshes
        result_array += validator.check_iterations(f, verbose, extension_pic)

        # this might still be a compatible file not fully respecting openPMD standards
        if result_array[0] != 0:
            try:
                if "/data" in f and f["/data"].keys()[0].isdigit():
                    self._nonstandard = True
                    mylog.info("Reading a file with the open_pmd frontend that does not respect standards? "
                                "Just understand that you're on your own for this!")
                    return True
            except:
                return False
        return True
