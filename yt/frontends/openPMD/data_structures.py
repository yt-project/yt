"""
openPMD data structures


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
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
from .fields import openPMDFieldInfo

from yt.utilities.file_handler import \
    HDF5FileHandler

import yt.frontends.openPMD.misc as validator

import h5py
import numpy as np
import os
import re
from math import ceil
from yt.utilities.logger import ytLogger as mylog
from .misc import get_component, is_const_component


class openPMDBasePathException(Exception):
    pass


class openPMDBasePath:
    def  _setBasePath(self, handle, filepath):
        """
        Set the base path for the first iteration found in the file.
        TODO implement into distinct methods:
            - __init__(self, handle)
            - getIterations(self)
            - getBasePath(self, iteration)
        """
        # basePath is fixed in openPMD 1.X to `/data/%T/`
        dataPath = u"/data"

        # if the file messed up the base path we avoid throwing a cluttered
        # exception below while looking for iterations:
        if handle.attrs["basePath"].decode("utf-8") != u"/data/%T/":
            raise openPMDBasePathException("openPMD: basePath is non-standard!")

        # does `/data/` exist?
        if not u"/data" in handle:
            raise openPMDBasePathException("openPMD: group for basePath does not exist!")

        # TODO Everything prior to this should (in theory) already be done by the validator
        # find iterations in basePath
        list_iterations = []
        if u"groupBased" in handle.attrs["iterationEncoding"]:
            for i in list(handle[dataPath].keys()):
                list_iterations.append(i)
            mylog.info("openPMD: found {} iterations in file".format(len(list_iterations)))
        elif u"fileBased" in handle.attrs["iterationEncoding"]:
            regex = u"^" + handle.attrs["iterationFormat"].replace('%T', '[0-9]+') + u"$"
            if filepath is '':
                mylog.warning("openPMD: For file based iterations, please use absolute file paths!")
                pass
            for filename in os.listdir(filepath):
                if re.match(regex, filename):
                    list_iterations.append(filename)
            mylog.info("openPMD: found {} iterations in directory".format(len(list_iterations)))
        else:
            mylog.warning(
                "openOMD: File does not have valid iteration encoding: {}".format(handle.attrs["iterationEncoding"]))

        if len(list_iterations) == 0:
            mylog.warning("openOMD: No iterations found!")

        # just handle the first iteration found
        if u"groupBased" in handle.attrs["iterationEncoding"] and len(list_iterations) > 1:
            mylog.warning("openPMD: only choose to load one iteration ({})".format(handle[dataPath].keys()[0]))
        self.basePath = "{}/{}/".format(dataPath, handle[dataPath].keys()[0])


class openPMDGrid(AMRGridPatch):
    # TODO THIS IS NOT TRUE FOR THE YT SENSE OF "GRID"
    """
    This class defines the characteristics of the grids
    Actually there is only one grid for the whole simulation box
    """
    _id_offset = 0
    __slots__ = ["_level_id"]
    part_ind = 0
    off_part = 0
    mesh_ind = 0
    off_mesh = 0

    def __init__(self, id, index, level=-1, pi=0, op=0, mi=0, om=0):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        self.part_ind = pi
        self.off_part = op
        self.mesh_ind = mi
        self.off_mesh = om
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "openPMDGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class openPMDHierarchy(GridIndex, openPMDBasePath):
    """
    Defines which fields and particles are created and read from the hard disk
    Furthermore it defines the characteristics of the grids
    """
    grid = openPMDGrid

    def __init__(self, ds, dataset_type='openPMD'):
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._setBasePath(self.dataset._handle, self.directory)
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        """
            Parses the dataset to define field names for yt.

            NOTE: Each should be a tuple, where the first element is the on-disk
            fluid type or particle type.  Convention suggests that the on-disk
            fluid type is usually the dataset_type and the on-disk particle type
            (for a single population of particles) is "io".
            look for fluid fields

            From yt doc:
            self.field_list must be populated as a list of strings corresponding to "native" fields in the data files.
        """
        # TODO This only parses one file
        f = self.dataset._handle
        bp = self.basePath
        mp = f.attrs["meshesPath"]
        pp = f.attrs["particlesPath"]
        output_fields = []

        for group in f[bp + mp].keys():
            try:
                for direction in f[bp + mp + group].keys():
                    output_fields.append(group + "_" + direction)
            except:
                # This is for dataSets, they do not have keys
                output_fields.append(group.replace("_","-"))
        self.field_list = [("openPMD", str(c)) for c in output_fields]

        particle_fields = []
        if bp + pp in f:
            for particleName in f[bp + pp].keys():
                for record in f[bp + pp + particleName].keys():
                    if is_const_component(f[bp + pp + particleName + "/" + record]):
                        # Record itself (e.g. particle_mass) is constant
                        particle_fields.append(particleName + "_" + record)
                    elif 'particlePatches' not in record:
                        try:
                            # Create a field for every axis (x,y,z) of every property (position)
                            # of every species (electrons)
                            keys = f[bp + pp + particleName + "/" + record].keys()
                            for axis in keys:
                                particle_fields.append(particleName + "_" + record + "_" + axis)
                        except:
                            # Record is a dataset, does not have axes (e.g. weighting)
                            particle_fields.append(particleName + "_" + record)
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

    def _count_grids(self):
        """
            Counts the number of grids in the dataSet. (Only one in the current standard)

            From yt doc:
            this must set self.num_grids to be the total number of grids (equiv AMRGridPatch'es) in the simulation
        """
        # TODO For the moment we only create grids if there are particles present
        # TODO every species might (read: does) have a different number of active particles
        try:
            f = self.dataset._handle
            bp = self.basePath
            pp = f.attrs["particlesPath"]
            spec = f[bp + pp].keys()[0]
            # We could do this by setting an upper bound for the size of each grid
            # Applying a (rough) 100MB restriction to particle data
            pos = f[bp + pp + spec + "/position"].keys()[0]
            if is_const_component(f[bp + pp + spec + "/position/" + pos]):
                self.np = f[bp + pp + spec + "/position/" + pos].attrs["shape"]
            else:
                self.np = f[bp + pp + spec + "/position/" + pos].len()
            gridsize = 100 * 10**6  # Bytes
            # For 3D: about 8 Mio. particles per grid
            # For 2D: about 12,5 Mio. particles per grid
            # For 1D: about 25 Mio. particles per grid
            self.ppg = int(gridsize/(self.dataset.dimensionality*4))  # 4 Byte per value per dimension (f32)
            # Use an upper bound of equally sized grids, last one might be smaller
            self.num_grids = int(ceil(self.np*self.ppg**-1))
        except:
            mylog.info("Could not detect particlePatch, falling back to single grid!")
            self.num_grids = 1

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

        remaining = self.dataset.domain_dimensions[0]
        for i in range(self.num_grids):
            self.grid_left_edge[i] = self.dataset.domain_left_edge.copy()  # (N, 3) <= float64
            self.grid_left_edge[i][0] = i * self.num_grids**-1
            self.grid_right_edge[i] = self.dataset.domain_right_edge.copy()  # (N, 3) <= float64
            self.grid_right_edge[i][0] = (i + 1) * self.num_grids**-1
            self.grid_dimensions[i] = self.dataset.domain_dimensions  # (N, 3) <= int
            if i != self.num_grids-1:
                self.grid_dimensions[i][0] = int(ceil(self.grid_dimensions[i][0] * self.num_grids**-1))
                remaining -= self.grid_dimensions[i][0]
                nap = self.ppg
            else:
                # The last grid need not be the same size as the previous ones
                self.grid_dimensions[i][0] = remaining
                nap = self.np % self.ppg
            self.grid_particle_count[i] = nap
            self.grids[i] = self.grid(
                i, self, self.grid_levels[i, 0],
                pi=i*self.ppg, op=nap,
                mi=int(ceil(self.dataset.domain_dimensions[0] * self.num_grids**-1)) * i, om=self.grid_dimensions[i][0])

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


class openPMDDataset(Dataset, openPMDBasePath):
    """
    A dataset object contains all the information of the simulation and
    is intialized with yt.load()
    
    TODO Ideally, a data set object should only contain a single data set.
         afaik, yt.load() can load multiple data sets and also supports
         multiple iteration-loading if done that way, e.g., from a prefix
         of files.
    """
    _index_class = openPMDHierarchy
    _field_info_class = openPMDFieldInfo

    def __init__(self, filename, dataset_type='openPMD',
                 storage_filename=None,
                 units_override=None,
                 unit_system="mks"):
        self._handle = HDF5FileHandler(filename)
        self._filepath = os.path.dirname(filename)
        self._setBasePath(self._handle, self._filepath)
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.fluid_types += ('openPMD',)
        parts = tuple(str(c) for c in self._handle[self.basePath + self._handle.attrs["particlesPath"]].keys())
        if len(parts) > 1:
            # Only use infile particle names if there is more than one species
            self.particle_types = parts
        mylog.debug("openPMD - self.particle_types: {}".format(self.particle_types))
        self.particle_types_raw = self.particle_types
        self.particle_types = tuple(self.particle_types)
        self.particle_types = tuple(self.particle_types)


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
        bp = self.basePath
        mp = f.attrs["meshesPath"]
        pp = f.attrs["particlesPath"]

        self.unique_identifier = 0
        self.parameters = 0  # no additional parameters  <= full of code-specific items of use

        # We set the highest dimensionality in the simulation as the global one
        fshape = None
        try:
            for mesh in f[bp + mp].keys():
                for axis in f[bp + mp + "/" + mesh].keys():
                    fshape = np.maximum(fshape, f[bp + mp + "/" + mesh + "/" + axis].shape)
        except:
            pass
        if len(fshape) < 1:
            self.dimensionality = 0
            for species in f[bp + pp].keys():
                self.dimensionality = max(
                    len(f[bp + pp + "/" + species].keys()),
                    self.dimensionality)
        else:
            self.dimensionality = len(fshape)

        # TODO fill me with actual start and end positions in reasonable units
        # This COULD be done with minimum/maximum particle positions in the simulation
        # or through extent of the field meshes
        self.domain_left_edge = np.zeros(3, dtype=np.float64)
        self.domain_right_edge = np.ones(3, dtype=np.float64)

        # gridding of the meshes (assumed all mesh entries are on the same mesh)
        self.domain_dimensions = np.ones(3, dtype=np.int64)
        self.domain_dimensions[:self.dimensionality] = fshape

        self.current_time = f[bp].attrs["time"]
        # We here assume non-periodic boundary conditions
        self.periodicity = np.zeros(3, dtype=np.bool)
        # Used for AMR, not applicable for us
        self.refine_by = 1
        # Not a cosmological simulation
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
        if result_array[0] != 0:
            return False
        return True
