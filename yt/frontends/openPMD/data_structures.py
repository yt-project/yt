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
from yt.utilities.logger import ytLogger as mylog

import sys

class openPMDBasePathException(Exception) :
    pass

class openPMDBasePath :
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
            mylog.warning("openOMD: File does not have valid iteration encoding:")
            mylog.warning(handle.attrs["iterationEncoding"])

        # TODO in the future (see above) this can be a mylog.warning instead of an error
        if len(list_iterations) == 0 :
            raise openPMDBasePathException("openPMD: no iterations found!")

        # just handle the first iteration found
        mylog.warning("openPMD: only choose to load the first iteration ({})".format(handle[dataPath].keys()[0]))
        self.basePath = "{}/{}/".format(dataPath, handle[dataPath].keys()[0])




class openPMDGrid(AMRGridPatch):
    """
    This class defines the characteristics of the grids
    Actually there is only one grid for the whole simolation box
    """
    _id_offset = 0
    __slots__ = ["_level_id"]

    def __init__(self, id, index, level=-1):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        # There is only one grid. So there are no parent or child grids
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

            Parameters
            ----------
            self:
                A reference to self
        """
        # TODO This only parses one file
        f = self.dataset._handle
        meshesPath = f.attrs["meshesPath"]
        particlesPath = f.attrs["particlesPath"]
        output_fields = []

        for group in f[self.basePath + meshesPath].keys():
            try:
                for direction in f[self.basePath + meshesPath + group].keys():
                    output_fields.append(group + "_" + direction)
            except:
                # This is for dataSets, they do not have keys
                output_fields.append(group.replace("_","-"))
        self.field_list = [("openPMD", str(c)) for c in output_fields]

        def is_const_component(record_component):
            return ("value" in record_component.attrs.keys())

        # TODO Pay attention to scalars
        particle_fields = []
        # WHY would this equal false?
        if self.basePath + particlesPath in f:
            for particleName in f[self.basePath + particlesPath].keys():
                for record in f[self.basePath + particlesPath + particleName].keys():
                    if is_const_component(f[self.basePath + particlesPath + particleName + "/" + record]):
                        # Record itself (eg particle_mass) is constant
                        mylog.info("%s", record)
                        particle_fields.append(particleName + "_" + record)
                    elif 'particlePatches' not in record:
                        try:
                            keys = f[self.basePath + particlesPath + particleName + "/" + record].keys()
                            for axis in keys:
                                mylog.info("%s_%s", record, axis)
                                particle_fields.append(particleName + "_" + record + "_" + axis)
                                pass
                        except:
                            mylog.info("%s", record)
                            particle_fields.append(particleName + "_" + record)
                            pass
                    else:
                        # We probably do not want particlePatches as accessible field lists
                        pass
            self.field_list.extend([("io", c) for c in particle_fields])
            #self.field_list.extend([(str(c).split("_")[0], 'particle' + c.lstrip(c.split("_")[0])) for c in particle_fields])

    def _count_grids(self):
        """
            Counts the number of grids in the dataSet. (Only one in the current standard)

            From yt doc:
            this must set self.num_grids to be the total number of grids (equiv AMRGridPatch'es) in the simulation

            Parameters
            ----------
            self:
                A reference to self
        """
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
            Each of these variables is an array, with an entry for each of the self.num_grids grids.
            Additionally, grids must be an array of AMRGridPatch objects that already know their IDs.

            Parameters
            ----------
            self:
                A reference to self
        """
        # This needs to fill the following arrays, where N is self.num_grids:
        meshesPath = self.dataset._handle.attrs["meshesPath"]
        particlesPath = self.dataset._handle.attrs["particlesPath"]

        dims = self.dataset.domain_left_edge.shape[0]#self.dimensionality

        # These objects are expecting 3 values so we pad with zeros in 1D adims 2D
        self.grid_left_edge[
            0,:dims] = self.dataset.domain_left_edge.copy()  # (N, 3) <= float64
        self.grid_right_edge[
            0,:dims] = self.dataset.domain_right_edge.copy()  # (N, 3) <= float64
        self.grid_dimensions[
            0][:dims] = self.dataset.domain_dimensions[:dims] # (N, 3) <= int

# TODO this disables particle reads for now
#      Should might be set in _read_particles for each species,
#      also each species might need its own grid (?)
#        self.grid_particle_count[0] = self.dataset._handle[
#            self.basePath + particlesPath + "/electrons/position/x"].shape[0]  # (N, 1) <= int
        # self.grid_levels = 1           #(N, 1) <= int
        # self.grids = np.empty(1, dtype='object') #(N, 1) <= grid objects

        self.grid_levels.flat[:] = 0
        self.grids = np.empty(self.num_grids, dtype='object')
        # You have to inalize the grids
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i, 0])

    def _populate_grid_objects(self):
        """
            This function initializes the grids

            From yt doc:
            this initializes the grids by calling _prepare_grid() and _setup_dx() on all of them.
            Additionally, it should set up Children and Parent lists on each grid object.

            Parameters
            ----------
            self:
                A reference to self
        """

        # self._reconstruct_parent_child()

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
                 units_override=None):
        # Opens a HDF5 file and stores its file handle in _handle
        # All _handle objects refers to the file
        self._handle = HDF5FileHandler(filename)
        self._filepath = os.path.dirname(filename)
        self._setBasePath(self._handle, self._filepath)
        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename
        self.fluid_types += ('openPMD',)
        #self.particle_types = self._handle[self.basePath + self._handle.attrs["particlesPath"]].keys()
        #self.particle_types = tuple(self.particle_types)
        #self.particle_types += ('all',)
        self.particle_types = ["io", "all"]
        self.particle_types = tuple(self.particle_types)
        self.particle_types_raw = self.particle_types


    def _set_code_unit_attributes(self):
        """
            From yt doc:
            handle conversion between the different physical units and the code units

            Parameters
            ----------
            self:
                A reference to self
        """
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        self.length_unit = self.quan(1.0, "m")
        self.mass_unit = self.quan(1.0, "kg")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "m/s")
        self.magnetic_unit = self.quan(1.0, "gauss")
        #self.magnetic_unit = self.quan(1.0, "T")

    def _parse_parameter_file(self):
        """
            From yt doc:
            read in metadata describing the overall data on disk

            Parameters
            ----------
            self:
                A reference to self
        """
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be updated to be in code units at a later time.  This includes
        # the cosmological parameters.

        f = self._handle
        meshesPath = f.attrs["meshesPath"].decode()
        particlesPath = f.attrs["particlesPath"].decode()
        positionPath = self.basePath + particlesPath + "/electrons/position/"

        # This defines the size of the simulaion box
        self.unique_identifier = 0  # no identifier
        self.parameters = 0  # no additional parameters  <= full of code-specific items of use

        # TODO At this point one assumes the whole file/simulation
        #      contains for all mesh records the same dimensionality and shapes
        # TODO This probably won't work for const records
        # TODO Support particle-only files
        # pick first field
        try :
            firstIteration = list(f["/data/"].keys())[0]
            meshes = f["/data/" + str(firstIteration) + "/" + meshesPath]
            firstMeshName = list(meshes.keys())[0]
            firstMesh = meshes[firstMeshName]
            if type(firstMesh) == h5py.Dataset :
                fshape = firstMesh.shape
            else :
                fshape = firstMesh[list(firstMesh.keys())[0]].shape
        except :
          print("ERROR: Can only read files that have at least one mesh entry!")

        # Usually 2D/3D for picongpu
        self.dimensionality = len(fshape)

        # TODO fill me with actual start and end positions in reasonable units
        self.domain_left_edge = np.zeros(3, dtype=np.float64)
        self.domain_right_edge = np.ones(3, dtype=np.float64)

        # gridding of the meshes (assumed all mesh entries are on the same mesh)
        self.domain_dimensions = np.ones(3, dtype=np.int64)
        self.domain_dimensions[:len(fshape)] = fshape

        # TODO assumes non-peridic boundary conditions
        self.periodicity = np.zeros(3, dtype=np.bool)

        self.current_time = f[self.basePath].attrs[
            "time"]  # <= simulation time in code units

        # Used for AMR, not applicable for us
        # TODO Probably 1
        self.refine_by = 1

        # Not a cosmological simulation
        self.cosmological_simulation = 0  # <= int, 0 or 1
        # Not necessary to set these according to Axel
        # self.current_redshift = 0  # <= float
        # self.omega_lambda = 0  # <= float
        # self.omega_matter = 0  # <= float
        # self.hubble_constant = 0  # <= float

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """
        This function test if the (with yt.load()) a file could be opened with
        this frontend
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
