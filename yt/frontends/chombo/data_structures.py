"""
Data structures for Chombo.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import re
import os
import weakref
import numpy as np

from six import string_types
from stat import \
    ST_CTIME

from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.extern import six
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.file_handler import \
    HDF5FileHandler
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.lib.misc_utilities import \
    get_box_grids_level

from .fields import ChomboFieldInfo, Orion2FieldInfo, \
    ChomboPICFieldInfo1D, ChomboPICFieldInfo2D, ChomboPICFieldInfo3D, \
    PlutoFieldInfo


def is_chombo_hdf5(fn):
    try:
        with h5py.File(fn, 'r') as fileh:
            valid = "Chombo_global" in fileh["/"]
    except (KeyError, IOError, ImportError):
        return False
    return valid


class ChomboGrid(AMRGridPatch):
    _id_offset = 0
    __slots__ = ["_level_id", "stop_index"]

    def __init__(self, id, index, level, start, stop):
        AMRGridPatch.__init__(self, id, filename = index.index_filename,
                              index = index)
        self._parent_id = []
        self._children_ids = []
        self.Level = level
        self.ActiveDimensions = stop - start + 1

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        if self.start_index is not None:
            return self.start_index
        if self.Parent is None:
            iLE = self.LeftEdge - self.ds.domain_left_edge
            start_index = iLE / self.dds
            return np.rint(start_index).astype('int64').ravel()
        pdx = self.Parent[0].dds
        start_index = (self.Parent[0].get_global_startindex()) + \
            np.rint((self.LeftEdge - self.Parent[0].LeftEdge)/pdx)
        self.start_index = (start_index*self.ds.refine_by).astype('int64').ravel()
        return self.start_index

    def _setup_dx(self):
        # has already been read in and stored in index
        self.dds = self.ds.arr(self.index.dds_list[self.Level], "code_length")

    @property
    def Parent(self):
        if len(self._parent_id) == 0:
            return None
        return [self.index.grids[pid - self._id_offset]
                for pid in self._parent_id]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset]
                for cid in self._children_ids]


class ChomboHierarchy(GridIndex):

    grid = ChomboGrid
    _data_file = None

    def __init__(self, ds, dataset_type='chombo_hdf5'):
        self.domain_left_edge = ds.domain_left_edge
        self.domain_right_edge = ds.domain_right_edge
        self.dataset_type = dataset_type
        self.field_indexes = {}
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = os.path.abspath(
            self.dataset.parameter_filename)
        self.directory = ds.fullpath
        self._handle = ds._handle

        self._levels = [
            key for key in self._handle.keys() if key.startswith('level')
        ]

        GridIndex.__init__(self, ds, dataset_type)

        self._read_particles()

    def _read_particles(self):

        # only do anything if the dataset contains particles
        if not any([f[1].startswith('particle_') for f in self.field_list]):
            return

        self.num_particles = 0
        particles_per_grid = []
        for key, val in self._handle.items():
            if key.startswith('level'):
                level_particles = val['particles:offsets'][:]
                self.num_particles += level_particles.sum()
                particles_per_grid = np.concatenate((particles_per_grid, level_particles))

        for i, grid in enumerate(self.grids):
            self.grids[i].NumberOfParticles = particles_per_grid[i]
            self.grid_particle_count[i] = particles_per_grid[i]

        assert(self.num_particles == self.grid_particle_count.sum())

    # Chombo datasets, by themselves, have no "known" fields. However, 
    # we will look for "fluid" fields by finding the string "component" in
    # the output file, and "particle" fields by finding the string "particle".
    def _detect_output_fields(self):

        # look for fluid fields
        output_fields = []
        for key, val in self._handle.attrs.items():
            if key.startswith("component"):
                output_fields.append(val.decode("ascii"))
        self.field_list = [("chombo", c) for c in output_fields]

        # look for particle fields
        particle_fields = []
        for key, val in self._handle.attrs.items():
            if key.startswith("particle"):
                particle_fields.append(val.decode("ascii"))
        self.field_list.extend([("io", c) for c in particle_fields])

    def _count_grids(self):
        self.num_grids = 0
        for lev in self._levels:
            d = self._handle[lev]
            if 'Processors' in d:
                self.num_grids += d['Processors'].len()
            elif 'boxes' in d:
                self.num_grids += d['boxes'].len()
            else:
                raise RuntimeError("Uknown file specification")

    def _parse_index(self):
        f = self._handle # shortcut
        self.max_level = f.attrs['num_levels'] - 1

        grids = []
        self.dds_list = []
        i = 0
        D = self.dataset.dimensionality
        for lev_index, lev in enumerate(self._levels):
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            try:
                boxes = f[lev]['boxes'].value
            except KeyError:
                boxes = f[lev]['particles:boxes'].value
            dx = f[lev].attrs['dx']
            self.dds_list.append(dx * np.ones(3))

            if D == 1:
                self.dds_list[lev_index][1] = 1.0
                self.dds_list[lev_index][2] = 1.0

            if D == 2:
                self.dds_list[lev_index][2] = 1.0

            for level_id, box in enumerate(boxes):
                si = np.array([box['lo_%s' % ax] for ax in 'ijk'[:D]])
                ei = np.array([box['hi_%s' % ax] for ax in 'ijk'[:D]])

                if D == 1:
                    si = np.concatenate((si, [0.0, 0.0]))
                    ei = np.concatenate((ei, [0.0, 0.0]))

                if D == 2:
                    si = np.concatenate((si, [0.0]))
                    ei = np.concatenate((ei, [0.0]))

                pg = self.grid(len(grids),self,level=level_number,
                               start = si, stop = ei)
                grids.append(pg)
                grids[-1]._level_id = level_id
                self.grid_levels[i] = level_number
                self.grid_left_edge[i] = self.dds_list[lev_index]*si.astype(self.float_type)
                self.grid_right_edge[i] = self.dds_list[lev_index]*(ei.astype(self.float_type)+1)
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        self.grids = np.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g

    def _populate_grid_objects(self):
        self._reconstruct_parent_child()
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _reconstruct_parent_child(self):
        mask = np.empty(len(self.grids), dtype='int32')
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[i,:],
                                self.grid_right_edge[i,:],
                                self.grid_levels[i] + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            ids = np.where(mask.astype("bool")) # where is a tuple
            grid._children_ids = ids[0] + grid._id_offset 
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids): # Second pass
            for child in grid.Children:
                child._parent_id.append(i + grid._id_offset)


class ChomboDataset(Dataset):
    _index_class = ChomboHierarchy
    _field_info_class = ChomboFieldInfo

    def __init__(self, filename, dataset_type='chombo_hdf5',
                 storage_filename = None, ini_filename = None,
                 units_override=None, unit_system="cgs"):
        self.fluid_types += ("chombo",)
        self._handle = HDF5FileHandler(filename)
        self.dataset_type = dataset_type

        self.geometry = "cartesian"
        self.ini_filename = ini_filename
        self.fullplotdir = os.path.abspath(filename)
        Dataset.__init__(self, filename, self.dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename
        self.cosmological_simulation = False

        # These are parameters that I very much wish to get rid of.
        self.parameters["HydroMethod"] = 'chombo'
        self.parameters["DualEnergyFormalism"] = 0 
        self.parameters["EOSType"] = -1 # default

    def _set_code_unit_attributes(self):
        if not hasattr(self, 'length_unit'):
            mylog.warning("Setting code length unit to be 1.0 cm")
        if not hasattr(self, 'mass_unit'):
            mylog.warning("Setting code mass unit to be 1.0 g")
        if not hasattr(self, 'time_unit'):
            mylog.warning("Setting code time unit to be 1.0 s")
        setdefaultattr(self, 'length_unit', self.quan(1.0, "cm"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "g"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'magnetic_unit',
                       self.quan(np.sqrt(4.*np.pi), "gauss"))
        setdefaultattr(self, 'velocity_unit',
                       self.length_unit / self.time_unit)

    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.directory, default)
        return f

    def _parse_parameter_file(self):

        self.unique_identifier = \
                               int(os.stat(self.parameter_filename)[ST_CTIME])
        self.dimensionality = self._handle['Chombo_global/'].attrs['SpaceDim']
        self.domain_left_edge = self._calc_left_edge()
        self.domain_right_edge = self._calc_right_edge()
        self.domain_dimensions = self._calc_domain_dimensions()

        # if a lower-dimensional dataset, set up pseudo-3D stuff here.
        if self.dimensionality == 1:
            self.domain_left_edge = np.concatenate((self.domain_left_edge, [0.0, 0.0]))
            self.domain_right_edge = np.concatenate((self.domain_right_edge, [1.0, 1.0]))
            self.domain_dimensions = np.concatenate((self.domain_dimensions, [1, 1]))

        if self.dimensionality == 2:
            self.domain_left_edge = np.concatenate((self.domain_left_edge, [0.0]))
            self.domain_right_edge = np.concatenate((self.domain_right_edge, [1.0]))
            self.domain_dimensions = np.concatenate((self.domain_dimensions, [1]))

        self.refine_by = self._handle['/level_0'].attrs['ref_ratio']
        self._determine_periodic()
        self._determine_current_time()

    def _determine_current_time(self):
        # some datasets will not be time-dependent, and to make
        # matters worse, the simulation time is not always
        # stored in the same place in the hdf file! Make
        # sure we handle that here.
        try:
            self.current_time = self._handle.attrs['time']
        except KeyError:
            try:
                self.current_time = self._handle['/level_0'].attrs['time']
            except KeyError:
                self.current_time = 0.0

    def _determine_periodic(self):
        # we default to true unless the HDF5 file says otherwise
        is_periodic = np.array([True, True, True])
        for dir in range(self.dimensionality):
            try:
                is_periodic[dir] = self._handle['/level_0'].attrs['is_periodic_%d' % dir]
            except KeyError:
                is_periodic[dir] = True
        self.periodicity = tuple(is_periodic)

    def _calc_left_edge(self):
        fileh = self._handle
        dx0 = fileh['/level_0'].attrs['dx']
        D = self.dimensionality
        LE = dx0*((np.array(list(fileh['/level_0'].attrs['prob_domain'])))[0:D])
        return LE

    def _calc_right_edge(self):
        fileh = self._handle
        dx0 = fileh['/level_0'].attrs['dx']
        D = self.dimensionality
        RE = dx0*((np.array(list(fileh['/level_0'].attrs['prob_domain'])))[D:] + 1)
        return RE

    def _calc_domain_dimensions(self):
        fileh = self._handle
        D = self.dimensionality
        L_index = ((np.array(list(fileh['/level_0'].attrs['prob_domain'])))[0:D])
        R_index = ((np.array(list(fileh['/level_0'].attrs['prob_domain'])))[D:] + 1)
        return R_index - L_index

    @classmethod
    def _is_valid(self, *args, **kwargs):

        if not is_chombo_hdf5(args[0]):
            return False

        pluto_ini_file_exists = False
        orion2_ini_file_exists = False

        if isinstance(args[0], six.string_types): 
            dir_name = os.path.dirname(os.path.abspath(args[0]))
            pluto_ini_filename = os.path.join(dir_name, "pluto.ini")
            orion2_ini_filename = os.path.join(dir_name, "orion2.ini")
            pluto_ini_file_exists = os.path.isfile(pluto_ini_filename)
            orion2_ini_file_exists = os.path.isfile(orion2_ini_filename)

        if not (pluto_ini_file_exists or orion2_ini_file_exists):
            try:
                fileh = h5py.File(args[0],'r')
                valid = "Chombo_global" in fileh["/"]
                # ORION2 simulations should always have this:
                valid = valid and not ('CeilVA_mass' in fileh.attrs.keys())
                valid = valid and not ('Charm_global' in fileh.keys())
                fileh.close()
                return valid
            except:
                pass
        return False

    @parallel_root_only
    def print_key_parameters(self):
        for a in ["current_time", "domain_dimensions", "domain_left_edge",
                  "domain_right_edge"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)


class PlutoHierarchy(ChomboHierarchy):

    def __init__(self, ds, dataset_type="chombo_hdf5"):
        ChomboHierarchy.__init__(self, ds, dataset_type)

    def _parse_index(self):
        f = self._handle # shortcut
        self.max_level = f.attrs['num_levels'] - 1

        grids = []
        self.dds_list = []
        i = 0
        D = self.dataset.dimensionality
        for lev_index, lev in enumerate(self._levels):
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            try:
                boxes = f[lev]['boxes'].value
            except KeyError:
                boxes = f[lev]['particles:boxes'].value
            dx = f[lev].attrs['dx']
            self.dds_list.append(dx * np.ones(3))

            if D == 1:
                self.dds_list[lev_index][1] = 1.0
                self.dds_list[lev_index][2] = 1.0

            if D == 2:
                self.dds_list[lev_index][2] = 1.0

            for level_id, box in enumerate(boxes):
                si = np.array([box['lo_%s' % ax] for ax in 'ijk'[:D]])
                ei = np.array([box['hi_%s' % ax] for ax in 'ijk'[:D]])

                if D == 1:
                    si = np.concatenate((si, [0.0, 0.0]))
                    ei = np.concatenate((ei, [0.0, 0.0]))

                if D == 2:
                    si = np.concatenate((si, [0.0]))
                    ei = np.concatenate((ei, [0.0]))

                pg = self.grid(len(grids),self,level=level_number,
                               start = si, stop = ei)
                grids.append(pg)
                grids[-1]._level_id = level_id
                self.grid_levels[i] = level_number
                self.grid_left_edge[i] = self.dds_list[lev_index]*si.astype(self.float_type)+self.domain_left_edge.value
                self.grid_right_edge[i] = self.dds_list[lev_index]*(ei.astype(self.float_type)+1)+self.domain_left_edge.value
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        self.grids = np.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g


class PlutoDataset(ChomboDataset):

    _index_class = PlutoHierarchy
    _field_info_class = PlutoFieldInfo

    def __init__(self, filename, dataset_type='chombo_hdf5',
                 storage_filename = None, ini_filename = None,
                 units_override=None, unit_system="cgs"):

        ChomboDataset.__init__(self, filename, dataset_type, 
                               storage_filename, ini_filename,
                               units_override=units_override,
                               unit_system=unit_system)

    def _parse_parameter_file(self):
        """
        Check to see whether a 'pluto.ini' file
        exists in the plot file directory. If one does, attempt to parse it.
        Otherwise grab the dimensions from the hdf5 file.
        """

        pluto_ini_file_exists = False
        dir_name = os.path.dirname(os.path.abspath(self.fullplotdir))
        pluto_ini_filename = os.path.join(dir_name, "pluto.ini")
        pluto_ini_file_exists = os.path.isfile(pluto_ini_filename)

        self.unique_identifier = \
                               int(os.stat(self.parameter_filename)[ST_CTIME])
        self.dimensionality = self._handle['Chombo_global/'].attrs['SpaceDim']
        self.domain_dimensions = self._calc_domain_dimensions()
        self.refine_by = self._handle['/level_0'].attrs['ref_ratio']

        if pluto_ini_file_exists:
            lines=[line.strip() for line in open(pluto_ini_filename)]
            domain_left_edge = np.zeros(self.dimensionality)
            domain_right_edge = np.zeros(self.dimensionality)
            for il,ll in enumerate(lines[lines.index('[Grid]')+2:lines.index('[Grid]')+2+self.dimensionality]):
                domain_left_edge[il] = float(ll.split()[2])
                domain_right_edge[il] = float(ll.split()[-1])
            self.periodicity = [0]*3
            for il,ll in enumerate(lines[lines.index('[Boundary]')+2:lines.index('[Boundary]')+2+6:2]):
                self.periodicity[il] = (ll.split()[1] == 'periodic')
            self.periodicity = tuple(self.periodicity)
            for il,ll in enumerate(lines[lines.index('[Parameters]')+2:]):
                if (ll.split()[0] == 'GAMMA'):
                    self.gamma = float(ll.split()[1])
            self.domain_left_edge = domain_left_edge
            self.domain_right_edge = domain_right_edge
        else:
            self.domain_left_edge = self._calc_left_edge()
            self.domain_right_edge = self._calc_right_edge()
            self.periodicity = (True, True, True)

        # if a lower-dimensional dataset, set up pseudo-3D stuff here.
        if self.dimensionality == 1:
            self.domain_left_edge = np.concatenate((self.domain_left_edge, [0.0, 0.0]))
            self.domain_right_edge = np.concatenate((self.domain_right_edge, [1.0, 1.0]))
            self.domain_dimensions = np.concatenate((self.domain_dimensions, [1, 1]))

        if self.dimensionality == 2:
            self.domain_left_edge = np.concatenate((self.domain_left_edge, [0.0]))
            self.domain_right_edge = np.concatenate((self.domain_right_edge, [1.0]))
            self.domain_dimensions = np.concatenate((self.domain_dimensions, [1]))

        self._determine_current_time()

    @classmethod
    def _is_valid(self, *args, **kwargs):

        if not is_chombo_hdf5(args[0]):
            return False

        pluto_ini_file_exists = False

        if isinstance(args[0], six.string_types):
            dir_name = os.path.dirname(os.path.abspath(args[0]))
            pluto_ini_filename = os.path.join(dir_name, "pluto.ini")
            pluto_ini_file_exists = os.path.isfile(pluto_ini_filename)

        if pluto_ini_file_exists:
            return True

        return False


class Orion2Hierarchy(ChomboHierarchy):

    def __init__(self, ds, dataset_type="orion_chombo_native"):
        ChomboHierarchy.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):

        # look for fluid fields
        output_fields = []
        for key, val in self._handle.attrs.items():
            if key.startswith("component"):
                output_fields.append(val.decode("ascii"))
        self.field_list = [("chombo", c) for c in output_fields]

        # look for particle fields
        self.particle_filename = self.index_filename[:-4] + 'sink'
        if not os.path.exists(self.particle_filename):
            return
        pfield_list = [("io", str(c)) for c in self.io.particle_field_index.keys()]
        self.field_list.extend(pfield_list)

    def _read_particles(self):
        if not os.path.exists(self.particle_filename): return
        with open(self.particle_filename, 'r') as f:
            lines = f.readlines()
            self.num_stars = int(lines[0].strip().split(' ')[0])
            for num, line in enumerate(lines[1:]):
                particle_position_x = float(line.split(' ')[1])
                particle_position_y = float(line.split(' ')[2])
                particle_position_z = float(line.split(' ')[3])
                coord = [particle_position_x, particle_position_y, particle_position_z]
                # for each particle, determine which grids contain it
                # copied from object_finding_mixin.py
                mask = np.ones(self.num_grids)
                for i in range(len(coord)):
                    np.choose(np.greater(self.grid_left_edge.d[:,i],coord[i]), (mask,0), mask)
                    np.choose(np.greater(self.grid_right_edge.d[:,i],coord[i]), (0,mask), mask)
                ind = np.where(mask == 1)
                selected_grids = self.grids[ind]
                # in orion, particles always live on the finest level.
                # so, we want to assign the particle to the finest of
                # the grids we just found
                if len(selected_grids) != 0:
                    grid = sorted(selected_grids, key=lambda grid: grid.Level)[-1]
                    ind = np.where(self.grids == grid)[0][0]
                    self.grid_particle_count[ind] += 1
                    self.grids[ind].NumberOfParticles += 1

                    # store the position in the *.sink file for fast access.
                    try:
                        self.grids[ind]._particle_line_numbers.append(num + 1)
                    except AttributeError:
                        self.grids[ind]._particle_line_numbers = [num + 1]


class Orion2Dataset(ChomboDataset):

    _index_class = Orion2Hierarchy
    _field_info_class = Orion2FieldInfo

    def __init__(self, filename, dataset_type='orion_chombo_native',
                 storage_filename = None, ini_filename = None,
                 units_override=None):

        ChomboDataset.__init__(self, filename, dataset_type,
                               storage_filename, ini_filename,
                               units_override=units_override)

    def _parse_parameter_file(self):
        """
        Check to see whether an 'orion2.ini' file
        exists in the plot file directory. If one does, attempt to parse it.
        Otherwise grab the dimensions from the hdf5 file.
        """

        orion2_ini_file_exists = False
        dir_name = os.path.dirname(os.path.abspath(self.fullplotdir))
        orion2_ini_filename = os.path.join(dir_name, "orion2.ini")
        orion2_ini_file_exists = os.path.isfile(orion2_ini_filename)

        if orion2_ini_file_exists: self._parse_inputs_file('orion2.ini')
        self.unique_identifier = \
                               int(os.stat(self.parameter_filename)[ST_CTIME])
        self.dimensionality = 3
        self.domain_left_edge = self._calc_left_edge()
        self.domain_right_edge = self._calc_right_edge()
        self.domain_dimensions = self._calc_domain_dimensions()
        self.refine_by = self._handle['/level_0'].attrs['ref_ratio']
        self._determine_periodic()
        self._determine_current_time()

    def _parse_inputs_file(self, ini_filename):
        self.fullplotdir = os.path.abspath(self.parameter_filename)
        self.ini_filename = self._localize( \
            self.ini_filename, ini_filename)
        self.unique_identifier = \
                               int(os.stat(self.parameter_filename)[ST_CTIME])
        lines = open(self.ini_filename).readlines()
        # read the file line by line, storing important parameters
        for lineI, line in enumerate(lines):
            try:
                param, sep, vals = line.partition('=')
                if not sep:
                    # No = sign present, so split by space instead
                    param, sep, vals = line.partition(' ')
                param = param.strip()
                vals = vals.strip()
                if not param:  # skip blank lines
                    continue
                if param[0] == '#':  # skip comment lines
                    continue
                if param[0] == '[':  # skip stanza headers
                    continue
                vals = vals.partition("#")[0] # strip trailing comments
                try:
                    self.parameters[param] = np.int64(vals)
                except ValueError:
                    try:
                        self.parameters[param] = np.float64(vals)
                    except ValueError:
                        self.parameters[param] = vals
            except ValueError:
                mylog.error("ValueError: '%s'", line)
            if param == "GAMMA":
                self.gamma = np.float64(vals)

    @classmethod
    def _is_valid(self, *args, **kwargs):

        if not is_chombo_hdf5(args[0]):
            return False

        pluto_ini_file_exists = False
        orion2_ini_file_exists = False

        if isinstance(args[0], string_types):
            dir_name = os.path.dirname(os.path.abspath(args[0]))
            pluto_ini_filename = os.path.join(dir_name, "pluto.ini")
            orion2_ini_filename = os.path.join(dir_name, "orion2.ini")
            pluto_ini_file_exists = os.path.isfile(pluto_ini_filename)
            orion2_ini_file_exists = os.path.isfile(orion2_ini_filename)

        if orion2_ini_file_exists:
            return True

        if not pluto_ini_file_exists:
            try:
                fileh = h5py.File(args[0],'r')
                valid = 'CeilVA_mass' in fileh.attrs.keys()
                valid = "Chombo_global" in fileh["/"] and "Charm_global" not in fileh["/"]
                valid = valid and 'CeilVA_mass' in fileh.attrs.keys()
                fileh.close()
                return valid
            except:
                pass
        return False


class ChomboPICHierarchy(ChomboHierarchy):

    def __init__(self, ds, dataset_type="chombo_hdf5"):
        ChomboHierarchy.__init__(self, ds, dataset_type)


class ChomboPICDataset(ChomboDataset):

    _index_class = ChomboPICHierarchy
    _field_info_class = ChomboPICFieldInfo3D

    def __init__(self, filename, dataset_type='chombo_hdf5',
                 storage_filename=None, ini_filename=None,
                 units_override=None):

        ChomboDataset.__init__(self, filename, dataset_type,
                               storage_filename, ini_filename,
                               units_override=units_override)

        if self.dimensionality == 1:
            self._field_info_class = ChomboPICFieldInfo1D

        if self.dimensionality == 2:
            self._field_info_class = ChomboPICFieldInfo2D

    @classmethod
    def _is_valid(self, *args, **kwargs):

        if not is_chombo_hdf5(args[0]):
            return False

        pluto_ini_file_exists = False
        orion2_ini_file_exists = False

        if isinstance(args[0], six.string_types):
            dir_name = os.path.dirname(os.path.abspath(args[0]))
            pluto_ini_filename = os.path.join(dir_name, "pluto.ini")
            orion2_ini_filename = os.path.join(dir_name, "orion2.ini")
            pluto_ini_file_exists = os.path.isfile(pluto_ini_filename)
            orion2_ini_file_exists = os.path.isfile(orion2_ini_filename)

        if orion2_ini_file_exists:
            return False

        if pluto_ini_file_exists:
            return False

        try:
            fileh = h5py.File(args[0],'r')
            valid = "Charm_global" in fileh["/"]
            fileh.close()
            return valid
        except:
            pass
        return False
