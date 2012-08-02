"""
Data structures for Castro.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2010 J. S. Oishi.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import re
import os
import weakref
import itertools
from collections import defaultdict
from string import strip, rstrip
from stat import ST_CTIME

import numpy as na

from yt.funcs import *
from yt.data_objects.field_info_container import FieldInfoContainer, NullFunc
from yt.data_objects.grid_patch import AMRGridPatch
from yt.data_objects.hierarchy import AMRHierarchy
from yt.data_objects.static_output import StaticOutput
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.lib import get_box_grids_level

from .definitions import \
    castro2enzoDict, \
    parameterDict, \
    yt2castroFieldsDict, \
    castro_FAB_header_pattern, \
    castro_particle_field_names, \
    boxlib_bool_to_int
from .fields import \
    CastroFieldInfo, \
    KnownCastroFields, \
    add_castro_field


class CastroGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, LeftEdge, RightEdge, index, level, filename, offset,
                 dimensions, start, stop, paranoia=False, **kwargs):
        super(CastroGrid, self).__init__(self, index, **kwargs)
        self.filename = filename
        self._offset = offset
        self._paranoid = paranoia  # TODO: Factor this behavior out in tests

        ### TODO: error check this (test)
        self.ActiveDimensions = (dimensions.copy()).astype('int32')#.transpose()
        self.start_index = start.copy()#.transpose()
        self.stop_index = stop.copy()#.transpose()
        self.LeftEdge  = LeftEdge.copy()
        self.RightEdge = RightEdge.copy()
        self.index = index
        self.level = level

    def get_global_startindex(self):
        return self.start_index

    def _prepare_grid(self):
        """ Copies all the appropriate attributes from the hierarchy. """
        # This is definitely the slowest part of generating the hierarchy.
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8

        h = self.hierarchy # cache it
        #self.StartIndices = h.gridStartIndices[self.id]
        #self.EndIndices = h.gridEndIndices[self.id]
        h.grid_levels[self.id,0] = self.Level
        h.grid_left_edge[self.id,:] = self.LeftEdge[:]
        h.grid_right_edge[self.id,:] = self.RightEdge[:]
        #self.Time = h.gridTimes[self.id,0]
        #self.NumberOfParticles = h.gridNumberOfParticles[self.id,0]
        self.field_indexes = h.field_indexes
        self.Children = h.gridTree[self.id]
        pIDs = h.gridReverseTree[self.id]

        if len(pIDs) > 0:
            self.Parent = [weakref.proxy(h.grids[pID]) for pID in pIDs]
        else:
            self.Parent = None

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if self.Parent is not None:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)

        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "CastroGrid_%04i" % (self.id)

class CastroHierarchy(AMRHierarchy):
    grid = CastroGrid

    def __init__(self, pf, data_style='castro_native'):
        super(CastroHierarchy, self).__init__(self, pf, self.data_style)

        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        header_filename = os.path.join(pf.fullplotdir, 'Header')
        self.directory = pf.fullpath
        self.data_style = data_style

        # This also sets up the grid objects
        self.read_global_header(header_filename,
                                self.parameter_file.paranoid_read) 
        self.read_particle_header()
        self._cache_endianness(self.levels[-1].grids[-1])
        self._setup_data_io()
        self._setup_field_list()
        self._populate_hierarchy()

    def read_global_header(self, filename, paranoid_read):
        """ Read the global header file for an Castro plotfile output. """
        counter = 0
        header_file = open(filename, 'r')
        self._global_header_lines = header_file.readlines()

        # parse the file
        self.castro_version = self._global_header_lines[0].rstrip()
        self.n_fields = int(self._global_header_lines[1])

        counter = self.n_fields + 2
        self.field_list = []
        for i, line in enumerate(self._global_header_lines[2:counter]):
            self.field_list.append(line.rstrip())

        # this is unused...eliminate it?
        #for f in self.field_indexes:
        #    self.field_list.append(castro2ytFieldsDict.get(f, f))

        self.dimension = int(self._global_header_lines[counter])
        if self.dimension != 3:
            raise RunTimeError("Castro must be in 3D to use yt.")

        counter += 1
        self.Time = float(self._global_header_lines[counter])
        counter += 1
        self.finest_grid_level = int(self._global_header_lines[counter])
        self.n_levels = self.finest_grid_level + 1
        counter += 1

        # quantities with _unnecessary are also stored in the inputs
        # file and are not needed.  they are read in and stored in
        # case in the future we want to enable a "backwards" way of
        # taking the data out of the Header file and using it to fill
        # in in the case of a missing inputs file
        self.domainLeftEdge_unnecessary = na.array(map(float, self._global_header_lines[counter].split()))
        counter += 1
        self.domainRightEdge_unnecessary = na.array(map(float, self._global_header_lines[counter].split()))
        counter += 1
        self.refinementFactor_unnecessary = self._global_header_lines[counter].split()
        #na.array(map(int, self._global_header_lines[counter].split()))
        counter += 1
        self.globalIndexSpace_unnecessary = self._global_header_lines[counter]
        #domain_re.search(self._global_header_lines[counter]).groups()
        counter += 1
        self.timestepsPerLevel_unnecessary = self._global_header_lines[counter]
        counter += 1

        self.dx = na.zeros((self.n_levels, 3))
        for i, line in enumerate(self.__global_header_lines[counter:counter+self.n_levels]):
            self.dx[i] = na.array(map(float, line.split()))
        counter += self.n_levels
        self.geometry = int(self._global_header_lines[counter])
        if self.geometry != 0:
            raise RunTimeError("yt only supports cartesian coordinates.")
        counter += 1

        # this is just to debug. eventually it should go away.
        linebreak = int(self._global_header_lines[counter])
        if linebreak != 0:
            raise RunTimeError("INTERNAL ERROR! Header is unexpected size")
        counter += 1

        # Each level is one group with ngrids on it. each grid has 3 lines of 2 reals
        # BoxLib madness
        self.levels = []
        grid_counter = 0
        file_finder_pattern = r"FabOnDisk: (\w+_D_[0-9]{4}) (\d+)\n"
        re_file_finder = re.compile(file_finder_pattern)
        dim_finder_pattern = r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)\n"
        re_dim_finder = re.compile(dim_finder_pattern)
        data_files_pattern = r"Level_[\d]/"
        data_files_finder = re.compile(data_files_pattern)

        for level in range(0, self.n_levels):
            tmp = self._global_header_lines[counter].split()
            # Should this be grid_time or level_time??
            lev, ngrids, grid_time = int(tmp[0]), int(tmp[1]), float(tmp[2])
            counter += 1
            nsteps = int(self._global_header_lines[counter])
            counter += 1
            self.levels.append(CastroLevel(lev, ngrids))
            # Open level header, extract file names and offsets for each grid.
            # Read slightly out of order here: at the end of the lo, hi pairs
            # for x, y, z is a *list* of files types in the Level directory. 
            # Each type has Header and a number of data files
            # (one per processor)
            tmp_offset = counter + 3*ngrids
            nfiles = 0
            key_off = 0
            files =   {} # dict(map(lambda a: (a,[]), self.field_list))
            offsets = {} # dict(map(lambda a: (a,[]), self.field_list))

            while (nfiles + tmp_offset < len(self._global_header_lines) and
                   data_files_finder.match(self._global_header_lines[nfiles+tmp_offset])):
                filen = os.path.join(self.parameter_file.fullplotdir,
                                     self._global_header_lines[nfiles+tmp_offset].strip())
                # open each "_H" header file, and get the number of
                # components within it
                level_header_file = open(filen+'_H','r').read()
                start_stop_index = re_dim_finder.findall(level_header_file) # just take the last one
                grid_file_offset = re_file_finder.findall(level_header_file)
                ncomp_this_file = int(level_header_file.split('\n')[2])

                for i in range(ncomp_this_file):
                    key = self.field_list[i+key_off]
                    f, o = zip(*grid_file_offset)
                    files[key] = f
                    offsets[key] = o
                    self.field_indexes[key] = i

                key_off += ncomp_this_file
                nfiles += 1

            # convert dict of lists to list of dicts
            fn = []
            off = []
            lead_path = os.path.join(self.parameter_file.fullplotdir,
                                     'Level_%i' % level)
            for i in range(ngrids):
                fi = [os.path.join(lead_path, files[key][i]) for key in self.field_list]
                of = [int(offsets[key][i]) for key in self.field_list]
                fn.append(dict(zip(self.field_list, fi)))
                off.append(dict(zip(self.field_list, of)))

            for grid in range(0, ngrids):
                gfn = fn[grid]  # filename of file containing this grid
                gfo = off[grid] # offset within that file
                xlo, xhi = map(float, self._global_header_lines[counter].split())
                counter += 1
                ylo, yhi = map(float, self._global_header_lines[counter].split())
                counter += 1
                zlo, zhi = map(float, self._global_header_lines[counter].split())
                counter += 1
                lo = na.array([xlo, ylo, zlo])
                hi = na.array([xhi, yhi, zhi])
                dims, start, stop = self._calculate_grid_dimensions(start_stop_index[grid])
                self.levels[-1].grids.append(self.grid(lo, hi, grid_counter,
                                                       level, gfn, gfo, dims,
                                                       start, stop,
                                                       paranoia=paranoid_read,  ### TODO: at least the code isn't schizophrenic paranoid
                                                       hierarchy=self))
                grid_counter += 1   # this is global, and shouldn't be reset
                                    # for each level

            # already read the filenames above...
            counter += nfiles
            self.num_grids = grid_counter
            self.float_type = 'float64'

        self.maxLevel = self.n_levels - 1
        self.max_level = self.n_levels - 1
        header_file.close()

    def read_particle_header(self):
        # We need to get particle offsets and particle counts
        if not self.parameter_file.use_particles:
            self.pgrid_info = na.zeros((self.num_grids, 3), dtype='int64')
            return

        self.field_list += castro_particle_field_names[:]
        header = open(os.path.join(self.parameter_file.fullplotdir, "DM",
                                   "Header"))
        version = header.readline()
        ndim = header.readline()
        nfields = header.readline()
        ntotalpart = int(header.readline())
        dummy = header.readline() # nextid
        maxlevel = int(header.readline()) # max level

        # Skip over how many grids on each level; this is degenerate
        for i in range(maxlevel+1): dummy = header.readline()
        grid_info = na.fromiter((int(i)
                                 for line in header.readlines()
                                 for i in line.split()),
                                dtype='int64',
                                count=3*self.num_grids).reshape((self.num_grids, 3))
        self.pgrid_info = grid_info

    def _cache_endianness(self, test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a test grid
        and assuming that all grids have the same endianness. This is a pretty
        safe assumption since Castro uses one file per processor, and if you're
        running on a cluster with different endian processors, then you're on
        your own!

        """
        # open the test file and grab the header
        in_file = open(os.path.expanduser(test_grid.filename[self.field_list[0]]), 'rb')
        header = in_file.readline()
        in_file.close()
        header.strip()
        # Parse it. The pattern is in castro.definitions.py
        header_re = re.compile(castro_FAB_header_pattern)
        bytes_per_real, endian, start, stop, centerType, n_components = header_re.search(header).groups()
        self._bytes_per_real = int(bytes_per_real)
        if self._bytes_per_real == int(endian[0]):
            dtype = '<'
        elif self._bytes_per_real == int(endian[-1]):
            dtype = '>'
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        dtype += ('f%i' % self._bytes_per_real) # always a floating point
        self._dtype = dtype

    def _calculate_grid_dimensions(self, start_stop):
        start = na.array(map(int, start_stop[0].split(',')))
        stop = na.array(map(int, start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension, start, stop

    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")

        self.grids = na.concatenate([level.grids for level in self.levels])
        basedir = self.parameter_file.fullplotdir

        for g, pg in itertools.izip(self.grids, self.pgrid_info):
            g.particle_filename = os.path.join(
                basedir, "DM", "Level_%s" % (g.Level), "DATA_%04i" % pg[0])
            g.NumberOfParticles = pg[1]
            g._particle_offset = pg[2]

        self.grid_particle_count[:,0] = self.pgrid_info[:,1]
        del self.pgrid_info

        gls = na.concatenate([level.ngrids * [level.level] for level in self.levels])
        self.grid_levels[:] = gls.reshape((self.num_grids,1))
        grid_dcs = na.concatenate([level.ngrids * [self.dx[level.level]]
                                  for level in self.levels], axis=0)

        self.grid_dxs = grid_dcs[:,0].reshape((self.num_grids,1))
        self.grid_dys = grid_dcs[:,1].reshape((self.num_grids,1))
        self.grid_dzs = grid_dcs[:,2].reshape((self.num_grids,1))

        left_edges = []
        right_edges = []
        dims = []
        for level in self.levels:
            left_edges += [g.LeftEdge for g in level.grids]
            right_edges += [g.RightEdge for g in level.grids]
            dims += [g.ActiveDimensions for g in level.grids]

        self.grid_left_edge = na.array(left_edges)
        self.grid_right_edge = na.array(right_edges)
        self.grid_dimensions = na.array(dims)
        self.gridReverseTree = [] * self.num_grids
        self.gridReverseTree = [ [] for i in range(self.num_grids)]
        self.gridTree = [ [] for i in range(self.num_grids)]

        mylog.debug("Done creating grid objects")

    def _populate_hierarchy(self):
        self._setup_grid_tree()
        #self._setup_grid_corners()

        for i, grid in enumerate(self.grids):
            if (i % 1e4) == 0:
                mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)

            grid._prepare_grid()
            grid._setup_dx()

    def _setup_grid_tree(self):
        mask = na.empty(self.grids.size, dtype='int32')
        for i, grid in enumerate(self.grids):
            get_box_grids_level(grid.LeftEdge, grid.RightEdge, grid.Level + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            children = self.grids[mask.astype("bool")]
            #assert(len(children) == len(self._get_grid_children(grid)))
            for child in children:
                self.gridReverseTree[child.id].append(i)
                self.gridTree[i].append(weakref.proxy(child))

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        dd["field_indexes"] = self.field_indexes
        AMRHierarchy._setup_classes(self, dd)
        #self._add_object_class('grid', "CastroGrid", CastroGridBase, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.grid_levels == (grid.Level+1)).flat)
        return self.grids[mask]

    def _setup_field_list(self):
        self.derived_field_list = []

        for field in self.field_info:
            try:
                fd = self.field_info[field].get_dependencies(pf=self.parameter_file)
            except:
                continue

            available = na.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)

        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

        if self.parameter_file.use_particles:
            # We know which particle fields will exist -- pending further
            # changes in the future.
            for field in castro_particle_field_names:
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
                # Note that we call add_castro_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D, 2D
                # and 3D fields.
                self.pf.field_info.add_castro_field(
                        field, lambda a, b: None,
                        convert_function=cf, take_log=False,
                        particle_type=True)

    ### TODO: check if this can be removed completely
    def _count_grids(self):
        """
        this is already provided in ???

        """
        pass

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = na.ones((self.num_grids,3), 'int32')
        self.grid_left_edge = na.zeros((self.num_grids,3), self.float_type)
        self.grid_right_edge = na.ones((self.num_grids,3), self.float_type)
        self.grid_levels = na.zeros((self.num_grids,1), 'int32')
        self.grid_particle_count = na.zeros((self.num_grids,1), 'int32')

    def _parse_hierarchy(self):
        pass

    def _detect_fields(self):
        pass

    def _setup_derived_fields(self):
        pass

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in AMRHierarchy.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}

class CastroLevel:
    def __init__(self, level, ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []

class CastroStaticOutput(StaticOutput):
    """
    This class is a stripped down class that simply reads and parses *filename*,
    without looking at the Castro hierarchy.

    """
    _hierarchy_class = CastroHierarchy
    _fieldinfo_fallback = CastroFieldInfo
    _fieldinfo_known = KnownCastroFields

    def __init__(self, plotname, paramFilename=None, fparamFilename=None,
                 data_style='castro_native', paranoia=False,
                 storage_filename = None):
        """
        Need to override for Castro file structure.

        the paramfile is usually called "inputs"
        and there may be a fortran inputs file usually called "probin"
        plotname here will be a directory name
        as per BoxLib, data_style will be one of
         * Native
         * IEEE (not implemented in yt)
         * ASCII (not implemented in yt)

        """
        super(CastroStaticOutput, self).__init__(self, plotname.rstrip("/"),
                                                 data_style='castro_native')
        self.storage_filename = storage_filename
        self.paranoid_read = paranoia
        self.parameter_filename = paramFilename
        self.fparameter_filename = fparamFilename
        self.__ipfn = paramFilename

        self.fparameters = {}

        # These should maybe not be hardcoded?
        ### TODO: this.
        self.parameters["HydroMethod"] = 'castro' # always PPM DE
        self.parameters["Time"] = 1.0 # default unit is 1...
        self.parameters["DualEnergyFormalism"] = 0 # always off.
        self.parameters["EOSType"] = -1 # default

        if self.fparameters.has_key("mu"):
            self.parameters["mu"] = self.fparameters["mu"]

    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.directory, default)
        return f

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        pname = args[0].rstrip("/")
        dn = os.path.dirname(pname)
        if len(args) > 1:
            kwargs['paramFilename'] = args[1]

        pfname = kwargs.get("paramFilename", os.path.join(dn, "inputs"))

        # We check for the job_info file's existence because this is currently
        # what distinguishes Castro data from MAESTRO data.
        ### ^ that is nuts
        pfn = os.path.join(pfname)
        if not os.path.exists(pfn):
            return False
        castro = any(("castro." in line for line in open(pfn)))
        nyx = any(("nyx." in line for line in open(pfn)))
        castro = castro and (not nyx) # it's only castro if it's not nyx
        maestro = os.path.exists(os.path.join(pname, "job_info"))
        orion = (not castro) and (not maestro)
        return castro

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various dictionaries.

        """
        # Boxlib madness
        self.fullplotdir = os.path.abspath(self.parameter_filename)
        self._parse_header_file()
        self.parameter_filename = self._localize(self.__ipfn, 'inputs')
        self.fparameter_filename = self._localize(self.fparameter_filename, 'probin')
        if os.path.isfile(self.fparameter_filename):
            self._parse_fparameter_file()
            for param in self.fparameters:
                if castro2enzoDict.has_key(param):
                    self.parameters[castro2enzoDict[param]] = self.fparameters[param]

        # Let's read the file
        self.unique_identifier = int(os.stat(self.parameter_filename)[ST_CTIME])
        lines = open(self.parameter_filename).readlines()
        self.use_particles = False

        for line in lines:
            if line.find("#") >= 1: # Keep the commented lines...
                line = line[:line.find("#")]
            line = line.strip().rstrip()
            if len(line) < 2 or line.find("#") == 0: # ...but skip comments
                continue

            try:
                param, vals = map(strip, map(rstrip, line.split("=")))
            except ValueError:
                mylog.error("ValueError: '%s'", line)

            if castro2enzoDict.has_key(param):
                paramName = castro2enzoDict[param]
                t = map(parameterDict[paramName], vals.split())
                if len(t) == 1:
                    self.parameters[paramName] = t[0]
                else:
                    if paramName == "RefineBy":
                        self.parameters[paramName] = t[0]
                    else:
                        self.parameters[paramName] = t
            elif param.startswith("geometry.prob_hi"):
                self.domain_right_edge = na.array([float(i) for i in vals.split()])
            elif param.startswith("geometry.prob_lo"):
                self.domain_left_edge = na.array([float(i) for i in vals.split()])
            elif param.startswith("particles.write_in_plotfile"):
                self.use_particles = boxlib_bool_to_int(vals)

        self.parameters["TopGridRank"] = len(self.parameters["TopGridDimensions"])
        self.dimensionality = self.parameters["TopGridRank"]
        self.domain_dimensions = self.parameters["TopGridDimensions"]
        self.refine_by = self.parameters.get("RefineBy", 2)

        if (self.parameters.has_key("ComovingCoordinates") and
            bool(self.parameters["ComovingCoordinates"])):
            self.cosmological_simulation = 1
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]

            # Stupid that we have to read a separate file for this :/
            a_file = open(os.path.join(self.fullplotdir, "comoving_a"))
            line = a_file.readline().strip()
            a_file.close()

            self.parameters["CosmologyCurrentRedshift"] = 1 / float(line) - 1
            self.cosmological_scale_factor = float(line)
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
        else:
            ### TODO: make these defaults automatic
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _parse_fparameter_file(self):
        """
        Parses the fortran parameter file for Castro. Most of this will be
        useless, but this is where it keeps mu = mass per particle/m_hydrogen.

        """
        lines = open(self.fparameter_filename).readlines()
        for line in lines:
            if line.count("=") == 1:
                param, vals = map(strip, map(rstrip, line.split("=")))
                if vals.count("'") == 0:
                    t = map(float, [a.replace('D','e').replace('d','e') for a in vals.split()]) # all are floating point.
                else:
                    t = vals.split()
                if len(t) == 1:
                    self.fparameters[param] = t[0]
                else:
                    self.fparameters[param] = t

    def _parse_header_file(self):
        """
        Parses the BoxLib header file to get any parameters stored there.
        Hierarchy information is read out of this file in CastroHierarchy. 

        Currently, only Time is read here.

        """
        header_file = open(os.path.join(self.fullplotdir, "Header"))
        lines = header_file.readlines()
        header_file.close()
        n_fields = int(lines[1])
        self.current_time = float(lines[3 + n_fields])

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file.

        """
        self.units = {}
        self.time_units = {}

        if len(self.parameters) == 0:
            self._parse_parameter_file()

        if self.cosmological_simulation:
            cf = 1e5 * self.cosmological_scale_factor   # Where does the 1e5 come from?
            for ax in 'xyz':
                self.units['particle_velocity_%s' % ax] = cf
            self.units['particle_mass'] = 1.989e33  ### TODO: Make a global solar mass def

        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units. Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = seconds / sec_conversion[unit]
        for key in yt2castroFieldsDict:
            self.conversion_factors[key] = 1.0
        for key in castro_particle_field_names:
            self.conversion_factors[key] = 1.0

    def _setup_nounits_units(self):
        z = 0
