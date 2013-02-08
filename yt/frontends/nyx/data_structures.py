"""
Data structures for Nyx.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Casey W. Stark, J. S. Oishi, Matthew Turk.  All Rights
  Reserved.

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

from collections import defaultdict
import itertools
import os
import re
from stat import ST_CTIME
from string import strip, rstrip
import weakref

import numpy as np

from yt.funcs import *
from yt.data_objects.grid_patch import AMRGridPatch
from yt.geometry.grid_geometry_handler import GridGeometryHandler
from yt.data_objects.static_output import StaticOutput
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.utilities.lib import get_box_grids_level
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion

from .definitions import parameter_type_dict, nyx_to_enzo_dict, \
                         fab_header_pattern, nyx_particle_field_names
from .utils import boxlib_bool_to_int
from .fields import NyxFieldInfo, add_nyx_field, KnownNyxFields


class NyxGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, left_edge, right_edge, index, level, filename, offset,
                 dimensions, start, stop, **kwargs):
        """ Build a grid that understands Nyx's boxlib format. """
        # passes index as the patch ``id``
        AMRGridPatch.__init__(self, index, **kwargs)
        self.filename = filename
        self._offset = offset

        # @todo: enzo-isms.
        # Must copy to avoid refs, in case the user alters the args later.
        self.ActiveDimensions = (dimensions.copy()).astype('int32')
        self.start_index = start.copy()
        self.stop_index = stop.copy()
        self.LeftEdge  = left_edge.copy()
        self.RightEdge = right_edge.copy()
        self.index = index
        self.Level = level

    def get_global_startindex(self):
        return self.start_index

    def _prepare_grid(self):
        """ Copies all the appropriate attributes from the hierarchy. """
        h = self.hierarchy  # alias
        h.grid_levels[self.id, 0] = self.Level
        h.grid_left_edge[self.id,:] = self.LeftEdge[:]
        h.grid_right_edge[self.id,:] = self.RightEdge[:]

        # Might still work
        #self.Time = h.gridTimes[self.id,0]
        #self.NumberOfParticles = h.gridNumberOfParticles[self.id,0]

        # @todo: enzo-isms
        self.field_indexes = h.field_indexes
        self.Children = h.gridTree[self.id]
        pIDs = h.gridReverseTree[self.id]

        if len(pIDs) > 0:
            self.Parent = [weakref.proxy(h.grids[pID]) for pID in pIDs]
        else:
            # must be root grid
            self.Parent = None

    def _setup_dx(self):
        # has already been read in and stored in hierarchy
        dx = self.hierarchy.grid_dxs[self.index][0]
        dy = self.hierarchy.grid_dys[self.index][0]
        dz = self.hierarchy.grid_dzs[self.index][0]
        self.dds = np.array([dx, dy, dz])
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "NyxGrid_%04i" % (self.id)

class NyxHierarchy(GridGeometryHandler):
    grid = NyxGrid

    def __init__(self, pf, data_style="nyx_native"):
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        self.directory = pf.path
        header_path = os.path.join(self.directory, "Header")  # make a kwarg?

        self.data_style = data_style
        #self._setup_classes()

        # This also sets up the grid objects
        self.read_global_header(header_path)
        self.read_particle_header()
        self.__cache_endianness(self.levels[-1].grids[-1])

        GridGeometryHandler.__init__(self, pf, self.data_style)
        self._setup_data_io()
        self._setup_field_list()
        self._populate_hierarchy()

    def read_global_header(self, header_path):
        """ Read the global header file for an Nyx plotfile output. """
        counter = 0
        header_file = open(header_path, 'r')
        self._global_header_lines = header_file.readlines()

        # parse the file
        self.nyx_pf_version = self._global_header_lines[0].rstrip()
        self.n_fields = int(self._global_header_lines[1])

        # why the 2?
        counter = self.n_fields + 2
        self.field_list = []
        for i, line in enumerate(self._global_header_lines[2:counter]):
            self.field_list.append(line.rstrip())

        # figure out dimensions and make sure it's 3D
        self.dimension = int(self._global_header_lines[counter])
        if self.dimension != 3:
            raise RunTimeError("Current data is %iD. yt only supports Nyx data in 3D" % self.dimension)

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
        self.domainLeftEdge_unnecessary = np.array(map(float, self._global_header_lines[counter].split()))
        counter += 1
        self.domainRightEdge_unnecessary = np.array(map(float, self._global_header_lines[counter].split()))
        counter += 1
        self.refinementFactor_unnecessary = self._global_header_lines[counter].split() #np.array(map(int, self._global_header_lines[counter].split()))
        counter += 1
        self.globalIndexSpace_unnecessary = self._global_header_lines[counter]
        counter += 1
        self.timestepsPerLevel_unnecessary = self._global_header_lines[counter]
        counter += 1

        self.dx = np.zeros((self.n_levels, 3))
        for i, line in enumerate(self._global_header_lines[counter:counter + self.n_levels]):
            self.dx[i] = np.array(map(float, line.split()))
        counter += self.n_levels
        self.geometry = int(self._global_header_lines[counter])
        if self.geometry != 0:
            raise RunTimeError("yt only supports cartesian coordinates.")
        counter += 1

        # @todo: this is just to debug. eventually it should go away.
        linebreak = int(self._global_header_lines[counter])
        if linebreak != 0:
            raise RunTimeError("INTERNAL ERROR! This should be a zero.")
        counter += 1

        # each level is one group with ngrids on it. each grid has 3 lines of 2
        # reals
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
            # should this be grid_time or level_time??
            lev, ngrids, grid_time = int(tmp[0]), int(tmp[1]), float(tmp[2])
            counter += 1
            nsteps = int(self._global_header_lines[counter])
            counter += 1
            self.levels.append(NyxLevel(lev, ngrids))

            # Open level header, extract file names and offsets for each grid
            # read slightly out of order here: at the end of the lo, hi pairs
            # for x, y, z is a *list* of files types in the Level directory.
            # Each type has Header and a number of data files (one per
            # processor)
            tmp_offset = counter + 3 * ngrids
            nfiles = 0
            key_off = 0
            files = {}
            offsets = {}
            while nfiles + tmp_offset < len(self._global_header_lines) \
                  and data_files_finder.match(self._global_header_lines[nfiles + tmp_offset]):
                filen = os.path.join(self.parameter_file.path,
                                     self._global_header_lines[nfiles + tmp_offset].strip())
                # open each "_H" header file, and get the number of
                # components within it
                level_header_file = open(filen + '_H', 'r').read()
                start_stop_index = re_dim_finder.findall(level_header_file) # just take the last one
                grid_file_offset = re_file_finder.findall(level_header_file)
                ncomp_this_file = int(level_header_file.split('\n')[2])

                for i in range(ncomp_this_file):
                    key = self.field_list[i + key_off]
                    f, o = zip(*grid_file_offset)
                    files[key] = f
                    offsets[key] = o
                    self.field_indexes[key] = i

                key_off += ncomp_this_file
                nfiles += 1

            # convert dict of lists to list of dicts
            fn = []
            off = []
            lead_path = os.path.join(self.parameter_file.path,
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
                lo = np.array([xlo, ylo, zlo])
                hi = np.array([xhi, yhi, zhi])
                dims, start, stop = self.__calculate_grid_dimensions(start_stop_index[grid])
                self.levels[-1].grids.append(self.grid(lo, hi, grid_counter,
                                             level, gfn, gfo, dims, start, stop,
                                             hierarchy=self))
                grid_counter += 1 # this is global, and shouldn't be reset
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
            self.pgrid_info = np.zeros((self.num_grids, 3), dtype='int64')
            return
        self.field_list += nyx_particle_field_names[:]
        header = open(os.path.join(self.parameter_file.path, "DM", "Header"))
        version = header.readline()
        ndim = header.readline()
        nfields = header.readline()
        ntotalpart = int(header.readline())
        dummy = header.readline() # nextid
        maxlevel = int(header.readline()) # max level

        # Skip over how many grids on each level; this is degenerate
        for i in range(maxlevel + 1):dummy = header.readline()

        grid_info = np.fromiter((int(i) for line in header.readlines()
                                 for i in line.split()),
                                dtype='int64',
                                count=3*self.num_grids).reshape((self.num_grids, 3))
        # we need grid_info in `populate_grid_objects`, so save it to self
        self.pgrid_info = grid_info

    def __cache_endianness(self, test_grid):
        """
        Cache the endianness and bytes per real of the grids by using a test
        grid and assuming that all grids have the same endianness. This is a
        pretty safe assumption since Nyx uses one file per processor (@todo:
        make sure this is still true, I don't think so). If you are running on a
        cluster with different endian processors, then you are disappoint.

        """
        # open the test file & grab the header
        inFile = open(os.path.expanduser(test_grid.filename[self.field_list[0]]), 'rb')
        header = inFile.readline()
        inFile.close()
        header.strip()

        headerRe = re.compile(fab_header_pattern)
        bytesPerReal, endian, start, stop, centerType, nComponents = \
            headerRe.search(header).groups()
        self._bytesPerReal = int(bytesPerReal)
        if self._bytesPerReal == int(endian[0]):
            dtype = '<'
        elif self._bytesPerReal == int(endian[-1]):
            dtype = '>'
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        dtype += ('f%i' % self._bytesPerReal) # always a floating point
        self._dtype = dtype

    def __calculate_grid_dimensions(self, start_stop):
        start = np.array(map(int, start_stop[0].split(',')))
        stop = np.array(map(int, start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension, start, stop

    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")

        self.grids = np.concatenate([level.grids for level in self.levels])
        basedir = self.parameter_file.path
        for g, pg in itertools.izip(self.grids, self.pgrid_info):
            g.particle_filename = os.path.join(basedir, "DM",
                                               "Level_%s" % (g.Level),
                                               "DATA_%04i" % pg[0])
            g.NumberOfParticles = pg[1]
            g._particle_offset = pg[2]

        self.grid_particle_count[:, 0] = self.pgrid_info[:, 1]
        del self.pgrid_info

        gls = np.concatenate([level.ngrids * [level.level] for level in self.levels])
        self.grid_levels[:] = gls.reshape((self.num_grids, 1))
        grid_dcs = np.concatenate([level.ngrids*[self.dx[level.level]]
                                   for level in self.levels], axis=0)

        self.grid_dxs = grid_dcs[:, 0].reshape((self.num_grids, 1))
        self.grid_dys = grid_dcs[:, 1].reshape((self.num_grids, 1))
        self.grid_dzs = grid_dcs[:, 2].reshape((self.num_grids, 1))

        left_edges = []
        right_edges = []
        dims = []
        for level in self.levels:
            left_edges += [g.LeftEdge for g in level.grids]
            right_edges += [g.RightEdge for g in level.grids]
            dims += [g.ActiveDimensions for g in level.grids]

        self.grid_left_edge = np.array(left_edges)
        self.grid_right_edge = np.array(right_edges)
        self.grid_dimensions = np.array(dims)
        self.gridReverseTree = [] * self.num_grids
        self.gridReverseTree = [ [] for i in range(self.num_grids)]  # why the same thing twice?
        self.gridTree = [ [] for i in range(self.num_grids)]

        mylog.debug("Done creating grid objects")

    def _populate_hierarchy(self):
        self.__setup_grid_tree()

        for i, grid in enumerate(self.grids):
            if (i % 1e4) == 0:
                mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)

            grid._prepare_grid()
            grid._setup_dx()

    def __setup_grid_tree(self):
        mask = np.empty(self.grids.size, dtype='int32')
        for i, grid in enumerate(self.grids):
            get_box_grids_level(grid.LeftEdge, grid.RightEdge, grid.Level + 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            children = self.grids[mask.astype("bool")]
            for child in children:
                self.gridReverseTree[child.id].append(i)
                self.gridTree[i].append(weakref.proxy(child))

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        dd["field_indexes"] = self.field_indexes
        GridGeometryHandler._setup_classes(self, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = np.logical_and(mask, (self.grid_levels == (grid.Level + 1)).flat)
        return self.grids[mask]

    def _setup_field_list(self):
        if self.parameter_file.use_particles:
            # We know which particle fields will exist -- pending further
            # changes in the future.
            for field in nyx_particle_field_names:
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
                # Note that we call add_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D,
                # 2D and 3D fields.
                self.pf.field_info.add_field(field, NullFunc,
                                             convert_function=cf,
                                             take_log=False, particle_type=True)

    def _count_grids(self):
        """ this is already provided in ??? """
        pass

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = np.ones((self.num_grids, 3), 'int32')
        self.grid_left_edge = np.zeros((self.num_grids, 3), self.float_type)
        self.grid_right_edge = np.ones((self.num_grids, 3), self.float_type)
        self.grid_levels = np.zeros((self.num_grids, 1), 'int32')
        self.grid_particle_count = np.zeros((self.num_grids, 1), 'int32')

    def _parse_hierarchy(self):
        pass

    def _detect_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []
        for field in self.parameter_file.field_info:
            try:
                fd = self.parameter_file.field_info[field].get_dependencies(
                            pf = self.parameter_file)
            except:
                continue
            available = np.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    def _initialize_state_variables(self):
        """
        Override not to re-initialize num_grids in GridGeometryHandler.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}

class NyxLevel:
    def __init__(self, level, ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []

class NyxStaticOutput(StaticOutput):
    """
    This class is a stripped down class that simply reads and parses *filename*,
    without looking at the Nyx hierarchy.

    """
    _hierarchy_class = NyxHierarchy
    _fieldinfo_fallback = NyxFieldInfo
    _fieldinfo_known = KnownNyxFields

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        pname = args[0].rstrip("/")
        dn = os.path.dirname(pname)
        if len(args) > 1:
            kwargs['paramFilename'] = args[1]

        pfname = kwargs.get("paramFilename", os.path.join(dn, "inputs"))

        # @todo: new Nyx output.
        # We check for the job_info file's existence because this is currently
        # what distinguishes Nyx data from MAESTRO data.
        pfn = os.path.join(pfname)
        if not os.path.exists(pfn): return False
        nyx = any(("nyx." in line for line in open(pfn)))
        maestro = os.path.exists(os.path.join(pname, "job_info"))
        orion = (not nyx) and (not maestro)
        return nyx

    def __init__(self, plotname, param_filename="inputs",
                 fparam_filename="probin", data_style="nyx_native",
                 storage_filename=None):
        """
        Need to override for Nyx file structure, for now.

        The paramfile is usually called "inputs" and there may be a fortran
        inputs file usually called "probin". `plotname` here will be a directory
        name as per BoxLib, data_style will be one of

         * Native
         * IEEE (not implemented in yt)
         * ASCII (not implemented in yt)

        """
        self.storage_filename = storage_filename
        self.parameter_filename = param_filename
        self.fparameter_filename = fparam_filename

        self.path = os.path.abspath(plotname)  # data folder

        # silly inputs and probin file thing (this is on the Nyx todo list)
        self.parameter_file_path = os.path.join(os.path.dirname(self.path),
                                                self.parameter_filename)

        self.fparameter_file_path = os.path.join(os.path.dirname(self.path),
                                                 self.fparameter_filename)

        self.fparameters = {}

        # @todo: quick fix...
        self.use_particles = False

        # @todo: first line
        # runs ``self._parse_parameter_file()``, ``self._set_units()``, and
        # ``self.print_key_parameters()``
        StaticOutput.__init__(self, plotname.rstrip("/"), data_style=data_style)

        # @todo: check all of these and hopefully factor out of the constructor.
        # These should maybe not be hardcoded?
        self.parameters["HydroMethod"] = "nyx"  # always PPM DE
        self.parameters["Time"] = 1.  # default unit is 1...
        self.parameters["DualEnergyFormalism"] = 0  # always off.
        self.parameters["EOSType"] = -1  # default

        # @todo: hopefully delete this after implementing new Nyx output
        if self.fparameters.has_key("mu"):
            self.parameters["mu"] = self.fparameters["mu"]

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various dictionaries.

        """
        self._parse_header_file()

        if os.path.isfile(self.fparameter_file_path):
            self._parse_fparameter_file()

        # Let's read the file
        self.unique_identifier = int(os.stat(self.parameter_filename)[ST_CTIME])
        lines = open(self.parameter_file_path).readlines()

        for line in lines:
            if line.find("#") >= 1:  # Keep the commented lines...
                line = line[:line.find("#")]
            line = line.strip()
            if len(line) < 2 or line.find("#") == 0: # ...but skip comments
                continue

            try:
                param, val_string = map(strip, line.split("="))
            except ValueError:
                mylog.error("ValueError: '%s'", line)

            vals = val_string.split()

            # @todo: don't do this here...
            if nyx_to_enzo_dict.has_key(param):
                param_name = nyx_to_enzo_dict[param]
                vals = map(parameter_type_dict[param_name], vals)
                if len(vals) == 1:
                    self.parameters[param_name] = vals[0]
                else:
                    # don't know why this is special
                    if param_name == "RefineBy":
                        self.parameters[param_name] = vals[0]
                    else:
                        self.parameters[param_name] = vals

            elif param.startswith("geometry.prob_hi"):
                self.domain_right_edge = np.array([float(i) for i in vals])
            elif param.startswith("geometry.prob_lo"):
                self.domain_left_edge = np.array([float(i) for i in vals])
            elif param.startswith("particles.write_in_plotfile"):
                self.use_particles = boxlib_bool_to_int(vals[0])
            elif param.startswith("nyx.lo_bc"):
                self.periodicity = ensure_tuple([i == 0 for i in vals])

        # aliases we need
        self.parameters["TopGridRank"] = len(self.parameters["TopGridDimensions"])
        self.dimensionality = self.parameters["TopGridRank"]
        self.domain_dimensions = self.parameters["TopGridDimensions"]
        self.refine_by = self.parameters.get("RefineBy", 2)  # 2 is silent default? Makes sense I suppose.

        # Nyx is always cosmological.
        self.cosmological_simulation = 1
        self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
        self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
        self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]

        # Read in the `comoving_a` file and parse the value. We should fix this
        # in the new Nyx output format...
        a_file = open(os.path.join(self.path, "comoving_a"))
        a_string = a_file.readline().strip()
        a_file.close()

        # Set the scale factor and redshift
        self.cosmological_scale_factor = float(a_string)
        self.parameters["CosmologyCurrentRedshift"] = 1 / float(a_string) - 1

        # alias
        self.current_redshift = self.parameters["CosmologyCurrentRedshift"]

    def _parse_header_file(self):
        """
        Parses the BoxLib header file to get any parameters stored there.
        Hierarchy information is read out of this file in NyxHierarchy.

        Currently, only Time is read here.

        """
        header_file = open(os.path.join(self.path, "Header"))
        lines = header_file.readlines()  # hopefully this is small
        header_file.close()

        n_fields = int(lines[1])  # this could change
        self.current_time = float(lines[3 + n_fields])  # fragile

    def _parse_fparameter_file(self):
        """
        Parses the fortran parameter file for Nyx. Most of this will be useless,
        but this is where it keeps mu = mass per particle/m_hydrogen. Also it
        contains the cosmological variables.

        """
        # @todo: delete after new Nyx output
        lines = open(self.fparameter_file_path).readlines()
        for line in lines:
            if line.count("=") == 1:
                nyx_param, val_string = map(strip, line.split("="))

                # Check if we care about this param. If so, translate it.
                if nyx_to_enzo_dict.has_key(nyx_param):
                    param = nyx_to_enzo_dict[nyx_param]
                else:
                    continue

                # parse vals string and correct for fortran double syntax
                vals = val_string.split()
                if val_string.count("'") == 0:  # must be floats
                    vals = map(float, [val.replace('D', 'e').replace('d', 'e')
                                       for val in vals])

                # single element or array?
                if len(vals) == 1:
                    self.parameters[param] = vals[0]
                else:
                    self.parameters[param] = vals

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the
        parameter file.

        """
        self.units = {}
        self.time_units = {}

        if len(self.parameters) == 0:  # don't think this is possible, but safe
            self._parse_parameter_file()

        # Masses are always in $ M_{\odot} $
        self.units["particle_mass"] = 1.989e33

        mylog.warning("Length units: setting 1.0 = 1.0 Mpc.")
        self.units.update(mpc_conversion)
        self.units["density"] = self.units["particle_mass"]/(self.units["cm"])**3
        self.units["particle_mass_density"] = self.units["density"]
        self.units["Density"] = 1

        # @todo: enzo-isms
        mylog.warning("Time units: setting 1.0 = Mpc/km s ~ 10^12 yr .")
        self.time_units["s"] = 1.0 / 3.08568025e19
        self.conversion_factors["Time"] = 1.0 / 3.08568025e19

        # velocities are not comoving!
        # Nyx is in km/s so we need to convert to cm/s, hence 1e5
        cf = 1e5 * (self.cosmological_scale_factor)
        for ax in "xyz":
            self.units["particle_velocity_%s" % ax] = cf

        # misc
        self.conversion_factors = defaultdict(lambda: 1.0)  # what is this for? - Steffen: this is to get 1.0 for values not in the dict
        self.time_units["1"] = 1
        self.units["1"] = 1.0
        self.units["unitary"] = 1.0 / (self.domain_right_edge -
                                       self.domain_left_edge).max()

        # time
        for unit in sec_conversion.keys():
            self.time_units[unit] = self.time_units["s"] / sec_conversion[unit]

        # not the most useful right now, but someday
        for key in nyx_particle_field_names:
            self.conversion_factors[key] = 1.0

    def _setup_nounits_units(self):
        z = 0

    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.directory, default)
        return f
