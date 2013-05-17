"""
Data structures for Boxlib Codes 

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 J. S. Oishi.  All Rights Reserved.

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

import os
import re
import weakref

from collections import defaultdict
from string import strip, rstrip
from stat import ST_CTIME
import glob

import numpy as np

from yt.funcs import *
from yt.data_objects.field_info_container import FieldInfoContainer, NullFunc
from yt.data_objects.grid_patch import AMRGridPatch
from yt.data_objects.hierarchy import AMRHierarchy
from yt.data_objects.static_output import StaticOutput
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.utilities.lib import \
    get_box_grids_level

from .definitions import \
    orion2enzoDict, \
    parameterDict
from .fields import \
    OrionFieldInfo, \
    add_orion_field, \
    KnownOrionFields
from .io import IOHandlerBoxlib

# This is what we use to find scientific notation that might include d's
# instead of e's.
_scinot_finder = re.compile(r"[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?")
# This is the dimensions in the Cell_H file for each level
_dim_finder = re.compile(r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)$")
# This is the line that prefixes each set of data for a FAB in the FAB file
_header_pattern = re.compile(r"^FAB \(\(\d+, \([0-9 ]+\)\),\((\d+), " +
                             r"\(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) " +
                             "\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n")

class BoxlibGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, grid_id, offset, filename = None,
                 hierarchy = None):
        super(BoxlibGrid, self).__init__(grid_id, filename, hierarchy)
        self._offset = offset
        self._parent_id = []
        self._children_ids = []

    def _prepare_grid(self):
        super(BoxlibGrid, self)._prepare_grid()
        my_ind = self.id - self._id_offset
        self.start_index = self.hierarchy.grid_start_index[my_ind]

    def get_global_startindex(self):
        return self.start_index

    def _setup_dx(self):
        # has already been read in and stored in hierarchy
        my_ind = self.id - self._id_offset
        self.dds = self.hierarchy.level_dds[self.Level,:]
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "BoxlibGrid_%04i" % (self.id)

    @property
    def Parent(self):
        if len(self._parent_id) == 0:
            return None
        return [self.hierarchy.grids[pid - self._id_offset]
                for pid in self._parent_id]

    @property
    def Children(self):
        return [self.hierarchy.grids[cid - self._id_offset]
                for cid in self._children_ids]

class BoxlibHierarchy(AMRHierarchy):
    grid = BoxlibGrid
    def __init__(self, pf, data_style='boxlib_native'):
        self.data_style = data_style
        self.header_filename = os.path.join(pf.fullplotdir, 'Header')

        AMRHierarchy.__init__(self, pf, data_style)
        self._cache_endianness(self.grids[-1])
        #self._read_particles()

    def _parse_hierarchy(self):
        """
        read the global header file for an Boxlib plotfile output.
        """
        self.max_level = self.parameter_file._max_level
        header_file = open(self.header_filename,'r')
        # We can now skip to the point in the file we want to start parsing at.
        header_file.seek(self.parameter_file._header_mesh_start)

        dx = []
        for i in range(self.max_level + 1):
            dx.append([float(v) for v in header_file.next().split()])
        self.level_dds = np.array(dx, dtype="float64")
        if int(header_file.next()) != 0:
            raise RunTimeError("yt only supports cartesian coordinates.")
        if int(header_file.next()) != 0:
            raise RunTimeError("INTERNAL ERROR! This should be a zero.")

        # each level is one group with ngrids on it. each grid has 3 lines of 2 reals
        self.grids = []
        grid_counter = 0
        for level in range(self.max_level + 1):
            vals = header_file.next().split()
            # should this be grid_time or level_time??
            lev, ngrids, grid_time = int(vals[0]),int(vals[1]),float(vals[2])
            assert(lev == level)
            nsteps = int(header_file.next())
            for gi in range(ngrids):
                xlo, xhi = [float(v) for v in header_file.next().split()]
                ylo, yhi = [float(v) for v in header_file.next().split()]
                zlo, zhi = [float(v) for v in header_file.next().split()]
                self.grid_left_edge[grid_counter + gi, :] = [xlo, ylo, zlo]
                self.grid_right_edge[grid_counter + gi, :] = [xhi, yhi, zhi]
            # Now we get to the level header filename, which we open and parse.
            fn = os.path.join(self.parameter_file.fullplotdir,
                              header_file.next().strip())
            level_header_file = open(fn + "_H")
            level_dir = os.path.dirname(fn)
            # We skip the first two lines, although they probably contain
            # useful information I don't know how to decipher.
            level_header_file.next()
            level_header_file.next()
            # Now we get the number of data files
            ncomp_this_file = int(level_header_file.next())
            # Skip the next line, and we should get the number of grids / FABs
            # in this file
            level_header_file.next()
            # To decipher this next line, we expect something like:
            # (8 0
            # where the first is the number of FABs in this level.
            ngrids = int(level_header_file.next().split()[0][1:])
            # Now we can iterate over each and get the indices.
            for gi in range(ngrids):
                # components within it
                start, stop = _dim_finder.match(level_header_file.next()).groups()
                start = np.array(start.split(","), dtype="int64")
                stop = np.array(stop.split(","), dtype="int64")
                dims = stop - start + 1
                self.grid_dimensions[grid_counter + gi,:] = dims
                self.grid_start_index[grid_counter + gi,:] = start
            # Now we read two more lines.  The first of these is a close
            # parenthesis.
            level_header_file.next()
            # This line I'm not 100% sure of, but it's either number of grids
            # or number of FABfiles.
            level_header_file.next()
            # Now we iterate over grids to find their offsets in each file.
            for gi in range(ngrids):
                # Now we get the data file, at which point we're ready to
                # create the grid.
                dummy, filename, offset = level_header_file.next().split()
                filename = os.path.join(level_dir, filename)
                go = self.grid(grid_counter + gi, int(offset), filename, self)
                go.Level = self.grid_levels[grid_counter + gi,:] = level
                self.grids.append(go)
            grid_counter += ngrids
            # already read the filenames above...
            self.float_type = 'float64'

    def _cache_endianness(self,test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Boxlib uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        with open(os.path.expanduser(test_grid.filename), 'rb') as f:
            header = f.readline()
        
        bpr, endian, start, stop, centering, nc = \
            _header_pattern.search(header).groups()
        # Note that previously we were using a different value for BPR than we
        # use now.  Here is an example set of information directly from BoxLib:
        #  * DOUBLE data
        #  * FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))((0,0) (63,63) (0,0)) 27
        #  * FLOAT data
        #  * FAB ((8, (32 8 23 0 1 9 0 127)),(4, (1 2 3 4)))((0,0) (63,63) (0,0)) 27
        if bpr == endian[0]:
            dtype = '<f%s' % bpr
        elif bpr == endian[-1]:
            dtype = '>f%s' % bpr
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        mylog.debug("FAB header suggests dtype of %s", dtype)
        self._dtype = np.dtype(dtype)

    def __calculate_grid_dimensions(self, start_stop):
        start = np.array(map(int,start_stop[0].split(',')))
        stop = np.array(map(int,start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension,start,stop
        
    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = np.array(self.grids, dtype='object')
        self._reconstruct_parent_child()
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid._prepare_grid()
            grid._setup_dx()
        mylog.debug("Done creating grid objects")

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

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = np.logical_and(mask, (self.grid_levels == (grid.Level+1)).flat)
        return np.where(mask)

    def _count_grids(self):
        # We can get everything from the Header file, but note that we're
        # duplicating some work done elsewhere.  In a future where we don't
        # pre-allocate grid arrays, this becomes unnecessary.
        header_file = open(self.header_filename, 'r')
        header_file.seek(self.parameter_file._header_mesh_start)
        # Skip over the level dxs, geometry and the zero:
        [header_file.next() for i in range(self.parameter_file._max_level + 3)]
        # Now we need to be very careful, as we've seeked, and now we iterate.
        # Does this work?  We are going to count the number of places that we
        # have a three-item line.  The three items would be level, number of
        # grids, and then grid time.
        self.num_grids = 0
        for line in header_file:
            if len(line.split()) != 3: continue
            self.num_grids += int(line.split()[1])
        
    def _initialize_grid_arrays(self):
        super(BoxlibHierarchy, self)._initialize_grid_arrays()
        self.grid_start_index = np.zeros((self.num_grids,3), 'int64')

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in AMRHierarchy.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}

    def _detect_fields(self):
        # This is all done in _parse_header_file
        self.field_list = self.parameter_file._field_list[:]
        self.field_indexes = dict((f, i)
                                for i, f in enumerate(self.field_list))

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)

class BoxlibStaticOutput(StaticOutput):
    """
    This class is a stripped down class that simply reads and parses
    *filename*, without looking at the Boxlib hierarchy.
    """
    _hierarchy_class = BoxlibHierarchy
    _fieldinfo_fallback = OrionFieldInfo
    _fieldinfo_known = KnownOrionFields

    def __init__(self, plotname, paramFilename=None, fparamFilename=None,
                 data_style='boxlib_native', paranoia=False,
                 storage_filename = None):
        """
        The paramfile is usually called "inputs"
        and there may be a fortran inputs file usually called "probin"
        plotname here will be a directory name
        as per BoxLib, data_style will be Native (implemented here), IEEE (not
        yet implemented) or ASCII (not yet implemented.)
        """
        self.storage_filename = storage_filename
        self.paranoid_read = paranoia
        self.parameter_filename = paramFilename
        self.fparameter_filename = fparamFilename
        self.__ipfn = paramFilename

        StaticOutput.__init__(self, plotname.rstrip("/"),
                              data_style='boxlib_native')

        # These should maybe not be hardcoded?
        self.parameters["HydroMethod"] = 'boxlib'
        self.parameters["Time"] = 1. # default unit is 1...
        self.parameters["EOSType"] = -1 # default

    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.directory, default)
        return f

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        pname = args[0].rstrip("/")
        dn = os.path.dirname(pname)
        if len(args) > 1: kwargs['paramFilename'] = args[1]
        pfname = kwargs.get("paramFilename", os.path.join(dn, "inputs"))

        # We check for the job_info file's existence because this is currently
        # what distinguishes Orion data from MAESTRO data.
        pfn = os.path.join(pfname)
        if not os.path.exists(pfn): return False
        castro = any(("castro." in line for line in open(pfn)))
        nyx = any(("nyx." in line for line in open(pfn)))
        maestro = os.path.exists(os.path.join(pname, "job_info"))
        really_orion = any(("geometry.prob_lo" in line for line in open(pfn)))
        orion = (not castro) and (not maestro) and (not nyx) and really_orion
        return orion
        
    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        self.fullplotdir = os.path.abspath(self.parameter_filename)
        self._parse_header_file()
        self.parameter_filename = self._localize(
                self.__ipfn, 'inputs')
        self.fparameter_filename = self._localize(
                self.fparameter_filename, 'probin')
        if os.path.isfile(self.fparameter_filename):
            self._parse_probin()
        # Let's read the file
        hfn = os.path.join(self.fullplotdir, 'Header')
        self.unique_identifier = int(os.stat(hfn)[ST_CTIME])
        # the 'inputs' file is now optional
        if not os.path.exists(self.parameter_filename):
            return
        for line in (line.split("#")[0].strip() for line in
                     open(self.parameter_filename)):
            if len(line) == 0: continue
            param, vals = [s.strip() for s in line.split("=")]
            if param == "amr.n_cell":
                vals = self.domain_dimensions = np.array(vals.split(), dtype='int32')
            elif param == "amr.ref_ratio":
                vals = self.refine_by = int(vals[0])
            elif param == "Prob.lo_bc":
                vals = self.periodicity = ensure_tuple([i == 0 for i in vals.split()])
            elif param == "castro.use_comoving":
                vals = self.cosmological_simulation = bool(vals)
            else:
                # Now we guess some things about the parameter and its type
                v = vals.split()[0] # Just in case there are multiple; we'll go
                                    # back afterward to using vals.
                try:
                    float(v.upper().replace("D","E"))
                except:
                    pcast = str
                else:
                    syms = (".", "D+", "D-", "E+", "E-")
                    if any(sym in v for sym in syms for v in vals.split()):
                        pcast = float
                    else:
                        pcast = int
                vals = [pcast(v) for v in vals.split()]
                if len(vals) == 1: vals = vals[0]
            self.parameters[param] = vals

        if self.parameters.has_key("ComovingCoordinates") and bool(self.parameters["ComovingCoordinates"]):
            self.cosmological_simulation = 1
            self.omega_lambda = self.parameters["comoving_OmL"]
            self.omega_matter = self.parameters["comoving_OmM"]
            self.hubble_constant = self.parameters["comoving_h"]
            a_file = open(os.path.join(self.fullplotdir,'comoving_a'))
            line = a_file.readline().strip()
            a_file.close()
            self.current_redshift = 1/float(line) - 1
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _parse_probin(self):
        """
        Parses the fortran parameter file for Orion. Most of this will
        be useless, but this is where it keeps mu = mass per
        particle/m_hydrogen.
        """
        for line in (l for l in open(self.fparameter_filename) if "=" in l):
            param, vals = [v.strip() for v in line.split("=")]
            # Now, there are a couple different types of parameters.
            # Some will be where you only have floating point values, others
            # will be where things are specified as string literals.
            # Unfortunately, we're also using Fortran values, which will have
            # things like 1.d-2 which is pathologically difficult to parse if
            # your C library doesn't include 'd' in its locale for strtod.
            # So we'll try to determine this.
            vals = vals.split()
            if any(_scinot_finder.match(v) for v in vals):
                vals = [float(v.replace("D","e").replace("d","e"))
                        for v in vals]
            if len(vals) == 1:
                vals = vals[0]
            self.parameters[param] = vals

    def _parse_header_file(self):
        """
        We parse the Boxlib header, which we use as our basis.  Anything in the
        inputs file will override this, but the inputs file is not strictly
        necessary for orientation of the data in space.
        """

        # Note: Python uses a read-ahead buffer, so using .next(), which would
        # be my preferred solution, won't work here.  We have to explicitly
        # call readline() if we want to end up with an offset at the very end.
        # Fortunately, elsewhere we don't care about the offset, so we're fine
        # everywhere else using iteration exclusively.
        header_file = open(os.path.join(self.fullplotdir,'Header'))
        self.orion_version = header_file.readline().rstrip()
        n_fields = int(header_file.readline())

        self._field_list = [header_file.readline().strip()
                           for i in range(n_fields)]

        self.dimensionality = int(header_file.readline())
        if self.dimensionality != 3:
            raise RunTimeError("Boxlib 1D and 2D support not currently available.")
        self.current_time = float(header_file.readline())
        # This is traditionally a hierarchy attribute, so we will set it, but
        # in a slightly hidden variable.
        self._max_level = int(header_file.readline()) 
        self.domain_left_edge = np.array(header_file.readline().split(),
                                         dtype="float64")
        self.domain_right_edge = np.array(header_file.readline().split(),
                                         dtype="float64")
        ref_factors = np.array([int(i) for i in
                                header_file.readline().split()])
        if ref_factors.size == 0:
            # We use a default of two, as Nyx doesn't always output this value
            ref_factors = [2]
        # We can't vary refinement factors based on dimension, or whatever else
        # they are vaied on.  In one curious thing, I found that some Castro 3D
        # data has only two refinement factors, which I don't know how to
        # understand.
        assert(np.unique(ref_factors).size == 1)
        self.refine_by = ref_factors[0]
        # Now we read the global index space, to get 
        index_space = header_file.readline()
        # This will be of the form:
        #  ((0,0,0) (255,255,255) (0,0,0)) ((0,0,0) (511,511,511) (0,0,0))
        # So note that if we split it all up based on spaces, we should be
        # fine, as long as we take the first two entries, which correspond to
        # the root level.  I'm not 100% pleased with this solution.
        root_space = index_space.replace("(","").replace(")","").split()[:2]
        start = np.array(root_space[0].split(","), dtype="int64")
        stop = np.array(root_space[1].split(","), dtype="int64")
        self.domain_dimensions = stop - start + 1
        # Skip timesteps per level
        header_file.readline()
        self._header_mesh_start = header_file.tell()

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

    def _setup_nounits_units(self):
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]
            
    @parallel_root_only
    def print_key_parameters(self):
        for a in ["current_time", "domain_dimensions", "domain_left_edge",
                  "domain_right_edge"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)

class OrionHierarchy(BoxlibHierarchy):
    def _read_particles(self):
        """
        reads in particles and assigns them to grids. Will search for
        Star particles, then sink particles if no star particle file
        is found, and finally will simply note that no particles are
        found if neither works. To add a new Orion particle type,
        simply add it to the if/elif/else block.

        """
        self.grid_particle_count = np.zeros(len(self.grids))

        for particle_filename in ["StarParticles", "SinkParticles"]:
            fn = os.path.join(self.pf.fullplotdir, particle_filename)
            if os.path.exists(fn): self._read_particle_file(fn)

    def _read_particle_file(self, fn):
        """actually reads the orion particle data file itself.

        """
        if not os.path.exists(fn): return
        with open(fn, 'r') as f:
            lines = f.readlines()
            self.num_stars = int(lines[0].strip()[0])
            for line in lines[1:]:
                particle_position_x = float(line.split(' ')[1])
                particle_position_y = float(line.split(' ')[2])
                particle_position_z = float(line.split(' ')[3])
                coord = [particle_position_x, particle_position_y, particle_position_z]
                # for each particle, determine which grids contain it
                # copied from object_finding_mixin.py
                mask=np.ones(self.num_grids)
                for i in xrange(len(coord)):
                    np.choose(np.greater(self.grid_left_edge[:,i],coord[i]), (mask,0), mask)
                    np.choose(np.greater(self.grid_right_edge[:,i],coord[i]), (0,mask), mask)
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
        return True
                
