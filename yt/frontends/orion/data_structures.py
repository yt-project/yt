"""
Data structures for Orion. 

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

import numpy as na

from yt.funcs import *
from yt.data_objects.field_info_container import FieldInfoContainer, NullFunc
from yt.data_objects.grid_patch import AMRGridPatch
from yt.data_objects.hierarchy import AMRHierarchy
from yt.data_objects.static_output import StaticOutput
from yt.utilities.definitions import mpc_conversion
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only

from .definitions import \
    orion2enzoDict, \
    parameterDict, \
    yt2orionFieldsDict, \
    orion_FAB_header_pattern
from .fields import \
    OrionFieldInfo, \
    add_orion_field, \
    KnownOrionFields


class OrionGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, LeftEdge, RightEdge, index, level, filename, offset,
                 dimensions, start, stop, paranoia=False, **kwargs):
        AMRGridPatch.__init__(self, index, **kwargs)
        self.filename = filename
        self._offset = offset
        self._paranoid = paranoia
        
        # should error check this
        self.ActiveDimensions = (dimensions.copy()).astype('int32')#.transpose()
        self.start_index = start.copy()#.transpose()
        self.stop_index = stop.copy()#.transpose()
        self.LeftEdge  = LeftEdge.copy()
        self.RightEdge = RightEdge.copy()
        self.index = index
        self.Level = level

    def get_global_startindex(self):
        return self.start_index

    def _prepare_grid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # This is definitely the slowest part of generating the hierarchy
        # Now we give it pointers to all of its attributes
        # Note that to keep in line with Enzo, we have broken PEP-8
        h = self.hierarchy # cache it
        #self.StartIndices = h.gridStartIndices[self.id]
        #self.EndIndices = h.gridEndIndices[self.id]
        h.grid_levels[self.id,0] = self.Level
        h.grid_left_edge[self.id,:] = self.LeftEdge[:]
        h.grid_right_edge[self.id,:] = self.RightEdge[:]
        #self.Time = h.gridTimes[self.id,0]
        self.NumberOfParticles = 0 # these will be read in later
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
        return "OrionGrid_%04i" % (self.id)

class OrionHierarchy(AMRHierarchy):
    grid = OrionGrid
    def __init__(self, pf, data_style='orion_native'):
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        header_filename = os.path.join(pf.fullplotdir,'Header')
        self.directory = pf.fullpath
        self.data_style = data_style

        self.readGlobalHeader(header_filename,self.parameter_file.paranoid_read) # also sets up the grid objects
        self.__cache_endianness(self.levels[-1].grids[-1])
        AMRHierarchy.__init__(self,pf, self.data_style)
        self._populate_hierarchy()
        self._read_particles()

    def _read_particles(self):
        """
        reads in particles and assigns them to grids. Will search for
        Star particles, then sink particles if no star particle file
        is found, and finally will simply note that no particles are
        found if neither works. To add a new Orion particle type,
        simply add it to the if/elif/else block.

        """
        self.grid_particle_count = na.zeros(len(self.grids))

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
                mask=na.ones(self.num_grids)
                for i in xrange(len(coord)):
                    na.choose(na.greater(self.grid_left_edge[:,i],coord[i]), (mask,0), mask)
                    na.choose(na.greater(self.grid_right_edge[:,i],coord[i]), (0,mask), mask)
                ind = na.where(mask == 1)
                selected_grids = self.grids[ind]
                # in orion, particles always live on the finest level.
                # so, we want to assign the particle to the finest of
                # the grids we just found
                if len(selected_grids) != 0:
                    grid = sorted(selected_grids, key=lambda grid: grid.Level)[-1]
                    ind = na.where(self.grids == grid)[0][0]
                    self.grid_particle_count[ind] += 1
                    self.grids[ind].NumberOfParticles += 1
        return True
                
    def readGlobalHeader(self,filename,paranoid_read):
        """
        read the global header file for an Orion plotfile output.
        """
        counter = 0
        header_file = open(filename,'r')
        self.__global_header_lines = header_file.readlines()

        # parse the file
        self.orion_version = self.__global_header_lines[0].rstrip()
        self.n_fields      = int(self.__global_header_lines[1])

        counter = self.n_fields+2
        self.field_list = []
        for i,line in enumerate(self.__global_header_lines[2:counter]):
            self.field_list.append(line.rstrip())

        # this is unused...eliminate it?
        #for f in self.field_indexes:
        #    self.field_list.append(orion2ytFieldsDict.get(f,f))

        self.dimension = int(self.__global_header_lines[counter])
        if self.dimension != 3:
            raise RunTimeError("Orion must be in 3D to use yt.")
        counter += 1
        self.Time = float(self.__global_header_lines[counter])
        counter += 1
        self.finest_grid_level = int(self.__global_header_lines[counter])
        self.n_levels = self.finest_grid_level + 1
        counter += 1
        # quantities with _unnecessary are also stored in the inputs
        # file and are not needed.  they are read in and stored in
        # case in the future we want to enable a "backwards" way of
        # taking the data out of the Header file and using it to fill
        # in in the case of a missing inputs file
        self.domainLeftEdge_unnecessary = na.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.domainRightEdge_unnecessary = na.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.refinementFactor_unnecessary = self.__global_header_lines[counter].split() #na.array(map(int,self.__global_header_lines[counter].split()))
        counter += 1
        self.globalIndexSpace_unnecessary = self.__global_header_lines[counter]
        #domain_re.search(self.__global_header_lines[counter]).groups()
        counter += 1
        self.timestepsPerLevel_unnecessary = self.__global_header_lines[counter]
        counter += 1
        self.dx = na.zeros((self.n_levels,3))
        for i,line in enumerate(self.__global_header_lines[counter:counter+self.n_levels]):
            self.dx[i] = na.array(map(float,line.split()))
        counter += self.n_levels
        self.geometry = int(self.__global_header_lines[counter])
        if self.geometry != 0:
            raise RunTimeError("yt only supports cartesian coordinates.")
        counter += 1

        # this is just to debug. eventually it should go away.
        linebreak = int(self.__global_header_lines[counter])
        if linebreak != 0:
            raise RunTimeError("INTERNAL ERROR! This should be a zero.")
        counter += 1

        # each level is one group with ngrids on it. each grid has 3 lines of 2 reals
        self.levels = []
        grid_counter = 0
        file_finder_pattern = r"FabOnDisk: (\w+_D_[0-9]{4}) (\d+)\n"
        re_file_finder = re.compile(file_finder_pattern)
        dim_finder_pattern = r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)\n"
        re_dim_finder = re.compile(dim_finder_pattern)
        data_files_pattern = r"Level_[\d]/"
        data_files_finder = re.compile(data_files_pattern)

        for level in range(0,self.n_levels):
            tmp = self.__global_header_lines[counter].split()
            # should this be grid_time or level_time??
            lev,ngrids,grid_time = int(tmp[0]),int(tmp[1]),float(tmp[2])
            counter += 1
            nsteps = int(self.__global_header_lines[counter])
            counter += 1
            self.levels.append(OrionLevel(lev,ngrids))
            # open level header, extract file names and offsets for
            # each grid
            # read slightly out of order here: at the end of the lo,hi
            # pairs for x,y,z is a *list* of files types in the Level
            # directory. each type has Header and a number of data
            # files (one per processor)
            tmp_offset = counter + 3*ngrids
            nfiles = 0
            key_off = 0
            files =   {} # dict(map(lambda a: (a,[]),self.field_list))
            offsets = {} # dict(map(lambda a: (a,[]),self.field_list))
            while nfiles+tmp_offset < len(self.__global_header_lines) and data_files_finder.match(self.__global_header_lines[nfiles+tmp_offset]):
                filen = os.path.join(self.parameter_file.fullplotdir, \
                                     self.__global_header_lines[nfiles+tmp_offset].strip())
                # open each "_H" header file, and get the number of
                # components within it
                level_header_file = open(filen+'_H','r').read()
                start_stop_index = re_dim_finder.findall(level_header_file) # just take the last one
                grid_file_offset = re_file_finder.findall(level_header_file)
                ncomp_this_file = int(level_header_file.split('\n')[2])
                for i in range(ncomp_this_file):
                    key = self.field_list[i+key_off]
                    f,o = zip(*grid_file_offset)
                    files[key] = f
                    offsets[key] = o
                    self.field_indexes[key] = i
                key_off += ncomp_this_file
                nfiles += 1
            # convert dict of lists to list of dicts
            fn = []
            off = []
            lead_path = os.path.join(self.parameter_file.fullplotdir,'Level_%i'%level)
            for i in range(ngrids):
                fi = [os.path.join(lead_path,files[key][i]) for key in self.field_list]
                of = [int(offsets[key][i]) for key in self.field_list]
                fn.append(dict(zip(self.field_list,fi)))
                off.append(dict(zip(self.field_list,of)))

            for grid in range(0,ngrids):
                gfn = fn[grid]  # filename of file containing this grid
                gfo = off[grid] # offset within that file
                xlo,xhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                ylo,yhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                zlo,zhi = map(float,self.__global_header_lines[counter].split())
                counter+=1
                lo = na.array([xlo,ylo,zlo])
                hi = na.array([xhi,yhi,zhi])
                dims,start,stop = self.__calculate_grid_dimensions(start_stop_index[grid])
                self.levels[-1].grids.append(self.grid(lo,hi,grid_counter,level,gfn, gfo, dims,start,stop,paranoia=paranoid_read,hierarchy=self))
                grid_counter += 1 # this is global, and shouldn't be reset
                                  # for each level

            # already read the filenames above...
            counter+=nfiles
            self.num_grids = grid_counter
            self.float_type = 'float64'

        self.maxLevel = self.n_levels - 1 
        self.max_level = self.n_levels - 1
        header_file.close()

    def __cache_endianness(self,test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Orion uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        inFile = open(os.path.expanduser(test_grid.filename[self.field_list[0]]),'rb')
        header = inFile.readline()
        inFile.close()
        header.strip()
        
        # parse it. the patter is in OrionDefs.py
        headerRe = re.compile(orion_FAB_header_pattern)
        bytesPerReal,endian,start,stop,centerType,nComponents = headerRe.search(header).groups()
        self._bytesPerReal = int(bytesPerReal)
        if self._bytesPerReal == int(endian[0]):
            dtype = '<'
        elif self._bytesPerReal == int(endian[-1]):
            dtype = '>'
        else:
            raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

        dtype += ('f%i' % self._bytesPerReal) # always a floating point
        self._dtype = dtype

    def __calculate_grid_dimensions(self,start_stop):
        start = na.array(map(int,start_stop[0].split(',')))
        stop = na.array(map(int,start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension,start,stop
        
    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = na.concatenate([level.grids for level in self.levels])
        self.grid_levels = na.concatenate([level.ngrids*[level.level] for level in self.levels])
        self.grid_levels = na.array(self.grid_levels.reshape((self.num_grids,1)),dtype='int32')
        grid_dcs = na.concatenate([level.ngrids*[self.dx[level.level]] for level in self.levels],axis=0)
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
        self.__setup_grid_tree()
        #self._setup_grid_corners()
        for i, grid in enumerate(self.grids):
            if (i%1e4) == 0: mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid._prepare_grid()
            grid._setup_dx()

    def __setup_grid_tree(self):
        for i, grid in enumerate(self.grids):
            children = self._get_grid_children(grid)
            for child in children:
                self.gridReverseTree[child.id].append(i)
                self.gridTree[i].append(weakref.proxy(child))

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        dd["field_indexes"] = self.field_indexes
        AMRHierarchy._setup_classes(self, dd)
        #self._add_object_class('grid', "OrionGrid", OrionGridBase, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.grid_levels == (grid.Level+1)).flat)
        return self.grids[mask]

    def _count_grids(self):
        """this is already provided in 

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
    
class OrionLevel:
    def __init__(self,level,ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []
    

class OrionStaticOutput(StaticOutput):
    """
    This class is a stripped down class that simply reads and parses
    *filename*, without looking at the Orion hierarchy.
    """
    _hierarchy_class = OrionHierarchy
    _fieldinfo_fallback = OrionFieldInfo
    _fieldinfo_known = KnownOrionFields

    def __init__(self, plotname, paramFilename=None, fparamFilename=None,
                 data_style='orion_native', paranoia=False,
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

        self.fparameters = {}

        StaticOutput.__init__(self, plotname.rstrip("/"),
                              data_style='orion_native')

        # These should maybe not be hardcoded?
        self.parameters["HydroMethod"] = 'orion' # always PPM DE
        self.parameters["Time"] = 1. # default unit is 1...
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
            self._parse_fparameter_file()
            for param in self.fparameters:
                if orion2enzoDict.has_key(param):
                    self.parameters[orion2enzoDict[param]]=self.fparameters[param]
        # Let's read the file
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[ST_CTIME])
        lines = open(self.parameter_filename).readlines()
        for lineI, line in enumerate(lines):
            if line.find("#") >= 1: # Keep the commented lines...
                line=line[:line.find("#")]
            line=line.strip().rstrip()
            if len(line) < 2 or line.find("#") == 0: # ...but skip comments
                continue
            try:
                param, vals = map(strip,map(rstrip,line.split("=")))
            except ValueError:
                mylog.error("ValueError: '%s'", line)
                continue
            if orion2enzoDict.has_key(param):
                paramName = orion2enzoDict[param]
                t = map(parameterDict[paramName], vals.split())
                if len(t) == 1:
                    self.parameters[paramName] = t[0]
                else:
                    if paramName == "RefineBy":
                        self.parameters[paramName] = t[0]
                    else:
                        self.parameters[paramName] = t
                
            elif param.startswith("geometry.prob_hi"):
                self.domain_right_edge = \
                    na.array([float(i) for i in vals.split()])
            elif param.startswith("geometry.prob_lo"):
                self.domain_left_edge = \
                    na.array([float(i) for i in vals.split()])

        self.parameters["TopGridRank"] = len(self.parameters["TopGridDimensions"])
        self.dimensionality = self.parameters["TopGridRank"]
        self.domain_dimensions = na.array(self.parameters["TopGridDimensions"],dtype='int32')
        self.refine_by = self.parameters["RefineBy"]

        if self.parameters.has_key("ComovingCoordinates") and bool(self.parameters["ComovingCoordinates"]):
            self.cosmological_simulation = 1
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]
            a_file = open(os.path.join(self.fullplotdir,'comoving_a'))
            line = a_file.readline().strip()
            a_file.close()
            self.parameters["CosmologyCurrentRedshift"] = 1/float(line) - 1
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def _parse_fparameter_file(self):
        """
        Parses the fortran parameter file for Orion. Most of this will
        be useless, but this is where it keeps mu = mass per
        particle/m_hydrogen.
        """
        lines = open(self.fparameter_filename).readlines()
        for line in lines:
            if line.count("=") == 1:
                param, vals = map(strip,map(rstrip,line.split("=")))
                if vals.count("'") == 0:
                    t = map(float,[a.replace('D','e').replace('d','e') for a in vals.split()]) # all are floating point.
                else:
                    t = vals.split()
                if len(t) == 1:
                    self.fparameters[param] = t[0]
                else:
                    self.fparameters[param] = t

    def _parse_header_file(self):
        """
        Parses the BoxLib header file to get any parameters stored
        there. Hierarchy information is read out of this file in
        OrionHierarchy. 

        Currently, only Time is read here.
        """
        header_file = open(os.path.join(self.fullplotdir,'Header'))
        lines = header_file.readlines()
        header_file.close()
        n_fields = int(lines[1])
        self.current_time = float(lines[3+n_fields])


                
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9
        for key in yt2orionFieldsDict:
            self.conversion_factors[key] = 1.0

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

