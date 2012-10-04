"""
Data structures for Maestro - borrows heavily from Orion frontend.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Chris Malone <chris.m.malone@gmail.com>
Affiliation: SUNY Stony Brook
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

import re
import os
import weakref
import numpy as np

from collections import \
    defaultdict
from string import \
    strip, \
    rstrip
from stat import \
    ST_CTIME

from yt.funcs import *
from yt.data_objects.grid_patch import \
           AMRGridPatch
from yt.geometry.grid_geometry_handler import \
           GridGeometryHandler
from yt.data_objects.static_output import \
           StaticOutput
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion

from .definitions import \
    maestro2enzoDict, \
    parameterTypes, \
    yt2maestroFieldsDict, \
    maestro_FAB_header_pattern

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from .fields import \
    MaestroFieldInfo, \
    add_maestro_field, \
    KnownMaestroFields


class MaestroGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, LeftEdge, RightEdge, index, level, filename, offset, dimensions,start,stop,paranoia=False,**kwargs):
        AMRGridPatch.__init__(self, index,**kwargs)
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
        h.grid_levels[self.id,0] = self.Level
        h.grid_left_edge[self.id,:] = self.LeftEdge[:]
        h.grid_right_edge[self.id,:] = self.RightEdge[:]
        self.field_indexes = h.field_indexes
        self.Children = h.gridTree[self.id]
        pIDs = h.gridReverseTree[self.id]
        if len(pIDs) > 0:
            self.Parent = [weakref.proxy(h.grids[pID]) for pID in pIDs]
        else:
            self.Parent = None

    def _setup_dx(self):
        # has already been read in and stored in hierarchy
        dx = self.hierarchy.grid_dxs[self.index][0]
        dy = self.hierarchy.grid_dys[self.index][0]
        dz = self.hierarchy.grid_dzs[self.index][0]
        self.dds = np.array([dx, dy, dz])
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "MaestroGrid_%04i" % (self.id)

class MaestroHierarchy(GridGeometryHandler):
    grid = MaestroGrid
    def __init__(self, pf, data_style='maestro'):
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        header_filename = os.path.join(pf.fullplotdir,'Header')
        self.directory = pf.fullpath
        self.data_style = data_style

        # this also sets up the grid objects
        self.readGlobalHeader(header_filename, 
                              self.parameter_file.paranoid_read) 

        pf.current_time = self.Time
        
        self.__cache_endianness(self.levels[-1].grids[-1])
        GridGeometryHandler.__init__(self,pf, self.data_style)
        self._setup_data_io()
        self._setup_field_list()
        self._populate_hierarchy()
        
    def readGlobalHeader(self,filename,paranoid_read):
        """
        read the global header file for an Maestro plotfile output.
        """
        counter = 0
        header_file = open(filename,'r')
        self.__global_header_lines = header_file.readlines()

        # parse the file
        self.maestro_version = self.__global_header_lines[0].rstrip()
        self.n_fields      = int(self.__global_header_lines[1])

        counter = self.n_fields+2
        self.field_list = []
        for i,line in enumerate(self.__global_header_lines[2:counter]):
            self.field_list.append(line.rstrip())

        self.dimension = int(self.__global_header_lines[counter])
        if self.dimension != 3:
            raise RunTimeError("Maestro must be in 3D to use yt.")
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
        self.domainLeftEdge_unnecessary = np.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.domainRightEdge_unnecessary = np.array(map(float,self.__global_header_lines[counter].split()))
        counter += 1
        self.refinementFactor_unnecessary = self.__global_header_lines[counter].split()
        counter += 1
        self.globalIndexSpace_unnecessary = self.__global_header_lines[counter]

        counter += 1 # unused line in Maestro BoxLib
        
        counter += 1
        self.dx = np.zeros((self.n_levels,3))
        for i,line in enumerate(self.__global_header_lines[counter:counter+self.n_levels]):
            self.dx[i] = np.array(map(float,line.split()))

        counter += self.n_levels # unused line in Maestro BoxLib
        
        counter += 1 # unused line in Maestro BoxLib
        
        counter += 1
        # each level is one group with ngrids on it. each grid has 3 lines of 
        # 2 reals
        self.levels = []
        grid_counter = 0
        file_finder_pattern = r"FabOnDisk: (\w+_D_[0-9]{5}) (\d+)\n"
        re_file_finder = re.compile(file_finder_pattern)
        dim_finder_pattern = r"\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \(\d+,\d+,\d+\)\)\n"
        re_dim_finder = re.compile(dim_finder_pattern)
        data_files_pattern = r"Level_[\d]+/"
        data_files_finder = re.compile(data_files_pattern)

        for level in range(0,self.n_levels):
            tmp = self.__global_header_lines[counter].split()
            lev,ngrids,unused = int(tmp[0]),int(tmp[1]),float(tmp[2])

            counter += 1 # unused line in Maestro BoxLib

            counter += 1
            self.levels.append(MaestroLevel(lev,ngrids))
            # open level header, extract file names and offsets for
            # each grid
            # read slightly out of order here: at the end of the lo,hi
            # pairs for x,y,z is a *list* of files types in the Level
            # directory. each type has Header and a number of data
            # files (one per processor)
            tmp_offset = counter + 3*ngrids
            nfiles = 0
            key_off = 0
            files =   {}
            offsets = {}
            while (nfiles+tmp_offset < len(self.__global_header_lines) and 
                   data_files_finder.match(
                    self.__global_header_lines[nfiles+tmp_offset])):
                filen = os.path.join(self.parameter_file.fullplotdir, 
                                     self.__global_header_lines[nfiles+tmp_offset].strip())
                # open each "_H" header file, and get the number of
                # components within it
                level_header_file = open(filen+'_H','r').read()
                start_stop_index = re_dim_finder.findall(level_header_file)
                grid_file_offset = re_file_finder.findall(level_header_file)
                ncomp_this_file = int(level_header_file.split('\n')[2])
                f,o = zip(*grid_file_offset)
                for i in range(ncomp_this_file):
                    key = self.field_list[i+key_off]
                    files[key] = f
                    offsets[key] = o
                    self.field_indexes[key] = i
                key_off += ncomp_this_file
                nfiles += 1
            # convert dict of lists to list of dicts
            fn = []
            off = []
            lead_path = os.path.join(self.parameter_file.fullplotdir,
                                     'Level_%02i'%level)
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
                lo = np.array([xlo,ylo,zlo])
                hi = np.array([xhi,yhi,zhi])
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
        Cache the endianness and bytes per real of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Maestro uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        inFile = open(os.path.expanduser(test_grid.filename[self.field_list[0]]),'rb')
        header = inFile.readline()
        inFile.close()
        header.strip()
        
        # parse it. the pattern is in definitions.py
        headerRe = re.compile(maestro_FAB_header_pattern)
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
        start = np.array(map(int,start_stop[0].split(',')))
        stop = np.array(map(int,start_stop[1].split(',')))
        dimension = stop - start + 1
        return dimension,start,stop
        
    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = np.concatenate([level.grids for level in self.levels])
        self.grid_levels = np.concatenate([level.ngrids*[level.level] for level in self.levels])
        self.grid_levels = self.grid_levels.reshape((self.num_grids,1))
        grid_dcs = np.concatenate([level.ngrids*[self.dx[level.level]] for level in self.levels],axis=0)
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
        self.grid_left_edge = np.array(left_edges)
        self.grid_right_edge = np.array(right_edges)
        self.grid_dimensions = np.array(dims)
        self.gridReverseTree = [] * self.num_grids
        self.gridReverseTree = [ [] for i in range(self.num_grids)]
        self.gridTree = [ [] for i in range(self.num_grids)]
        mylog.debug("Done creating grid objects")

    def _populate_hierarchy(self):
        self.__setup_grid_tree()
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
        GridGeometryHandler._setup_classes(self, dd)
        self.object_types.sort()

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        mask = np.logical_and(mask, (self.grid_levels == (grid.Level+1)).flat)
        return self.grids[mask]

    def _setup_field_list(self):
        self.derived_field_list = []
        for field in self.field_info:
            try:
                fd = self.field_info[field].get_dependencies(pf = self.parameter_file)
            except:
                continue
            available = np.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    def _count_grids(self):
        """this is already provided in 

        """
        pass

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = np.ones((self.num_grids,3), 'int32')
        self.grid_left_edge = np.zeros((self.num_grids,3), self.float_type)
        self.grid_right_edge = np.ones((self.num_grids,3), self.float_type)
        self.grid_levels = np.zeros((self.num_grids,1), 'int32')
        self.grid_particle_count = np.zeros((self.num_grids,1), 'int32')

    def _parse_hierarchy(self):
        pass
    
    def _detect_fields(self):
        pass

    def _setup_derived_fields(self):
        pass

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in GridGeometryHandler.__init__

        """
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None
        self._max_locations = {}
    
class MaestroLevel:
    def __init__(self,level,ngrids):
        self.level = level
        self.ngrids = ngrids
        self.grids = []
    

class MaestroStaticOutput(StaticOutput):
    """
    This class is a stripped down class that simply reads and parses
    *filename*, without looking at the Maestro hierarchy.
    """
    _hierarchy_class = MaestroHierarchy
    _fieldinfo_fallback = MaestroFieldInfo
    _fieldinfo_known = KnownMaestroFields

    def __init__(self, plotname, paramFilename=None, 
                 data_style='maestro', paranoia=False,
                 storage_filename = None):
        """need to override for Maestro file structure.

        plotname here will be a directory name
        
        paramFilename is usually named "job_info" and lives in the plotname
        directory.  This file contains the "probin" namelist, which contains
        the simulation parameters.

        """
        self.storage_filename = storage_filename
        self.paranoid_read = paranoia
        self.__ipfn = paramFilename

        StaticOutput.__init__(self, plotname.rstrip("/"),
                              data_style='maestro')

        # this is the unit of time; NOT the current time
        self.parameters["Time"] = 1 # second

        self._parse_header_file()


    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.fullplotdir, default)
        return f

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # fill our args
        pname = args[0].rstrip("/")
        return os.path.exists(os.path.join(pname, "job_info"))
                
        
    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        _n_cellx = _n_celly = _n_cellz = 0
        _prob_hi_x = _prob_hi_y = _prob_hi_z = 0.0
        _prob_lo_x = _prob_lo_y = _prob_lo_z = 0.0

        local_opts = {"n_cellx": int, 
                      "n_celly": int,
                      "n_cellz": int,
                      "prob_hi_x": float, 
                      "prob_hi_y": float, 
                      "prob_hi_z": float,
                      "prob_lo_x": float, 
                      "prob_lo_y": float, 
                      "prob_lo_z": float
                      }

        self.fullplotdir = os.path.abspath(self.parameter_filename)
        self.parameter_filename = self._localize(
                self.__ipfn, 'job_info')
        # Let's read the file
        self.unique_identifier = \
            int(os.stat(self.fullplotdir)[ST_CTIME])

        _parameters = NameList("probin",self.parameter_filename)
        for param, val in _parameters:
            if local_opts.has_key(param):
                exec("_%s = %s" % (param, val))
            elif maestro2enzoDict.has_key(param):
                paramName = maestro2enzoDict[param]
                t = parameterTypes[paramName](val)
                exec("self.%s = %s" % (paramName,t))

        self.domain_dimensions = np.array([_n_cellx,_n_celly,_n_cellz])
        self.domain_left_edge = np.array([_prob_lo_x,_prob_lo_y,_prob_lo_z])
        self.domain_right_edge = np.array([_prob_hi_x,_prob_hi_y,_prob_hi_z])
        
        self.cosmological_simulation = self.current_redshift = \
            self.omega_matter = self.omega_lambda = self.hubble_constant = 0

        # set the current time to zero.  we will update this when we read the
        # header file
        self.current_time = 0.0


    def _parse_header_file(self):
        """
        Parses the BoxLib header file to get any parameters stored
        there. Hierarchy information is read out of this file in
        MaestroHierarchy. 

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
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]
        for key in yt2maestroFieldsDict:
            self.conversion_factors[key] = 1.0

    def _setup_nounits_units(self):
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]


class NameList(object):
    """A simple class for storing the data in a FORTRAN namelist as a dict.
    """
    
    def __init__(self, name=None, filename=None):
        """Initialize the class.
        
        Arguments:
        - `name`: the name of the namelist as seen in FORTRAN
        - `filename`: the filename which contains namelist name
        """
        self._name = name.lower()
        self._filename = filename

        self._namelist = self._build_namelist()

    def __getitem__(self, key):
        """
        Returns the value of variable key from the namelist.
        """
        if key in self._namelist: return self._namelist[key]
        raise KeyError(key)

    def __iter__(self):
        return self._namelist.iteritems()

    def _build_namelist(self):
        """
        Parse self.filename and grab the information from the namelist
        """
        lines = open(self._filename).read().lower()

        namelist_name_pattern = "&%s" % self._name
        namelist_name_finder = re.compile(namelist_name_pattern)
        namelist_end_pattern = r"/"
        namelist_end_finder = re.compile(namelist_end_pattern)

        start = namelist_name_finder.search(lines).start()
        end = namelist_end_finder.search(lines,start).start()

        namelist_data_pattern = r"\s*(\w+ *= *[\w\.\-%\+]+)"
        namelist_data_finder = re.compile(namelist_data_pattern)
        
        key_value_pairs = namelist_data_finder.findall(lines,start,end)

        namelist = {}

        for key_value_pair in key_value_pairs:
            key, value = map(strip,map(rstrip,key_value_pair.split("=")))
            namelist[key] = value

        return namelist
