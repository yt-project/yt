"""
ART-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as na
import stat
import weakref
import struct

import os
import os.path


from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from .fields import ARTFieldInfo, KnownARTFields
from .definitions import ART_header
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.amr_utils import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
import yt.utilities.fortran_utils as fpu
from yt.geometry.oct_container import \
    RAMSESOctreeContainer

from yt.frontends.art.definitions import * 
from yt.frontends.art.io import _read_child_mask_level
from yt.frontends.art.io import read_particles
from yt.frontends.art.io import read_stars
from yt.frontends.art.io import _count_art_octs
from yt.frontends.art.io import _read_art_level_info
from yt.frontends.art.io import _read_art_child
from yt.frontends.art.io import _skip_record
from yt.frontends.art.io import _read_record
from yt.frontends.art.io import _read_frecord
from yt.frontends.art.io import _read_record_size
from yt.frontends.art.io import _read_struct
from yt.frontends.art.io import b2t

#Model these classes after the RAMSES data structures, 
#but we don't have domains

class ARTFile(object):
    def __init__(self,file_amr,discover_particles=True,
                 file_amr_template="10MpcBox_csf512_%s.d",
                 file_particle_header_template="PMcrd%s.DAT",
                 file_particle_data_template="PMcrs0%s.DAT",
                 file_star_data_template="stars_%s.dat",
                 file_particle_header=None, 
                 file_particle_data=None,
                 file_star_data=None,
                 limit_level = na.inf):
        """Discover all of the necesary files. Using the supplied
           file_amr the '%s' in file_amr_template is found and, 
           when supplied template file names, particle files
           are found.
        
        Arguments:
        file_amr -- The file holding the hydro octree

        Keyword arguments:
        discover_particles   -- If True, tries to find particle files
        file_amr_template    -- This glob is used to discover other files
        file_*template       -- Particle filename conventions
        file_particle_header -- Particle header file location
        file_particle_data   -- Particle data file location
        file_star_data       -- Star data file location
        limit_level          -- Ignore AMR past this level,
                                good for faster testing

        """ 
        if limit_level is None:
            self.limit_level = na.inf
        else:
            mylog.info("Using maximum level: %i",limit_level)
            self.limit_level = limit_level
        
        amr_prefix  = file_amr_template[:file_amr_template.find('%s')]
        amr_postfix = file_amr_template[file_amr_template.find('%s')+1:]
        file_uid    = file_amr.replace(amr_prefix,'').replace(amr_postfix,'')
        dirn = os.path.dirname(amr_filename)
        base = os.path.basename(amr_filename)
        
        self.file_amr= file_amr
        self.file_particle_header = file_particle_header
        self.file_particle_data = file_particle_data
        self.file_star_data = file_star_data

        if discover_particles:
            if file_particle_header is None:
                loc = dirn+file_particle_header_template%file_uid
                if os.path.exists(loc):
                    self.file_particle_header = loc
                    mylog.info("Discovered particle header: %s",
                            os.path.basename(loc))
            if file_particle_data is None:
                loc = dirn+file_particle_data_template%file_uid
                if os.path.exists(loc):
                    self.file_particle_data = loc
                    mylog.info("Discovered particle data:   %s",
                                os.path.basename(loc))
            if file_star_data is None:
                loc = dirn+file_star_data_template%file_uid
                if os.path.exists(loc):
                    self.file_star_data = loc
                    mylog.info("Discovered stellar data:    %s",
                                os.path.basename(loc))

    def _read_amr_header(self):
        """Read the ART header into the self.parameters dictionary"""
        f = open(self.file_amr,'rb')
        self.parameters = {}
        for format,name in  amr_header_struct.iteritems():
            size = struct.calcsize(format)
            output = struct.unpac(format,f.read(size))[0]
            self.parameters[name] = output
        for k,v in parameters:
            self.parameters[k]=v
        for to_skip in ['tl','dtl','tlold','dtlold','iSO']:
            _skip_record(f)
        
        #file info structural variables
        self.ncell, = struct.unpack('>l', _read_record(f))
        self.root_grid_mask_offset = f.tell()
        self.domain_dimensions  = na.ones(3,dtype='int64')
        self.domain_dimensions *= na.rint(self.ncell**(1./3.))
        
        root_cells = self.domain_dimensions.prod()
        self.root_iOctCh = _read_frecord(f,'>i')[:root_cells]
        self.root_iOctCh = self.root_iOctCh.reshape(self.domain_dimensions,
                                                    order='F')
        self.root_grid_offset = f.tell()
        _skip_record(f) # hvar
        _skip_record(f) # var
        self.iOctFree, self.nOct = struct.unpack('>ii', _read_record(f))
        self.child_grid_offset = f.tell()
        
        f.close()
        self._update_parameters()
        
    def _read_amr(self, oct_handler):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain 
           are needed - it's position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        f = open(self.file_amr, "rb")
        mylog.debug("Reading AMR file %s"%self.file_amr)
        
        (self.nhydro_vars, self.level_info,
        self.level_oct_offsets, 
        self.level_offsets) = \
                         _count_art_octs(f, 
                          self.child_grid_offset,
                          self.min_level, self.max_level)
        self.level_info = na.array(self.level_info)        
        self.level_info[0]=self.ncell
        self.level_offsets = na.array(self.level_offsets, dtype='int64')
        self.level_offsets[0] = self.root_grid_offset

        #iterate over all of the level, skip the root grid
        for level in xrange(1,len(self.level_info)):
            if self.level_info[level]==0:
                continue
            if level > self.limit_level: 
                mylog.debug("Reached the level limit")
                continue

            if level == 0:
                #create the root grid a priori
                nd = self.domain_dimensions[0]
                nocts = nd.prod()
                indices = na.arange(nocts,dtype='int64')
                #put the indices in fortran order
                #matching that of hydrovars read off disk
                indices = na.reshape(indices,(nd**3,3),order='F') 
                indices = na.ravel(indices)
                left_index = na.mgrid[0:nd,0:nd,0:nd] 
                left_index = na.reshape(left_index.T,(nd**3,3),order='F')
            else: 
                #find the LE for every oct
                left_index, file_loc, nocts = \
                        _read_art_level_info(f,self.level_oct_offsets,level)

            #the first arg is the domain ID, which for RAMSES is separated
            #by CPU  boundaries. It's always zero for us.
            cpu_map = na.empty((nocts,8),dtype="int64")
            oct_handler.add(0,level,nocts,left_index,indices,cpu_map)
        f.close()

    def _update_parameters(self):
        """Update the parameters, fill in the expected properties"""
        self.parameters = update_parameters(self.parameters)

        self.current_redshift = self.parameters['aexpn']**-1.0-1.0
        self.current_time = self.parameters['t']
        self.omega_lambda = self.parameters['Oml0']
        self.omega_matter = self.parameters['Om0']
        self.hubble_constant = header_vals['hubble']
        self.min_level = header_vals['min_level']
        self.max_level = header_vals['max_level']
        self.hubble_time  = 1.0/(self.hubble_constant*100/3.08568025e19)
        self.current_time = b2t(self.parameters['t'])        

class ARTGeometryHandler(OctreeGeometryHandler):

    def __init__(self,pf,data_style="art"):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        self.pf= weakref.proxy(pf)
        self.hierarchy_filename = self.parameter_file.file_amr
        self.directory = os.path.dirname(self.hiearchy_filename) 

        #initializes octs, reads in tree
        super(ARTGeometryHandler,self).__init__(pf,data_style)

        #read in particles & headers

    def _initialize_oct_handler(self):
        self.amr = ARTFile(self.pf) 
        
        self.amr._read_amr_header()
        self.amr._update_parameters()
        nocts = self.amr.level_info.sum()
        
        #allocate octs
        mylog.debug("Allocating %s octs", nocts)
        le = na.array([0.,0.,0.]) #left edge
        re = self.amr.domain_dimensions #right edge
        dx = self.amr.domain_dimensions[0] #number of cells on root side
        self.oct_handler = RAMSESOctreeContainer(dx,le,re,nocts)
        
        #this actually reads every oct and loads it into the octree
        self.amr._read_amr(self.oct_handler)

    def _detect_fields(self):
        self.field_list = field_list #imported from definitions
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _count_selection(self, dobj, mask = None, oct_indices = None):
        if mask is None:
            mask = dobj.selector.select_octs(self.oct_handler)
        if oct_indices is None:
            oct_indices = self.oct_handler.count(mask, split = True) 
        domains = getattr(dobj, "_domains", None)
        if domains is None:
            count = [i.size for i in oct_indices]
            nocts = sum(count)
            domains = [dom for dom in self.domains
                       if count[dom.domain_id - 1] > 0]
        count = self.oct_handler.count_cells(dobj.selector, mask)
        return count

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_domains", None) is None:
            mask = dobj.selector.select_octs(self.oct_handler)
            indices = self.oct_handler.count(mask, split = True) 
            count = [i.size for i in indices]
            nocts = sum(count)
            domains = [dom for dom in self.domains
                       if count[dom.domain_id - 1] > 0]
            dobj._grids = na.array(domains, dtype='object')
            dobj.size = self._count_selection(dobj, mask, indices)
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._domains)
        yield YTDataChunk(dobj, "all", oobjs, dobj.size)

    def _chunk_spatial(self, dobj, ngz):
        raise NotImplementedError

    def _chunk_io(self, dobj):
        pass


class ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTGeometryHandler
    _fieldinfo_fallback = ARTFieldInfo
    _fieldinfo_known = KnownARTFields
    
    def __init__(self, filename, data_style='ramses',
                 storage_filename = None):
        # Here we want to initiate a traceback, if the reader is not built.
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
        
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
        seconds = self.parameters['unit_t']
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9
        self.conversion_factors["Density"] = self.parameters['unit_d']
        vel_u = self.parameters['unit_l'] / self.parameters['unit_t']
        self.conversion_factors["x-velocity"] = vel_u
        self.conversion_factors["y-velocity"] = vel_u
        self.conversion_factors["z-velocity"] = vel_u

    def _setup_nounits_units(self):
        for unit in mpc_conversion.keys():
            self.units[unit] = self.parameters['unit_l'] * mpc_conversion[unit] / mpc_conversion["cm"]

    def _parse_parameter_file(self):
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'art'
        self.parameters["Time"] = 1. # default unit is 1...

        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self.current_time = self.parameters['time'] * self.parameters['unit_t']
        self.domain_right_edge = na.ones(3, dtype='float64') \
                                           * rheader['boxlen']
        self.domain_left_edge = na.zeros(3, dtype='float64')
        self.domain_dimensions = na.ones(3, dtype='int32') * 2
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / rheader["aexp"]) - 1.0
        self.omega_lambda = rheader["omega_l"]
        self.omega_matter = rheader["omega_m"]
        self.hubble_constant = rheader["H0"]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        fn = args[0]
        if not (fn.starts_with("10MpcBox") and fn.endswith(".d")):
            return False
        return os.path.exists(fn)




        


