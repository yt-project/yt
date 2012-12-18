"""
ARTIO-specific data structures

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
import numpy as np
import stat
import weakref
import cStringIO

from _artio_caller import \
    artio_is_valid, artio_fileset , artio_fileset_grid, \
    artio_grid_routines#, read_header 
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion 
from .fields import ARTIOFieldInfo, KnownARTIOFields
#######

###############################
from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput
############################

############################
from yt.utilities.lib import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
############################

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

import yt.utilities.fortran_utils as fpu

from yt.geometry.oct_container import \
    ARTIOOctreeContainer




class ARTIODomainFile(object):
#    _last_mask = None
#    _last_selector_id = None
#    nvar = 6
#    _handle = None

    def __init__(self, pf, domain_id):
        self.pf = pf
        self.domain_id = domain_id
        self.part_fn = "%s.p%03i" % (pf.parameter_filename[:-4],domain_id)
        self.grid_fn = "%s.g%03i" % (pf.parameter_filename[:-4],domain_id)
        self._fileset_prefix = pf.parameter_filename[:-4]
        self._handle = artio_fileset_grid(self._fileset_prefix) 
        
        self.artiogrid = artio_grid_routines(self._handle)
        self._read_grid_header()
        self._read_particle_header()
        self.level_count

    ###########################################
    def _read_particle_header(self):
        if not os.path.exists(self.part_fn):
            self.local_particle_count = 0
            self.particle_field_offsets = {}
            return

    def _read_grid_header(self):
        sfc_max = self._handle.parameters['grid_file_sfc_index'][1]-1
        self.local_oct_count = self.artiogrid.count_octs() #count_octs(self._handle)
        self.level_count = self.artiogrid.level_count
        print 'snl data_structures local_octs count',self.local_oct_count
    ###########################################
        
    def _read_grid(self, oct_handler):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain 
           are needed - it's position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        #oct_handler may not be attributed yet
        self.artiogrid=artio_grid_routines(self._handle) 
        self.artiogrid.grid_pos_fill(oct_handler)
        
    def select(self, selector):
        if id(selector) == self._last_selector_id:
            return self._last_mask
        self._last_mask = selector.fill_mask(self)
        self._last_selector_id = id(selector)
        return self._last_mask

    def count(self, selector):
        if id(selector) == self._last_selector_id:
            if self._last_mask is None: return 0
            return self._last_mask.sum()
        self.select(selector)
        return self.count(selector)

class ARTIODomainSubset(object):

    def __init__(self, domain, mask, cell_count):
        self.mask = mask
        self.domain = domain
        self.oct_handler = domain.pf.h.oct_handler
        self.cell_count = cell_count
        level_counts = self.oct_handler.count_levels(
            self.domain.pf.max_level, self.domain.domain_id, mask)
        level_counts[1:] = level_counts[:-1]
        level_counts[0] = 0
        self.level_counts = np.add.accumulate(level_counts)
        print 'cumulative masked level counts',self.level_counts
        
        self.artiogrid = self.domain.artiogrid
        #self.artiogrid=artio_grid_routines(self._handle) 
        
    def icoords(self, dobj):
        return self.oct_handler.icoords(self.domain.domain_id, self.mask,
                                        self.cell_count,
                                        self.level_counts.copy())

    def fcoords(self, dobj):
        return self.oct_handler.fcoords(self.domain.domain_id, self.mask,
                                        self.cell_count,
                                        self.level_counts.copy())

    def fwidth(self, dobj):
        # Recall domain_dimensions is the number of cells, not octs
        base_dx = 1.0/self.domain.pf.domain_dimensions
        widths = np.empty((self.cell_count, 3), dtype="float64")
        dds = (2**self.ires(dobj))
        for i in range(3):
            widths[:,i] = base_dx[i] / dds
        return widths

    def ires(self, dobj):
        return self.oct_handler.ires(self.domain.domain_id, self.mask,
                                     self.cell_count,
                                     self.level_counts.copy())

    def fill(self, fields):
        oct_handler = self.oct_handler
        all_fields = self.domain.pf.h.fluid_field_list
        min_level = self.domain.pf.min_level
        fields = [f for ft, f in fields]
        print 'all_fields:', all_fields 
        print 'fields:', fields
        


        tr = {}
        for field in fields: 
            tr[field] = np.zeros(self.cell_count, 'float64')
        temp = {}
        nc = self.domain.level_count.sum()
        print 'ncells total',self.cell_count, 'ncells masked', nc
        for field in fields:
            #temp[field] = np.empty((no,8), dtype="float64") 
            temp[field] = np.empty(nc, dtype="float32") 

        #buffer variables 
        self.artiogrid.grid_var_fill(temp, fields)
        
        #mask unused cells 
        oct_handler = self.oct_handler
        level_offset = 0
        level_offset += oct_handler.fill_mask(
             self.domain.domain_id,
            tr, temp, self.mask, level_offset) #[oct_container.pyx]

        return tr
    
    

    def fill_old(self, content, fields):
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        oct_handler = self.oct_handler
        all_fields = self.domain.pf.h.fluid_field_list
        fields = [f for ft, f in fields]
        tr = {}
        filled = pos = level_offset = 0
        min_level = self.domain.pf.min_level
        print 'hydro offset', self.domain.hydro_offset
        
        for field in fields:
            tr[field] = np.zeros(self.cell_count, 'float64')
        for level, offset in enumerate(self.domain.hydro_offset):
            print 'level', level, 'offset', offset
            if offset == -1: continue
            content.seek(offset)
            nc = self.domain.level_count[level]
            temp = {}
            for field in all_fields:
                temp[field] = np.empty((nc, 8), dtype="float64")
            for i in range(8):
                for field in all_fields:
                    if field not in fields:
                        #print "Skipping %s in %s : %s" % (field, level,
                        #        self.domain.domain_id)
                        fpu.skip(content)
                    else:
                        #print "Reading %s in %s : %s" % (field, level,
                        #        self.domain.domain_id)
                        temp[field][:,i] = fpu.read_vector(content, 'd') # cell 1
                        #print field, temp[field][:,i] 
            level_offset += oct_handler.fill_level(self.domain.domain_id, level,
                                   tr, temp, self.mask, level_offset)
            #print "FILL (%s : %s) %s" % (self.domain.domain_id, level, level_offset)
        #print "DONE (%s) %s of %s" % (self.domain.domain_id, level_offset, self.cell_count)
        raise NotImplementedError #snl
        return tr

class ARTIOGeometryHandler(OctreeGeometryHandler):

    def __init__(self, pf, data_style='artio'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.max_level = pf.max_level

        self.float_type = np.float64
        super(ARTIOGeometryHandler, self).__init__(pf, data_style)

    def _initialize_oct_handler(self):
        #domains are the class object ... ncpu == 1 currently 
        #only one file/"domain" for ART
        self.domains = [ARTIODomainFile(self.parameter_file, i + 1)
                        for i in range(self.parameter_file['ncpu'])]
        #this allocates space for the oct tree
        self.oct_handler = ARTIOOctreeContainer(
            self.parameter_file.domain_dimensions, 
            self.parameter_file.domain_left_edge,
            self.parameter_file.domain_right_edge) 
        mylog.debug("Allocating octs")
        self.oct_handler.allocate_domains(
            [dom.local_oct_count for dom in self.domains])
        for dom in self.domains:
            dom._read_grid(self.oct_handler)

    def _detect_fields(self):
        # snl Add additional fields and translator from artio <-> yt
        self.fluid_field_list = [ "Density", "x-velocity", "y-velocity",
	                        "z-velocity", "Pressure", "Metallicity" ]
        pfl = set([])
        for domain in self.domains:
            pfl.update(set(domain.particle_field_offsets.keys()))
        self.particle_field_list = list(pfl)
        self.field_list = self.fluid_field_list + self.particle_field_list
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTIOGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            mask = dobj.selector.select_octs(self.oct_handler)
            counts = self.oct_handler.count_cells(dobj.selector, mask)
            subsets = [ARTIODomainSubset(d, mask, c)
                       for d, c in zip(self.domains, counts) if c > 0]
            dobj._chunk_info = subsets
            dobj.size = sum(counts)
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, dobj.size)

    def _chunk_spatial(self, dobj, ngz):
        raise NotImplementedError

    def _chunk_io(self, dobj):
        # _current_chunk is made from identify_base_chunk 
        #object = dobj._current_chunk.objs or dobj._current_chunk.${dobj._chunk_info}
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], subset.cell_count)


class ARTIOStaticOutput(StaticOutput):
    _handle = None
    _hierarchy_class = ARTIOGeometryHandler
    _fieldinfo_fallback = ARTIOFieldInfo
    _fieldinfo_known = KnownARTIOFields
    
    def __init__(self, filename, data_style='artio',
                 storage_filename = None):
        if self._handle is not None : return
        self._filename = filename
        self._fileset_prefix = filename[:-4]
        self._handle = artio_fileset(self._fileset_prefix) 
        # Here we want to initiate a traceback, if the reader is not built.
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename 

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0: 
            self._parse_parameter_file() 
        for unit in mpc_conversion.keys():
            self.units[unit] = self.parameters['unit_l'] * mpc_conversion[unit] / mpc_conversion["cm"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = self.parameters['unit_t'] / sec_conversion[unit]
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        self.conversion_factors["Density"] = self.parameters['unit_d']
        vel_u = self.parameters['unit_l'] / self.parameters['unit_t']
        self.conversion_factors["x-velocity"] = vel_u
        self.conversion_factors["y-velocity"] = vel_u
        self.conversion_factors["z-velocity"] = vel_u

    def _parse_parameter_file(self):
        # hard-coded -- not provided by headers 
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'artio'
        self.parameters["Time"] = 1. # default unit is 1...
        self.min_level = 0  # ART has min_level=0. No self._handle.parameters['grid_min_level']

        # read header
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        num_grid = (int)(self._handle.parameters["num_root_cells"][0]**(1./3.)+0.5)
        self.domain_dimensions = np.ones(3,dtype='int32') * num_grid
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype='float64')*num_grid

        self.min_level = 0
        self.max_level = self._handle.parameters["max_refinement_level"][0]

        self.current_time = self._handle.parameters["tl"][0]
  
        # detect cosmology
        if self._handle.parameters["abox"] :
            self.cosmological_simulation = True
            self.omega_lambda = self._handle.parameters["OmegaL"][0]
            self.omega_matter = self._handle.parameters["OmegaM"][0]
            self.hubble_constant = self._handle.parameters["hubble"][0]
            self.current_redshift = 1.0/self._handle.parameters["abox"][0] - 1.

            self.parameters["initial_redshift"] = \
                self._handle.parameters["auni_init"][0]

        self.parameters['unit_l'] = self._handle.parameters["length_unit"][0]
        self.parameters['unit_t'] = self._handle.parameters["time_unit"][0]
        self.parameters['unit_m'] = self._handle.parameters["mass_unit"][0]
        self.parameters['unit_d'] = self.parameters['unit_m']/self.parameters['unit_l']**(3.0)

        # hard coded number of domains in ART = 1 ... that may change for parallelization 
        self.parameters['ncpu'] = 1

    @classmethod
    def _is_valid(self, *args, **kwargs) :
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"): return False
        return artio_is_valid(args[0][:-4])

