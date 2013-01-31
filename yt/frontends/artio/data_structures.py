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
    artio_is_valid, artio_fileset 
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

    def __init__(self, pf, domain_id):
        self.pf = pf
        self.domain_id = domain_id
        self._fileset_prefix = pf.parameter_filename[:-4]
        self.grid_fn = "%s.g%03i" % (pf.parameter_filename[:-4],domain_id)
	self.part_fn = "%s.p%03i" % (pf.parameter_filename[:-4],domain_id)
        #self._handle = artio_fileset(self._fileset_prefix)
        self._handle = self.pf._handle
        self.local_oct_count = self._handle.count_refined_octs() 
        print "Local oct count = ", self.local_oct_count

    #    if self._handle.parameters.has_key('num_particle_files') :
       # self.artioparticle = artio_particle_routines(self._handle)
 #           self._read_particle_header()
           # print "Particle support not implemented"
        self.local_particle_count = 0
        self.particle_field_offsets = {}                                                      
        
    def _read_grid(self, oct_handler):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain 
           are needed - it's position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        self._handle.grid_pos_fill(oct_handler)
        
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

    def __init__(self, domain, mask, masked_cell_count):
        print 'initing domain subset in data_structures.py'
        self.mask = mask
        self.domain = domain
        self.oct_handler = domain.pf.h.oct_handler
        self.masked_cell_count = masked_cell_count
        print 'counting levels in data_structures.py'
        ncum_masked_level = self.oct_handler.count_levels(
            self.domain.pf.max_level, self.domain.domain_id, mask)
        print 'compiling level mask in data_structures.py'
        ncum_masked_level[1:] = ncum_masked_level[:-1]
        ncum_masked_level[0] = 0
        self.ncum_masked_level = np.add.accumulate(ncum_masked_level)
        print 'cumulative masked level counts',self.ncum_masked_level
        
    def icoords(self, dobj):
        return self.oct_handler.icoords(self.domain.domain_id, self.mask,
                                        self.masked_cell_count,
                                        self.ncum_masked_level.copy())

    def fcoords(self, dobj):
        return self.oct_handler.fcoords(self.domain.domain_id, self.mask,
                                        self.masked_cell_count,
                                        self.ncum_masked_level.copy())

    def fwidth(self, dobj):
        # Recall domain_dimensions is the number of cells, not octs
        # snl FIX: please don't hardcode this here 
#        DRE = self.oct_handler.parameter_file.domain_right_edge 
#        DLE = self.oct_handler.parameter_file.domain_left_edge
#        nn = self.oct_handler.parameter_file.domain_dimension
#        for i in range(3):
#            base_dx = (DRE[i] - DLE[i])/nn[i]
#        print 'in fwidth in data_structures.py', DRE, DLE, nn
#        print base_dx
        base_dx = [1.0,1.0,1.0]
        widths = np.empty((self.masked_cell_count, 3), dtype="float64")
        dds = (2**self.ires(dobj))
        for i in range(3):
            widths[:,i] = base_dx[i] / dds
        return widths

    def ires(self, dobj):
        return self.oct_handler.ires(self.domain.domain_id, self.mask,
                                     self.masked_cell_count,
                                     self.ncum_masked_level.copy())

    def fill(self, fields):
        # translate fields into ARTIO names (this dict should be moved to fields.py)
        yt_to_art = { 
            'Density':'HVAR_GAS_DENSITY',
            'TotalEnergy':'HVAR_GAS_ENERGY',
            'GasEnergy':'HVAR_INTERNAL_ENERGY',
            'Pressure':'HVAR_PRESSURE',
            'XMomentumDensity':'HVAR_MOMENTUM_X',
            'YMomentumDensity':'HVAR_MOMENTUM_Y',
            'ZMomentumDensity':'HVAR_MOMENTUM_Z',
            'Gamma':'HVAR_GAMMA',
            'MetalDensitySNIa':'HVAR_METAL_DENSITY_Ia',
            'MetalDensitySNII':'HVAR_METAL_DENSITY_II',
            'Potential':'VAR_POTENTIAL',
            'PotentialHydro':'VAR_POTENTIAL_HYDRO' }

        tr = {}
        for fieldtype, fieldname in fields: 
            tr[fieldname] = np.zeros(self.masked_cell_count, 'float64')

        temp = {}
        for fieldtype, fieldname in fields:
            temp[yt_to_art[fieldname]] = np.empty(8*self.domain.local_oct_count, dtype="float32")  

        #buffer variables 
        self.domain._handle.grid_var_fill(temp, [yt_to_art[f[1]] for f in fields])
        
        # dhr - make sure these are shallow copies 
        temp2 = {}
        for fieldtype, fieldname in fields :
            temp2[fieldname] = temp[yt_to_art[fieldname]]
 
        #mask unused cells (all at once, not level-by-level)
        self.oct_handler.fill_mask( self.domain.domain_id,
                tr, temp2, self.mask, 0 ) 
        
        return tr

    def fill_particles(self,accessed_species, selector, fields):

        yt_to_art = {
            'particle_position_x': 'POSITION_X',
            'particle_position_y': 'POSITION_Y',
            'particle_position_z': 'POSITION_Z',
            'particle_velocity_x': 'VELOCITY_X',
            'particle_velocity_y': 'VELOCITY_Y',
            'particle_velocity_z': 'VELOCITY_Z' }

	yt_to_art_static = {
	'particle_mass': 'particle_species_mass'
	}

        masked_particles = {}
	self.domain._handle.particle_var_fill(accessed_species, masked_particles, selector, [yt_to_art[f[1]] for f in fields], )

	# dhr - make sure these are shallow copies
        tr = {}
        for fieldtype, fieldname in fields :
            tr[fieldname] = temp[yt_to_art[fieldname]]
	
        return tr

 #   def fill_particles(self, fields, accessed_species, masked_particle_count, particle_mask):
  #      fieldnames = [f for ft, f in fields]
  #      print 'all_fields:', all_fields 
  #      print 'fieldnames:', fieldnames
  #      tr = {}
     #   for onefieldname in fieldnames: 
     #       tr[onefieldname] = np.zeros(masked_particle_count, 'float64')
        # temp = {}
        # nc = self.domain.art_level_count.sum()
        # print 'ncells total',self.masked_cell_count, 'ncells masked', nc
        # for onefieldname in fieldnames:
        #    #temp[onefieldname] = np.empty((no,8), dtype="float64") 
        #    temp[onefieldname] = np.empty(nc, dtype="float32") 
        #     
        # #buffer variables 

     #   self.domain.particle_var_fill(tr, fieldnames, accessed_species, particle_mask)

        #     
        #mask unused cells 
        # level_offset += oct_handler.fill_mask(
        #      self.domain.domain_id,
        #     tr, temp, self.mask, level_offset) #[oct_container.pyx] RISM level_offset

     #   return tr

class ARTIOGeometryHandler(OctreeGeometryHandler):

    def __init__(self, pf, data_style='artio'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)

        self.max_level = pf.max_level
        print "max level: ", self.max_level
        self.float_type = np.float64
        super(ARTIOGeometryHandler, self).__init__(pf, data_style)

    def _initialize_oct_handler(self):
        #domains are the class object ... ncpu == 1 currently 
        #only one file/"domain" for ART
        print "Initializing oct container"
        self.domains = [ARTIODomainFile(self.parameter_file, i + 1)
                        for i in range(self.parameter_file['ncpu'])]
        # this allocates space for the oct tree note that 
        # nn is number of root-level OCTS. These don't exist in memory. 
        print 'domain_left_edge, domain_right_edge', self.parameter_file.domain_left_edge, self.parameter_file.domain_right_edge
            
        self.oct_handler = ARTIOOctreeContainer(
            self.parameter_file.domain_dimensions/2, 
            self.parameter_file.domain_left_edge,
            self.parameter_file.domain_right_edge) 
        mylog.debug("Allocating octs")
        self.oct_handler.allocate_domains(
            [dom.local_oct_count for dom in self.domains])
        for dom in self.domains:
            dom._read_grid(self.oct_handler)

    def _detect_fields(self):
        # snl Add additional fields and translator from artio <-> yt
        self.fluid_field_list = [ 'Density', 'TotalEnergy',
                                  'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
                                  'Pressure','Gamma','GasEnergy',
                                  'MetalDensitySNII', 'MetalDensitySNIa',
                                  'Potential','PotentialHydro']

	self.particle_field_list =['particle_position_x', 'particle_position_y',
				   'particle_position_z', 'particle_velocity_x',
				   'particle_velocity_y', 'particle_velocity_z',
				   'particle_mass']
      #  pfl = set([])
     #   for domain in self.domains:
     #       pfl.update(set(domain.particle_field_offsets.keys()))
     #   self.particle_field_list = list(pfl)
        self.field_list = self.fluid_field_list + self.particle_field_list
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTIOGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            mask = dobj.selector.select_octs(self.oct_handler)
            print 'cell count masked in called from data_structures.py'
            masked_cell_count = self.oct_handler.count_cells(dobj.selector, mask)
            print 'done cell count masked in called from data_structures.py'
            print 'calling ARTIODomainSubset from data_structures.py'
            subsets = [ARTIODomainSubset(d, mask, c)
                       for d, c in zip(self.domains, masked_cell_count) if c > 0]
            print 'done with domain subset'
            dobj._chunk_info = subsets
            dobj.size = sum(masked_cell_count)
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]
        print 'done with base chunk'

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
            yield YTDataChunk(dobj, "io", [subset], subset.masked_cell_count)

class ARTIOconstants():
    def __init__(self) : 
        self.yr = 365.25*86400
        self.Myr = 1.0e6*self.yr
        self.Gyr = 1.0e9*self.yr

        self.pc = 3.0856775813e18
        self.kpc = 1.0e3*self.pc
        self.Mpc = 1.0e6*self.pc
        
        self.kms = 1.0e5
        
        self.mp = 1.672621637e-24
        self.k = 1.3806504e-16
        self.G = 6.67428e-8
        self.c = 2.99792458e10
        
        self.eV = 1.602176487e-12
        self.amu = 1.660538782e-24
        self.mH  = 1.007825*self.amu
        self.mHe = 4.002602*self.amu

        self.Msun = 1.32712440018e26/self.G
        self.Zsun = 0.0199
        
        self.Yp    = 0.24
        self.wmu   = 4.0/(8.0-5.0*self.Yp)
        self.wmu_e = 1.0/(1.0-0.5*self.Yp)
        self.XH    = 1.0 - self.Yp
        self.XHe   = 0.25*self.Yp
        self.gamma = 5.0/3.0
        
        self.sigmaT = 6.6524e-25
    
        
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
            
        constants = ARTIOconstants()
        mb = constants.XH*constants.mH + constants.XHe*constants.mHe;
        
        self.parameters['unit_d'] = self.parameters['unit_m']/self.parameters['unit_l']**3.0
        self.parameters['unit_v'] = self.parameters['unit_l']/self.parameters['unit_t']
        self.parameters['unit_E'] = self.parameters['unit_m'] * self.parameters['unit_v']**2.0
        self.parameters['unit_T'] = self.parameters['unit_v']**2.0*mb/constants.k
        self.parameters['unit_rhoE'] = self.parameters['unit_E']/self.parameters['unit_l']**3.0
        self.parameters['unit_nden'] = self.parameters['unit_d']/mb
        self.parameters['Gamma'] = constants.gamma
        
        #         if self.cosmological_simulation :
        #             units_internal.length_in_chimps = unit_factors.length*cosmology->h/constants.Mpc
       
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        self.conversion_factors["Density"] = self.parameters['unit_d']
        self.conversion_factors["x-velocity"] = self.parameters['unit_v']
        self.conversion_factors["y-velocity"] = self.parameters['unit_v']
        self.conversion_factors["z-velocity"] = self.parameters['unit_v']
        self.conversion_factors["Temperature"] = self.parameters['unit_T']*constants.wmu*(constants.gamma-1) #*cell_gas_internal_energy(cell)/cell_gas_density(cell);
        print 'note temperature conversion is currently using fixed gamma not variable'
       
    def _parse_parameter_file(self):
        # hard-coded -- not provided by headers 
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'artio'
        self.parameters["Time"] = 1. # default unit is 1...

        # read header
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # dhr - replace floating point math
        num_grid = self._handle.num_grid
        self.domain_dimensions = np.ones(3,dtype='int32') * num_grid
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype='float64')*num_grid

        self.min_level = 0  # ART has min_level=0. non-existent self._handle.parameters['grid_min_level']
        self.max_level = self._handle.parameters["max_refinement_level"][0]

        self.current_time = self._handle.parameters["tl"][0]
  
        # detect cosmology
        if self._handle.parameters.has_key("abox") :
            abox = self._handle.parameters["abox"][0] 
            self.cosmological_simulation = True
            self.omega_lambda = self._handle.parameters["OmegaL"][0]
            self.omega_matter = self._handle.parameters["OmegaM"][0]
            self.hubble_constant = self._handle.parameters["hubble"][0]
            self.current_redshift = 1.0/self._handle.parameters["abox"][0] - 1.0

            self.parameters["initial_redshift"] = 1.0/self._handle.parameters["auni_init"][0] - 1.0
        else :
            self.cosmological_simulation = False
 
        #units
        if self.cosmological_simulation : 
            self.parameters['unit_m'] = self._handle.parameters["mass_unit"][0]
            self.parameters['unit_t'] = self._handle.parameters["time_unit"][0]*abox**2
            self.parameters['unit_l'] = self._handle.parameters["length_unit"][0]*abox
        else :
            self.parameters['unit_l'] = self._handle.parameters["length_unit"][0]
            self.parameters['unit_t'] = self._handle.parameters["time_unit"][0]
            self.parameters['unit_m'] = self._handle.parameters["mass_unit"][0]

        # hard coded number of domains in ART = 1 ... that may change for parallelization 
        self.parameters['ncpu'] = 1

    @classmethod
    def _is_valid(self, *args, **kwargs) :
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"): return False
        return artio_is_valid(args[0][:-4])

