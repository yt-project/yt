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

from .definitions import yt_to_art, ARTIOconstants,\
   fluid_fields, particle_fields, particle_star_fields
from _artio_caller import \
    artio_is_valid, artio_fileset
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion 
from .fields import ARTIOFieldInfo, KnownARTIOFields
from .definitions import codetime_fields

from yt.funcs import *
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

import yt.utilities.fortran_utils as fpu

from io import b2t

class ARTIOChunk(object) :

    def __init__(self, pf, selector, sfc_start, sfc_end ):
        self.pf = pf
        self.selector = selector
        self.sfc_start = sfc_start
        self.sfc_end = sfc_end

    _fcoords = None
    def fcoords(self, dobj):
        if self._icoords is None :
            print "Error: ARTIOChunk.fcoords called before fill"
            raise RuntimeError
        return self._fcoords

    _ires = None
    def ires(self, dobj):
        if self._icords is None :
            print "Error: ARTIOChunk.ires called before fill"
            raise RuntimeError
        return self._ires

    def fwidth(self, dobj):
        if self._ires is None :
            print "Error: ARTIOChunk.fwidth called before fill"
            raise RuntimeError
        return 2.**-self.ires

    _icoords = None
    def icoords(self, dobj):
        if self._fcoords is None or self._ires is None or self.fwidth is None :
            print "Error: ARTIOChunk.icoords called before fill fcoords/level"
            raise RuntimeError
        else : 
            self._icoords = (int) (self._fcoords/self.fwidth())
        return self._icoords 

    def fill(self, fields):
        # populate 
        (self._fcoords,self._ires, data) = \
                self.pf._handle.read_grid_chunk( self.selector,
                    self.sfc_min, self.sfc_max, fields )
        return data

    def fill_particles(self,accessed_species, selector, fields):
        art_fields = []
        for f in fields :
            assert (yt_to_art.has_key(f[1])) #fields must exist in ART
            art_fields.append(yt_to_art[f[1]])

        masked_particles = {}
        assert ( art_fields != None )
        self.domain._handle.particle_var_fill(accessed_species, \
                masked_particles, selector, art_fields )

        #convert time variables from code units
        for fieldtype, fieldname in fields :
            if fieldname in codetime_fields : 
                print 'convert time variables from code units'
                masked_particles[yt_to_art[fieldname]] = \
                    self.interpb2t(masked_particles[yt_to_art[fieldname]])
                print 'convert time variables from code units'

        # dhr - make sure these are shallow copies
        tr = {}
        for fieldtype, fieldname in fields :
            tr[fieldname] = masked_particles[yt_to_art[fieldname]]
        return tr

class ARTIOGeometryHandler(GeometryHandler):

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

    def _setup_geometry(self):
        mylog.debug("Initializing Geometry Handler empty for now.")

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return (self.parameter_file.domain_width /(2**self.max_level))
    
    def convert(self, unit):
        return self.parameter_file.conversion_factors[unit]

    def find_max(self, field, finest_levels = 3):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        if (field, finest_levels) in self._max_locations:
            return self._max_locations[(field, finest_levels)]
        mv, pos = self.find_max_cell_location(field, finest_levels)
        self._max_locations[(field, finest_levels)] = (mv, pos)
        return mv, pos

    def find_max_cell_location(self, field, finest_levels = 3):
        source = self.all_data()
        if finest_levels is not False:
            source.min_level = self.max_level - finest_levels
        mylog.debug("Searching for maximum value of %s", field)
        max_val, maxi, mx, my, mz = \
            source.quantities["MaxLocation"](field)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f",
              max_val, mx, my, mz)
        self.pf.parameters["Max%sValue" % (field)] = max_val
        self.pf.parameters["Max%sPos" % (field)] = "%s" % ((mx,my,mz),)
        return max_val, np.array((mx,my,mz), dtype='float64')

    def _detect_fields(self):
        self.fluid_field_list = fluid_fields
        self.particle_field_list = particle_fields
        self.field_list = self.fluid_field_list + self.particle_field_list
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTIOGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            print "Running selector on base grid"
            list_sfc_ranges = self.pf._handle.root_sfc_ranges(dobj.selector)
            print list_sfc_ranges
            print "creating chunks"
            chunks = [ARTIOChunk(d, mask, c)
                       for d, c in zip(self.domains, masked_cell_count) if c > 0]
            dobj._chunk_info = chunks
#            dobj.size = sum(masked_cell_count)
#            dobj.shape = (dobj.size,)
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
        for chunk in oobjs:
            yield YTDataChunk(dobj, "io", [subset], subset.masked_cell_count)

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

#        for particle_field in particle_fields:
#            self.conversion_factors[particle_field] =  1.0
        for ax in 'xyz':
            self.conversion_factors["particle_velocity_%s"%ax] = self.parameters['unit_v']
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]
        self.conversion_factors['particle_mass'] = self.parameters['unit_m']
        self.conversion_factors['particle_creation_time'] =  self.parameters['unit_t']
        self.conversion_factors['particle_mass_msun'] = self.parameters['unit_m']/constants.Msun

        #for mult_halo_profiler.py:
        self.parameters['TopGridDimensions'] = [128,128,128]
        self.parameters['RefineBy'] = 2
        self.parameters['DomainLeftEdge'] = [0,0,0]
        self.parameters['DomainRightEdge'] =  [128,128,128]
        self.parameters['TopGridRank'] = 3 #number of dimensions
       
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
            self.parameters["CosmologyInitialRedshift"] =  self.parameters["initial_redshift"] #for sfr
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

        # hard coded assumption of 3D periodicity (until added to parameter file)
        self.periodicity = (True,True,True)

    @classmethod
    def _is_valid(self, *args, **kwargs) :
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"): return False
        return artio_is_valid(args[0][:-4])

