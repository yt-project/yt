"""
ART-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Christopher Moody <cemoody@ucsc.edu>
Affiliation: UCSC
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
.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import os.path
import glob
import stat
import weakref
import cStringIO

from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.data_objects.hierarchy import \
      AMRHierarchy
from yt.data_objects.static_output import \
      StaticOutput
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from .fields import \
    ARTFieldInfo, add_art_field, KnownARTFields
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.lib import \
    get_box_grids_level
import yt.utilities.lib as amr_utils

from .definitions import *
from io import _read_child_mask_level
from io import _count_art_octs
from io import _read_art_level_info
from io import _skip_record
from io import _read_record
from io import _read_frecord
from io import _read_struct
from io import b2t

import yt.frontends.ramses._ramses_reader as _ramses_reader

from .fields import ARTFieldInfo, KnownARTFields
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.lib import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from yt.utilities.physical_constants import \
    mass_hydrogen_cgs, sec_per_Gyr

class ARTDomainFile(object):
    _last_mask = None
    _last_selector_id = None

    def __init__(self,pf,domain_id,nvar):
        self.nvar = nvar
        self.pf = pf
        self.domain_id =  domain_id
        #establish file neames

class ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTHierarchy
    _fieldinfo_fallback = ARTFieldInfo
    _fieldinfo_known = KnownARTFields

    def __init__(self,filename,data_style='art',
                 fields = None, storage_filename = None
                 skip_particles=False,skip_stars=False,
                 limit_level=None,spread_age=True):
        if fields is None:
            fields = fluid_fields
        self._fields_in_file = fields
        self._find_files(filename)
        self.file_amr = filename
        self.parameter_filename = filename
        self.skip_particles = skip_particles
        self.skip_stars = skip_stars
        self.limit_level = limit_level
        self.spread_age = spread_age
        StaticOutput.__init__(self,filename,data_style)
        self.storage_filename = storage_filename

    def __repr__(self):
        return self.basename.rsplit(".",1)[0]

    def _set_units(self):
        """
        Generates the conversion to various physical units based 
		on the parameters from the header
        """
        self.units = {}
        self.time_units = {}
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0

        #spatial units
        z   = self.current_redshift
        h   = self.hubble_constant
        boxcm_cal = self.parameters["boxh"]
        boxcm_uncal = boxcm_cal / h
        box_proper = boxcm_uncal/(1+z)
        aexpn = self["aexpn"]
        for unit in mpc_conversion:
            self.units[unit] = mpc_conversion[unit] * box_proper
            self.units[unit+'h'] = mpc_conversion[unit] * box_proper * h
            self.units[unit+'cm'] = mpc_conversion[unit] * boxcm_uncal
            self.units[unit+'hcm'] = mpc_conversion[unit] * boxcm_cal

        #all other units
        wmu = self.parameters["wmu"]
        Om0 = self.parameters['Om0']
        ng  = self.parameters['ng']
        wmu = self.parameters["wmu"]
        boxh   = self.parameters['boxh'] 
        aexpn  = self.parameters["aexpn"]
        hubble = self.parameters['hubble']

        cf = defaultdict(lambda: 1.0)
        r0 = boxh/ng
        P0= 4.697e-16 * Om0**2.0 * r0**2.0 * hubble**2.0
        T_0 = 3.03e5 * r0**2.0 * wmu * Om0 # [K]
        S_0 = 52.077 * wmu**(5.0/3.0)
        S_0 *= hubble**(-4.0/3.0)*Om0**(1.0/3.0)*r0**2.0
        #v0 =  r0 * 50.0*1.0e5 * np.sqrt(self.omega_matter)  #cm/s
        v0 = 50.0*r0*np.sqrt(Om0)
        t0 = r0/v0
        rho1 = 1.8791e-29 * hubble**2.0 * self.omega_matter
        rho0 = 2.776e11 * hubble**2.0 * Om0
        tr = 2./3. *(3.03e5*r0**2.0*wmu*self.omega_matter)*(1.0/(aexpn**2))     
        aM0 = rho0 * (boxh/hubble)**3.0 / ng**3.0
        cf['r0']=r0
        cf['P0']=P0
        cf['T_0']=T_0
        cf['S_0']=S_0
        cf['v0']=v0
        cf['t0']=t0
        cf['rho0']=rho0
        cf['rho1']=rho1
        cf['tr']=tr
        cf['aM0']=aM0

        #factors to multiply the native code units to CGS
        cf['Pressure'] = P0 #already cgs
        cf['Velocity'] = v0/aexpn*1.0e5 #proper cm/s
        cf["Mass"] = aM0 * 1.98892e33
        cf["Density"] = rho1*(aexpn**-3.0)
        cf["GasEnergy"] = rho0*v0**2*(aexpn**-5.0)
        cf["Potential"] = 1.0
        cf["Entropy"] = S_0
        cf["Temperature"] = tr
        self.cosmological_simulation = True
        self.conversion_factors = cf
        
        for particle_field in particle_fields:
            self.conversion_factors[particle_field] =  1.0
        for ax in 'xyz':
            self.conversion_factors["%s-velocity" % ax] = 1.0
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

    def _parse_parameter_file(self):
        """
        Get the various simulation parameters & constants.
        """
        self.dimensionality = 3
        self.refine_by = 2
        self.cosmological_simulation = True
        self.parameters = {}
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self.parameters.update(constants)
        with open(self.file_amr,'rb') as f:
            amr_header_vals = _read_struct(f,amr_header_struct)
            for to_skip in ['tl','dtl','tlold','dtlold','iSO']:
                _skip_record(f)
            (self.ncell,) = struct.unpack('>l', _read_record(f))
            # Try to figure out the root grid dimensions
            est = int(np.rint(self.ncell**(1.0/3.0)))
            # Note here: this is the number of *cells* on the root grid.
            # This is not the same as the number of Octs.
            self.domain_dimensions = np.ones(3, dtype='int64')*est 
            self.root_grid_mask_offset = f.tell()
            root_cells = self.domain_dimensions.prod()
            self.root_iOctCh = _read_frecord(f,'>i')[:root_cells]
            self.root_iOctCh = self.root_iOctCh.reshape(self.domain_dimensions,
                 order='F')
            self.root_grid_offset = f.tell()
            _skip_record(f) # hvar
            _skip_record(f) # var
            self.iOctFree, self.nOct = struct.unpack('>ii', _read_record(f))
            self.child_grid_offset = f.tell()
        self.parameters.update(amr_header_vals)
        if not self.skip_particles and self.file_particle_header:
            with open(self.file_particle_header,"rb") as fh:
                particle_header_vals = _read_struct(fh,particle_header_struct)
                fh.seek(seek_extras)
                n = particle_header_vals['Nspecies']
                wspecies = np.fromfile(fh,dtype='>f',count=10)
                lspecies = np.fromfile(fh,dtype='>i',count=10)
            self.parameters['wspecies'] = wspecies[:n]
            self.parameters['lspecies'] = lspecies[:n]
            ls_nonzero = np.diff(lspecies)[:n-1]
            mylog.info("Discovered %i species of particles",len(ls_nonzero))
            mylog.info("Particle populations: "+'%1.1e '*len(ls_nonzero),
                *ls_nonzero)
            for k,v in particle_header_vals.items():
                if k in self.parameters.keys():
                    if not self.parameters[k] == v:
                        mylog.info("Inconsistent parameter %s %1.1e  %1.1e",k,v,
                                   self.parameters[k])
                else:
                    self.parameters[k]=v
            self.parameters_particles = particle_header_vals
    
        #setup standard simulation yt expects to see
        self.current_redshift = self.parameters["aexpn"]**-1.0 - 1.0
        self.omega_lambda = amr_header_vals['Oml0']
        self.omega_matter = amr_header_vals['Om0']
        self.hubble_constant = amr_header_vals['hubble']
        self.min_level = amr_header_vals['min_level']
        self.max_level = amr_header_vals['max_level']
        self.hubble_time  = 1.0/(self.hubble_constant*100/3.08568025e19)
        self.current_time = b2t(self.parameters['t']) * sec_per_Gyr

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """
        Defined for the NMSU file naming scheme.
        This could differ for other formats.
        """
        fn = ("%s" % (os.path.basename(args[0])))
        f = ("%s" % args[0])
        prefix, suffix = filename_pattern['amr'].split('%s')
        if fn.endswith(suffix) and fn.startswith(prefix) and\
                os.path.exists(f): 
                return True
        return False

