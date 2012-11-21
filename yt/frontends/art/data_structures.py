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
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from .definitions import *
from .io import *
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
	#This is a bit silly since NMSU ART is single-domain
    #and therefore this class could be merged with ARTStaticOuput
	#But we'll stay analogous to ramses and keep the AMR
	#info here and the rest of the simulation info in 
	#static output
    _last_mask = None
    _last_selector_id = None
    nvar = 6
    def __init__(self, pf, domain_id):
        self.pf = pf
        self.domain_id = domain_id
        (self.nhydro_vars, self.level_info, self.level_oct_offsets,
            self.level_offsets) = \
                _count_art_octs(self.pf.file_amr,self.pf.child_grid_offset,
                        self.pf.level_min,self.pf.level_max)
        self.level_offsets[0] = self.pf.root_grid_offset
    
    def _read_amr(self, oct_handler):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain 
           are needed - it's position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        self.nocts = self.pf.domain_dimensions.prod()
        with open(self.pf.file_amr,"rb") as f:
            for level in range(self.pf.max_level):
                le, locts= _read_art_level_info(f, self.level_oct_offsets,level)
                self.nocts += locts 
                #note that because we're adapting to the RAMSES 
                #architecture, cpu=0 and domain id will be fixed as
                #there is only 1 domain
                oct_handler.add(0,level,loct,le,self.domain_id)
        mylog.info("read in %1.1e octs",self.nocts)

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

class RAMSESDomainSubset(object):
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

    def fill(self, content, fields):
        #we are given a file (content)
        #and we are to fill it with the requested fields
        oct_handler = self.oct_handler
        fields = [f for ft, f in fields]
        min_level = self.domain.pf.min_level
        tr = {}
        filled = pos = offset = 0
        for field in fields:
            tr[field] = np.zeros(self.cell_count, 'float64')
        for level in range(self.pf.max_level):
            nc = self.domain.level_count[level]
            temp = {}
            for field in all_fields:
                temp[field] = np.empty((nc, 8), dtype="float64")
            #now fill temp for the requested fields
            hydro_data = load_level(self.domain.pf.file_amr,
                    self.domain.level_offsets,
                    self.domain.level_info, level, 
                    self.domain.nhydro_vars)
            offset += oct_handler.fill_level(self.domain.domain_id, level,
                                   tr, temp, self.mask, offset)
        return tr


def ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTGeometryHandler
    _fieldinfo_fallback = ARTFieldInfo
    _fieldinfo_known = KnownARTFields
    
    def __init__(self, file_amr, storage_filename = None,
            skip_particles=False,skip_stars=False):
        self._find_files(file_amr)
        self.skip_particles = skip_particles
        self.skip_stars = skip_stars
        self.file_amr = file_amr
        StaticOutput.__init__(self, file_amr, data_style)

    def _find_files(self,file_amr):
        """
        Given the AMR base filename, attempt to find the
        particle header, star files, etc.
        """
        prefix,suffix = filename_pattern['file_amr'].split('%s')
        affix = file_amr.replace(prefix,'')
        affix = affix.replace(suffix,'')
        affix = affix.replace('_','')
        for filetype, pattern in filename_pattern.items():
            #sometimes the affix is surrounded by an extraneous _
            #so check for an extra character on either side
            check_filename = pattern%('?%s?'%affix)
            filenames = glob.glob(check_filename)
            if len(filenames)==1:
                setattr(self,"file_"+filetype,filename)
                mylog.info('discovered %s',filetype)
            else:
                setattr(self,"file_"+filetype,None)
                mylog.info("Ambiguous number of files found for %s",
                        check_filename)
    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
        
    def _set_units(self):
        """
        Generates the conversion to various physical units based 
		on the parameters from the header
        """
        self.units = {}
        self.time_units = {}
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge 
									 - self.domain_left_edge).max()
		self._parse_parameter_file()

        #spatial units
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
        z   = self.current_redshift
        h   = self.hubble_constant
        wmu = self.parameters["wmu"]
        Om0 = self.parameters['Om0']
        ng  = self.parameters['ng']
        wmu = self.parameters["wmu"]
        boxh   = self.parameters['boxh'] 
        aexpn  = self.parameters["aexpn"]
        hubble = self.parameters['hubble']

        r0 = boxh/ng
        P0= 4.697e-16 * Om0**2.0 * r0**2.0 * hubble**2.0
        T_0 = 3.03e5 * r0**2.0 * wmu * Om0 # [K]
        S_0 = 52.077 * wmu**(5.0/3.0) * 
        S_0 *= hubble**(-4.0/3.0)*Om0**(1.0/3.0)*r0**2.0
        v0 =  r0 * 50.0*1.0e5 * np.sqrt(self.omega_matter)  #cm/s
        t0 = r0/v0
        rho0 = 1.8791e-29 * hubble**2.0 * self.omega_matter
        tr = 2./3. *(3.03e5*r0**2.0*wmu*self.omega_matter)*(1.0/(aexpn**2))     

        #factors to multiply the native code units to CGS
        cf = defaultdict(lambda: 1.0)
        cf['Pressure'] = P0 #already cgs
        cf['Velocity'] = v0*1e3 #km/s -> cm/s
        cf["Mass"] = self.parameters["aM0"] * 1.98892e33
        cf["Density"] = rho0*(aexpn**-3.0)
        cf["GasEnergy"] = rho0*v0**2*(aexpn**-5.0)
        cf["Potential"] = 1.0
        cf["Entropy"] = S_0
        cf["Temperature"] = tr
        self.cosmological_simulation = True
        self.conversion_factors = cf
        
        for ax in 'xyz':
            self.conversion_factors["%s-velocity" % ax] = v0/aexpn
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
        header_vals = {}
		with open(self.file_amr,'rb') as fh:
            amr_header_vals = _read_struct(fh,amr_header_struct)
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
        if self.file_particle_header is None:
            with open(self.file_particle_header,"rb") as fh:
                particle_header_vals = _read_struct(fh,header_struct)
                fh.seek(seek_extras)
                n = particle_header_vals['Nspecies']
                wspecies = np.fromfile(fh,dtype='>f',count=10)
                lspecies = np.fromfile(fh,dtype='>i',count=10)
            self.parameters['wspecies'] = wspecies[:n]
            self.parameters['lspecies'] = lspecies[:n]
            ls_nonzero = [ls for ls in self.parameters['lspecies'] if ls>0 ]
            mylog.info("Discovered %i species of particles",len(ls_nonzero))
            mylog.info("Particle populations: "+'%1.1e '*len(ls_nonzero),
                ls_nonzero)
        self.parameters.update(amr_header_vals)
        self.parameters.update(particle_header_vals)
        self.parameters.update(constants)
    
        #setup standard simulation yt expects to see
        self.current_redshift = self.parameters["aexpn"]**-1.0 - 1.0
        self.omega_lambda = header_vals['Oml0']
        self.omega_matter = header_vals['Om0']
        self.hubble_constant = header_vals['hubble']
        self.min_level = header_vals['min_level']
        self.max_level = header_vals['max_level']
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
        if fn.endswith(".d") and fn.startswith('10Mpc') and\
                os.path.exists(f): 
                return True
        return False
