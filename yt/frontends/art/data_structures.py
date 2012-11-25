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
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from .definitions import *
from io import _read_struct
from io import _read_art_level_info
from io import _read_record
from io import _read_frecord
from io import _skip_record
from io import _count_art_octs
from io import b2t
from io import load_level

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
from yt.geometry.oct_container import \
    RAMSESOctreeContainer

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
        assert domain_id > 0
        with open(self.pf.file_amr,"rb") as f:
            (self.nhydro_vars, self.iNOLL, self.level_oct_offsets,
                self.level_offsets) = \
                    _count_art_octs(f,self.pf.child_grid_offset,
                        self.pf.min_level,self.pf.max_level)
            self.level_info = [self.pf.domain_dimensions.prod()/8]
            for level in range(1,self.pf.max_level+1):
                le, fl,locts,root_level = _read_art_level_info(f, 
                    self.level_oct_offsets,level)
                self.level_info.append(locts)
                assert locts == le.shape[0]
        self.level_offsets[0] = self.pf.root_grid_offset
        self.local_oct_count  = int(np.sum(self.level_info))

    def _read_amr(self, oct_handler):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain 
           are needed - it's position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        #jump through hoops to get the root grid added indpendently
        #root grid is 128^3 *cells* so 64 octs
        NX = self.pf.domain_dimensions/2
        LE = self.pf.domain_left_edge
        RE = self.pf.domain_right_edge
        root_dx = (RE - LE) / NX
        #the zero level is not rounded to the nearest integer
        #so we adjust it slightly
        LL = LE + root_dx/2.0
        RL = RE - root_dx/2.0
        rpos = np.mgrid[ LL[0]:RL[0]:NX[0]*1j,
                         LL[1]:RL[1]:NX[1]*1j,
                         LL[2]:RL[2]:NX[2]*1j ]
        rpos = np.vstack([p.ravel() for p in rpos]).T
        nassigned = oct_handler.add(1, 0, rpos.shape[0], rpos, 1)
        assert(oct_handler.nocts == rpos.shape[0])
        assert(nassigned == rpos.shape[0])
        #now add the octs on other levels
        prev_nocts = oct_handler.nocts
        with open(self.pf.file_amr,"rb") as f:
            for level in range(1,self.pf.max_level+1):
                le, fl,locts,root_level = _read_art_level_info(f, 
                    self.level_oct_offsets,level)
                #fle = le/(NX*2.0**(level+1))
                fle = le/(NX*2.0**(level+1))
                #fc  = fle + 0.5/(NX*2**(level))
                #note that because we're adapting to the RAMSES 
                #architecture, cpu=0 and domain id will be fixed as
                #there is only 1 domain
                oct_handler.add(1,level,-1,fle,self.domain_id)
                assert oct_handler.nocts - prev_nocts == fle.shape[0]
                prev_nocts = oct_handler.nocts
        assert oct_handler.nocts == self.local_oct_count

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

class ARTDomainSubset(object):
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
        offset = 0
        for field in fields:
            tr[field] = np.zeros(self.cell_count, 'float64')
        for level in range(0,self.domain.pf.max_level+1):
            no = self.domain.level_info[level]
            temp = {}
            #now fill temp for the requested fields
            #request the cell data (not oct)
            hydro_data = load_level(self.domain.pf.file_amr,
                    self.domain.level_offsets,
                    self.domain.level_info, level, 
                    self.domain.nhydro_vars)
            for field in fields:
                temp[field] = np.empty((no, 8), dtype="float64")
            for field in fields:
                fi = tuple(fluid_fields).index(field)
                t = hydro_data[fi]
                nc=t.shape[0]
                assert nc%8==0
                assert nc/8==no
                for ci in range(8):
                    #t = np.reshape(t,(nc/8,8)).astype('float64')
                    temp[field][:,ci] = t[ci::8]
                    #np.reshape(t,(nc/8,8)).astype('float64')
            offset += oct_handler.fill_level(self.domain.domain_id, level,
                                   tr, temp, self.mask, offset)
        return tr

class ARTGeometryHandler(OctreeGeometryHandler):
    def __init__(self, pf, data_style='art'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.max_level = pf.max_level
        self.float_type = np.float64
        super(ARTGeometryHandler, self).__init__(pf, data_style)

    def _initialize_oct_handler(self):
        #mirroring ramses but we only one domain
        self.domains = [ARTDomainFile(self.parameter_file, 1)]
        total_octs = sum(dom.local_oct_count for dom in self.domains)
        self.num_grids = total_octs
        #this merely allocates space for the oct tree
        #and nothing else
        #yes, we use the RAMSES oct container
        self.oct_handler = RAMSESOctreeContainer(
            self.parameter_file.domain_dimensions/2,
            self.parameter_file.domain_left_edge,
            self.parameter_file.domain_right_edge)
        self.oct_handler.allocate_domains(
            [dom.local_oct_count for dom in self.domains])
        mylog.info("Allocated %s octs", total_octs)
        #this actually reads every oct and loads it into the octree
        for dom in self.domains:
            dom._read_amr(self.oct_handler)

    def _detect_fields(self):
        self.particle_field_list = particle_fields
        self.field_list = fluid_fields + particle_fields
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            mask = dobj.selector.select_octs(self.oct_handler)
            counts = self.oct_handler.count_cells(dobj.selector, mask)
            subsets = [ARTDomainSubset(d, mask, c)
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
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], subset.cell_count)

class ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTGeometryHandler
    _fieldinfo_fallback = ARTFieldInfo
    _fieldinfo_known = KnownARTFields
    
    def __init__(self, file_amr, storage_filename = None,
            skip_particles=False,skip_stars=False,data_style='art'):
        self.data_style = data_style
        self._find_files(file_amr)
        self.skip_particles = skip_particles
        self.skip_stars = skip_stars
        self.file_amr = file_amr
        self.parameter_filename = file_amr
        self.domain_left_edge  = np.zeros(3,dtype='float64')
        self.domain_right_edge = np.ones(3,dtype='float64') 
        StaticOutput.__init__(self, file_amr, data_style)
        self.storage_filename = storage_filename

    def _find_files(self,file_amr):
        """
        Given the AMR base filename, attempt to find the
        particle header, star files, etc.
        """
        prefix,suffix = filename_pattern['amr'].split('%s')
        affix = os.path.basename(file_amr).replace(prefix,'')
        affix = affix.replace(suffix,'')
        affix = affix.replace('_','')
        affix = affix[1:-1]
        dirname = os.path.dirname(file_amr)
        for filetype, pattern in filename_pattern.items():
            #sometimes the affix is surrounded by an extraneous _
            #so check for an extra character on either side
            check_filename = dirname+'/'+pattern%('?%s?'%affix)
            filenames = glob.glob(check_filename)
            if len(filenames)==1:
                setattr(self,"file_"+filetype,filenames[0])
                mylog.info('discovered %s',filetype)
            elif len(filenames)>1:
                setattr(self,"file_"+filetype,None)
                mylog.info("Ambiguous number of files found for %s",
                        check_filename)
            else:
                setattr(self,"file_"+filetype,None)

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
        self.units['unitary'] = 1.0
        self._parse_parameter_file()

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

        r0 = boxh/ng
        P0= 4.697e-16 * Om0**2.0 * r0**2.0 * hubble**2.0
        T_0 = 3.03e5 * r0**2.0 * wmu * Om0 # [K]
        S_0 = 52.077 * wmu**(5.0/3.0)
        S_0 *= hubble**(-4.0/3.0)*Om0**(1.0/3.0)*r0**2.0
        v0 =  r0 * 50.0*1.0e5 * np.sqrt(self.omega_matter)  #cm/s
        t0 = r0/v0
        rho0 = 1.8791e-29 * hubble**2.0 * self.omega_matter
        tr = 2./3. *(3.03e5*r0**2.0*wmu*self.omega_matter)*(1.0/(aexpn**2))     
        aM0 = rho0 * (boxh/hubble)**3.0 / ng**3.0

        #factors to multiply the native code units to CGS
        cf = defaultdict(lambda: 1.0)
        cf['Pressure'] = P0 #already cgs
        cf['Velocity'] = v0*1e3 #km/s -> cm/s
        cf["Mass"] = aM0 * 1.98892e33
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
        if self.file_particle_header is None:
            with open(self.file_particle_header,"rb") as fh:
                particle_header_vals = _read_struct(fh,particle_header_struct)
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
            self.parameters.update(particle_header_vals)
    
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
        if fn.endswith(".d") and fn.startswith('10Mpc') and\
                os.path.exists(f): 
                return True
        return False
