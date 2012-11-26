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
from io import read_particles
from io import read_stars
from io import spread_ages
from io import _count_art_octs
from io import _read_art_level_info
from io import _read_art_child
from io import _skip_record
from io import _read_record
from io import _read_frecord
from io import _read_record_size
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

class ARTGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, hierarchy, level, locations,start_index, le,re,gd,
            child_mask=None,nop=0):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        start_index =start_index 
        self.Level = level
        self.Parent = []
        self.Children = []
        self.locations = locations
        self.start_index = start_index.copy()
        
        self.LeftEdge = le
        self.RightEdge = re
        self.ActiveDimensions = gd
        self.NumberOfParticles=nop
        for particle_field in particle_fields:
            setattr(self,particle_field,np.array([]))

    def _setup_dx(self):
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = np.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] \
                = self.dds

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.
        """
        if self.start_index != None:
            return self.start_index
        if len(self.Parent) == 0:
            start_index = self.LeftEdge / self.dds
            return np.rint(start_index).astype('int64').ravel()
        pdx = self.Parent[0].dds
        start_index = (self.Parent[0].get_global_startindex()) + \
                       np.rint((self.LeftEdge - self.Parent[0].LeftEdge)/pdx)
        self.start_index = (start_index*self.pf.refine_by)\
                           .astype('int64').ravel()
        return self.start_index

    def __repr__(self):
        return "ARTGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ARTHierarchy(AMRHierarchy):
    grid = ARTGrid
    _handle = None
    
    def __init__(self, pf, data_style='art'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        self.max_level = pf.max_level
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.float_type = np.float64
        AMRHierarchy.__init__(self,pf,data_style)
        if not self.pf.skip_particles:
            self._setup_particle_grids()
        self._setup_field_list()

    def _initialize_data_storage(self):
        pass
    
    def _detect_fields(self):
        self.field_list = []
        self.field_list += fluid_fields
        self.field_list += particle_fields
        
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()
            
    def _count_grids(self):
        LEVEL_OF_EDGE = 7
        MAX_EDGE = (2 << (LEVEL_OF_EDGE- 1))
        min_eff = 0.30
        vol_max = 128**3
        with open(self.pf.parameter_filename,'rb') as f:
            (self.pf.nhydro_vars, self.pf.level_info,
            self.pf.level_oct_offsets, 
            self.pf.level_child_offsets) = \
                             _count_art_octs(f, 
                              self.pf.child_grid_offset,
                              self.pf.min_level, self.pf.max_level)
            self.pf.level_info[0]=self.pf.ncell
            self.pf.level_info = np.array(self.pf.level_info)
            self.pf.level_offsets = self.pf.level_child_offsets
            self.pf.level_offsets = np.array(self.pf.level_offsets, 
                                             dtype='int64')
            self.pf.level_offsets[0] = self.pf.root_grid_offset
            self.pf.level_art_child_masks = {}
            cm = self.pf.root_iOctCh>0
            cm_shape = (1,)+cm.shape 
            self.pf.level_art_child_masks[0] = \
                    cm.reshape(cm_shape).astype('uint8')        
            del cm
            root_psg = _ramses_reader.ProtoSubgrid(
                            np.zeros(3, dtype='int64'), # left index of PSG
                            self.pf.domain_dimensions, # dim of PSG
                            np.zeros((1,3), dtype='int64'),# left edges of grids
                            np.zeros((1,6), dtype='int64') # empty
                            )
            self.proto_grids = [[root_psg],]
            for level in xrange(1, len(self.pf.level_info)):
                if self.pf.level_info[level] == 0:
                    self.proto_grids.append([])
                    continue
                psgs = []
                effs,sizes = [], []
                if self.pf.limit_level:
                    if level > self.pf.limit_level : continue
                #refers to the left index for the art octgrid
                left_index, fl, nocts,root_level = _read_art_level_info(f, 
                        self.pf.level_oct_offsets,level,
                        coarse_grid=self.pf.domain_dimensions[0])
                if level>1:
                    assert root_level == last_root_level
                last_root_level = root_level
                #left_index_gridpatch = left_index >> LEVEL_OF_EDGE
                #read in the child masks for this level and save them
                idc, art_child_mask = _read_child_mask_level(f, 
                        self.pf.level_child_offsets,
                    level,nocts*8,nhydro_vars=self.pf.nhydro_vars)
                art_child_mask = art_child_mask.reshape((nocts,2,2,2))
                self.pf.level_art_child_masks[level]=art_child_mask
                #child_mask is zero where child grids exist and
                #thus where higher resolution data is available
                #compute the hilbert indices up to a certain level
                #the indices will associate an oct grid to the nearest
                #hilbert index?
                base_level = int( np.log10(self.pf.domain_dimensions.max()) /
                                  np.log10(2))
                hilbert_indices = _ramses_reader.get_hilbert_indices(
                                        level + base_level, left_index)
                #print base_level, hilbert_indices.max(),
                hilbert_indices = hilbert_indices >> base_level + LEVEL_OF_EDGE
                #print hilbert_indices.max()
                # Strictly speaking, we don't care about the index of any
                # individual oct at this point.  So we can then split them up.
                unique_indices = np.unique(hilbert_indices)
                mylog.info("Level % 2i has % 10i unique indices for %0.3e octs",
                            level, unique_indices.size, hilbert_indices.size)
                #use the hilbert indices to order oct grids so that consecutive
                #items on a list are spatially near each other
                #this is useful because we will define grid patches over these
                #octs, which are more efficient if the octs are spatially close
                #split into list of lists, with domains containing 
                #lists of sub octgrid left indices and an index
                #referring to the domain on which they live
                pbar = get_pbar("Calc Hilbert Indices ",1)
                locs, lefts = _ramses_reader.get_array_indices_lists(
                            hilbert_indices, unique_indices, left_index, fl)
                pbar.finish()
                #iterate over the domains    
                step=0
                pbar = get_pbar("Re-gridding  Level %i "%level, len(locs))
                psg_eff = []
                for ddleft_index, ddfl in zip(lefts, locs):
                    #iterate over just the unique octs
                    #why would we ever have non-unique octs?
                    #perhaps the hilbert ordering may visit the same
                    #oct multiple times - review only unique octs 
                    #for idomain in np.unique(ddfl[:,1]):
                    #dom_ind = ddfl[:,1] == idomain
                    #dleft_index = ddleft_index[dom_ind,:]
                    #dfl = ddfl[dom_ind,:]
                    dleft_index = ddleft_index
                    dfl = ddfl
                    initial_left = np.min(dleft_index, axis=0)
                    idims = (np.max(dleft_index, axis=0) - initial_left).ravel()
                    idims +=2
                    #this creates a grid patch that doesn't cover the whole leve
                    #necessarily, but with other patches covers all the regions
                    #with octs. This object automatically shrinks its size
                    #to barely encompass the octs inside of it.
                    psg = _ramses_reader.ProtoSubgrid(initial_left, idims,
                                    dleft_index, dfl)
                    if psg.efficiency <= 0: continue
                    #because grid patches maybe mostly empty, and with octs
                    #that only partially fill the grid, it may be more efficient
                    #to split large patches into smaller patches. We split
                    #if less than 10% the volume of a patch is covered with octs
                    if idims.prod() > vol_max or psg.efficiency < min_eff:
                        psg_split = _ramses_reader.recursive_patch_splitting(
                            psg, idims, initial_left, 
                            dleft_index, dfl,min_eff=min_eff,use_center=True,
                            split_on_vol=vol_max)
                        psgs.extend(psg_split)
                        psg_eff += [x.efficiency for x in psg_split] 
                    else:
                        psgs.append(psg)
                        psg_eff =  [psg.efficiency,]
                    tol = 1.00001
                    step+=1
                    pbar.update(step)
                eff_mean = np.mean(psg_eff)
                eff_nmin = np.sum([e<=min_eff*tol for e in psg_eff])
                eff_nall = len(psg_eff)
                mylog.info("Average subgrid efficiency %02.1f %%",
                            eff_mean*100.0)
                mylog.info("%02.1f%% (%i/%i) of grids had minimum efficiency",
                            eff_nmin*100.0/eff_nall,eff_nmin,eff_nall)
                
                mylog.info("Done with level % 2i; max LE %i", level,
                           np.max(left_index))
                pbar.finish()
                self.proto_grids.append(psgs)
                if len(self.proto_grids[level]) == 1: continue
        self.num_grids = sum(len(l) for l in self.proto_grids)
        
    def _parse_hierarchy(self):
        grids = []
        gi = 0
        dd=self.pf.domain_dimensions
        for level, grid_list in enumerate(self.proto_grids):
            dds = ((2**level) * dd).astype("float64")
            for g in grid_list:
                fl = g.grid_file_locations
                props = g.get_properties()
                start_index = props[0,:]
                le = props[0,:].astype('float64')/dds
                re = props[1,:].astype('float64')/dds
                gd = props[2,:].astype('int64')
                if level==0:
                    le = np.zeros(3,dtype='float64')
                    re = np.ones(3,dtype='float64')
                    gd = dd
                self.grid_left_edge[gi,:] = le
                self.grid_right_edge[gi,:] = re
                self.grid_dimensions[gi,:] = gd
                assert np.all(self.grid_left_edge[gi,:]<=1.0)    
                self.grid_levels[gi,:] = level
                child_mask = np.zeros(props[2,:],'uint8')
                amr_utils.fill_child_mask(fl,start_index,
                    self.pf.level_art_child_masks[level],
                    child_mask)
                grids.append(self.grid(gi, self, level, fl, 
                    start_index,le,re,gd))
                gi += 1
        self.grids = np.empty(len(grids), dtype='object')
        if not self.pf.skip_particles and self.pf.file_particle_data:
            lspecies = self.pf.parameters['lspecies']
            wspecies = self.pf.parameters['wspecies']
            um  = self.pf.conversion_factors['Mass'] #mass units in g
            uv  = self.pf.conversion_factors['Velocity'] #mass units in g
            self.pf.particle_position,self.pf.particle_velocity = \
                read_particles(self.pf.file_particle_data,
                        self.pf.parameters['Nrow'])
            nparticles = lspecies[-1]
            if not np.all(self.pf.particle_position[nparticles:]==0.0):
                mylog.info('WARNING: unused particles discovered from lspecies')
            self.pf.particle_position = self.pf.particle_position[:nparticles]
            self.pf.particle_velocity = self.pf.particle_velocity[:nparticles]
            self.pf.particle_position  /= self.pf.domain_dimensions 
            self.pf.particle_velocity   = self.pf.particle_velocity
            self.pf.particle_velocity  *= uv #to proper cm/s
            self.pf.particle_star_index = len(wspecies)-1
            self.pf.particle_type = np.zeros(nparticles,dtype='int')
            self.pf.particle_mass = np.zeros(nparticles,dtype='float32')
            a=0
            for i,(b,m) in enumerate(zip(lspecies,wspecies)):
                if i == self.pf.particle_star_index:
                    sa,sb = a,b
                self.pf.particle_type[a:b] = i #particle type
                self.pf.particle_mass[a:b] = m*um #mass in grams
                a=b
            if not self.pf.skip_stars and self.pf.file_particle_stars: 
                (nstars_rs,), mass, imass, tbirth, metallicity1, metallicity2, \
                        ws_old,ws_oldi,tdum,adum \
                     = read_stars(self.pf.file_particle_stars)
                self.pf.nstars_rs = nstars_rs     
                self.pf.nstars_pa = b-a
                inconsistent=self.pf.particle_type==self.pf.particle_star_index
                if not nstars_rs==np.sum(inconsistent):
                    mylog.info('WARNING!: nstars is inconsistent!')
                del inconsistent
                if nstars_rs > 0 :
                    n=min(1e2,len(tbirth))
                    birthtimes= b2t(tbirth,n=n)
                    birthtimes = birthtimes.astype('float64')
                    assert birthtimes.shape == tbirth.shape    
                    birthtimes*= 1.0e9 #from Gyr to yr
                    birthtimes*= 365*24*3600 #to seconds
                    ages = self.pf.current_time-birthtimes
                    spread = self.pf.spread_age
                    if type(spread)==type(5.5):
                        ages = spread_ages(ages,spread=spread)
                    elif spread:
                        ages = spread_ages(ages)
                    idx = self.pf.particle_type == self.pf.particle_star_index
                    for psf in particle_star_fields:
                        setattr(self.pf,psf,
                                np.zeros(nparticles,dtype='float32'))
                    self.pf.particle_age[sa:sb] = ages
                    self.pf.particle_mass[sa:sb] = mass
                    self.pf.particle_mass_initial[sa:sb] = imass
                    self.pf.particle_creation_time[sa:sb] = birthtimes
                    self.pf.particle_metallicity1[sa:sb] = metallicity1
                    self.pf.particle_metallicity2[sa:sb] = metallicity2
                    self.pf.particle_metallicity[sa:sb]  = metallicity1\
                                                          + metallicity2
        for gi,g in enumerate(grids):    
            self.grids[gi]=g
                    
    def _get_grid_parents(self, grid, LE, RE):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(LE, RE)
        mask[grid_ind] = True
        mask = np.logical_and(mask, (self.grid_levels == (grid.Level-1)).flat)
        return self.grids[mask]

    def _populate_grid_objects(self):
        mask = np.empty(self.grids.size, dtype='int32')
        pb = get_pbar("Populating grids", len(self.grids))
        for gi,g in enumerate(self.grids):
            pb.update(gi)
            amr_utils.get_box_grids_level(self.grid_left_edge[gi,:],
                                self.grid_right_edge[gi,:],
                                g.Level - 1,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            parents = self.grids[mask.astype("bool")]
            if len(parents) > 0:
                g.Parent.extend((p for p in parents.tolist()
                        if p.locations[0,0] == g.locations[0,0]))
                for p in parents: p.Children.append(g)
            #Now we do overlapping siblings; note that one has to "win" with
            #siblings, so we assume the lower ID one will "win"
            amr_utils.get_box_grids_level(self.grid_left_edge[gi,:],
                                self.grid_right_edge[gi,:],
                                g.Level,
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask, gi)
            mask[gi] = False
            siblings = self.grids[mask.astype("bool")]
            if len(siblings) > 0:
                g.OverlappingSiblings = siblings.tolist()
            g._prepare_grid()
            g._setup_dx()
            #instead of gridding particles assign them all to the root grid
            if gi==0:
                for particle_field in particle_fields:
                    source = getattr(self.pf,particle_field,None)
                    if source is None:
                        for i,ax in enumerate('xyz'):
                            pf = particle_field.replace('_%s'%ax,'')
                            source = getattr(self.pf,pf,None)
                            if source is not None:
                                source = source[:,i]
                                break
                    if source is not None:
                        mylog.info("Attaching %s to the root grid",
                                    particle_field)
                        g.NumberOfParticles = source.shape[0]
                        setattr(g,particle_field,source)
        pb.finish()
        self.max_level = self.grid_levels.max()

    def _setup_field_list(self):
        if not self.parameter_file.skip_particles:
            # We know which particle fields will exist -- pending further
            # changes in the future.
            for field in particle_fields:
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
                # Note that we call add_field on the field_info directly.  This
                # will allow the same field detection mechanism to work for 1D,
                # 2D and 3D fields.
                self.pf.field_info.add_field(field, NullFunc,
                                             convert_function=cf,
                                             take_log=True, particle_type=True)

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](
            self.pf.parameter_filename,
            self.pf.nhydro_vars,
            self.pf.level_info,
            self.pf.level_offsets)

class ARTStaticOutput(StaticOutput):
    _hierarchy_class = ARTHierarchy
    _fieldinfo_fallback = ARTFieldInfo
    _fieldinfo_known = KnownARTFields
    
    def __init__(self, file_amr, storage_filename = None,
            skip_particles=False,skip_stars=False,limit_level=None,
            spread_age=True,data_style='art'):
        self.data_style = data_style
        self._find_files(file_amr)
        self.skip_particles = skip_particles
        self.skip_stars = skip_stars
        self.file_amr = file_amr
        self.parameter_filename = file_amr
        self.limit_level = limit_level
        self.spread_age = spread_age
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
        cf['Velocity'] = v0/aexpn*1.0e5 #proper cm/s
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
        for particle_field in particle_fields:
            self.conversion_factors[particle_field] =  1.0
        self.conversion_factors['particle_creation_time'] =  31556926.0
        self.conversion_factors['Msun'] = 5.027e-34 

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

