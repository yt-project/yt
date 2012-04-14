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
import cPickle
import os
import struct

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
import yt.utilities.amr_utils as amr_utils

try:
    import yt.frontends.ramses._ramses_reader as _ramses_reader
except ImportError:
    _ramses_reader = None

from yt.utilities.physical_constants import \
    mass_hydrogen_cgs
    
from yt.frontends.art.definitions import art_particle_field_names

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

def num_deep_inc(f):
    def wrap(self, *args, **kwargs):
        self.num_deep += 1
        rv = f(self, *args, **kwargs)
        self.num_deep -= 1
        return rv
    return wrap

class ARTGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, hierarchy, level, locations, props,child_mask=None):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        start_index = props[0]
        self.Level = level
        self.Parent = []
        self.Children = []
        self.locations = locations
        self.start_index = start_index.copy()
        
        self.LeftEdge = props[0]
        self.RightEdge = props[1]
        self.ActiveDimensions = props[2] 
        #if child_mask is not None:
        #    self._set_child_mask(child_mask)

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.
        """
        if self.start_index != None:
            return self.start_index
        if len(self.Parent) == 0:
            start_index = self.LeftEdge / self.dds
            return na.rint(start_index).astype('int64').ravel()
        pdx = self.Parent[0].dds
        start_index = (self.Parent[0].get_global_startindex()) + \
                       na.rint((self.LeftEdge - self.Parent[0].LeftEdge)/pdx)
        self.start_index = (start_index*self.pf.refine_by).astype('int64').ravel()
        return self.start_index

    def __repr__(self):
        return "ARTGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ARTHierarchy(AMRHierarchy):

    grid = ARTGrid
    _handle = None
    
    def __init__(self, pf, data_style='art'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)
        self._setup_field_list()
        
    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        # This will need to be generalized to be used elsewhere.
        self.field_list = [ 'Density','TotalEnergy',
             'XMomentumDensity','YMomentumDensity','ZMomentumDensity',
             'Pressure','Gamma','GasEnergy',
             'MetalDensitySNII', 'MetalDensitySNIa',
             'PotentialNew','PotentialOld']
        self.field_list += art_particle_field_names

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        LEVEL_OF_EDGE = 7
        MAX_EDGE = (2 << (LEVEL_OF_EDGE- 1))
        
        min_eff = 0.30
        
        vol_max = 128**3
        
        f = open(self.pf.parameter_filename,'rb')
        
        
        (self.pf.nhydro_vars, self.pf.level_info,
        self.pf.level_oct_offsets, 
        self.pf.level_child_offsets) = \
                         _count_art_octs(f, 
                          self.pf.child_grid_offset,
                          self.pf.min_level, self.pf.max_level)
        self.pf.level_info[0]=self.pf.ncell
        self.pf.level_info = na.array(self.pf.level_info)        
        self.pf.level_offsets = self.pf.level_child_offsets
        self.pf.level_offsets = na.array(self.pf.level_offsets, dtype='int64')
        self.pf.level_offsets[0] = self.pf.root_grid_offset
        
        self.pf.level_art_child_masks = {}
        cm = self.pf.root_iOctCh>0
        cm_shape = (1,)+cm.shape 
        self.pf.level_art_child_masks[0] = cm.reshape(cm_shape).astype('uint8')        
        del cm
        
        root_psg = _ramses_reader.ProtoSubgrid(
                        na.zeros(3, dtype='int64'), # left index of PSG
                        self.pf.domain_dimensions, # dim of PSG
                        na.zeros((1,3), dtype='int64'), # left edges of grids
                        na.zeros((1,6), dtype='int64') # empty
                        )
        
        self.proto_grids = [[root_psg],]
        for level in xrange(1, len(self.pf.level_info)):
            if self.pf.level_info[level] == 0:
                self.proto_grids.append([])
                continue
            psgs = []
            effs,sizes = [], []

            if level > self.pf.limit_level : continue
            
            #refers to the left index for the art octgrid
            left_index, fl, nocts = _read_art_level_info(f, self.pf.level_oct_offsets,level)
            #left_index_gridpatch = left_index >> LEVEL_OF_EDGE
            
            #read in the child masks for this level and save them
            idc, art_child_mask = _read_child_mask_level(f, self.pf.level_child_offsets,
                level,nocts*8,nhydro_vars=self.pf.nhydro_vars)
            art_child_mask = art_child_mask.reshape((nocts,2,2,2))
            self.pf.level_art_child_masks[level]=art_child_mask
            #child_mask is zero where child grids exist and
            #thus where higher resolution data is available
            
            
            #compute the hilbert indices up to a certain level
            #the indices will associate an oct grid to the nearest
            #hilbert index?
            base_level = int( na.log10(self.pf.domain_dimensions.max()) /
                              na.log10(2))
            hilbert_indices = _ramses_reader.get_hilbert_indices(
                                    level + base_level, left_index)
            #print base_level, hilbert_indices.max(),
            hilbert_indices = hilbert_indices >> base_level + LEVEL_OF_EDGE
            #print hilbert_indices.max()
            
            # Strictly speaking, we don't care about the index of any
            # individual oct at this point.  So we can then split them up.
            unique_indices = na.unique(hilbert_indices)
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
                #for idomain in na.unique(ddfl[:,1]):
                #dom_ind = ddfl[:,1] == idomain
                #dleft_index = ddleft_index[dom_ind,:]
                #dfl = ddfl[dom_ind,:]
                
                dleft_index = ddleft_index
                dfl = ddfl
                initial_left = na.min(dleft_index, axis=0)
                idims = (na.max(dleft_index, axis=0) - initial_left).ravel()+2
                #this creates a grid patch that doesn't cover the whole level
                #necessarily, but with other patches covers all the regions
                #with octs. This object automatically shrinks its size
                #to barely encompass the octs inside of it.
                psg = _ramses_reader.ProtoSubgrid(initial_left, idims,
                                dleft_index, dfl)
                if psg.efficiency <= 0: continue
                
                #because grid patches may still be mostly empty, and with octs
                #that only partially fill the grid,it  may be more efficient
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
            eff_mean = na.mean(psg_eff)
            eff_nmin = na.sum([e<=min_eff*tol for e in psg_eff])
            eff_nall = len(psg_eff)
            mylog.info("Average subgrid efficiency %02.1f %%",
                        eff_mean*100.0)
            mylog.info("%02.1f%% (%i/%i) of grids had minimum efficiency",
                        eff_nmin*100.0/eff_nall,eff_nmin,eff_nall)
            
        
            mylog.debug("Done with level % 2i", level)
            pbar.finish()
            self.proto_grids.append(psgs)
            #print sum(len(psg.grid_file_locations) for psg in psgs)
            #mylog.info("Final grid count: %s", len(self.proto_grids[level]))
            if len(self.proto_grids[level]) == 1: continue
        self.num_grids = sum(len(l) for l in self.proto_grids)
                    
            
            

    num_deep = 0

        
    def _parse_hierarchy(self):
        """ The root grid has no octs except one which is refined.
        Still, it is the size of 128 cells along a length.
        Ignore the proto subgrid created for the root grid - it is wrong.
        """
        grids = []
        gi = 0
        
        for level, grid_list in enumerate(self.proto_grids):
            #The root level spans [0,2]
            #The next level spans [0,256]
            #The 3rd Level spans up to 128*2^3, etc.
            #Correct root level to span up to 128
            correction=1L
            if level == 0:
                correction=64L
            for g in grid_list:
                fl = g.grid_file_locations
                props = g.get_properties()*correction
                dds = ((2**level) * self.pf.domain_dimensions).astype("float64")
                self.grid_left_edge[gi,:] = props[0,:] / dds
                self.grid_right_edge[gi,:] = props[1,:] / dds
                self.grid_dimensions[gi,:] = props[2,:]
                self.grid_levels[gi,:] = level
                child_mask = na.zeros(props[2,:],'uint8')
                amr_utils.fill_child_mask(fl,props[0],
                    self.pf.level_art_child_masks[level],
                    child_mask)
                grids.append(self.grid(gi, self, level, fl, 
                    props*na.array(correction).astype('int64')))
                gi += 1
        self.grids = na.empty(len(grids), dtype='object')
        

        if self.pf.file_particle_data:
            lspecies = self.pf.parameters['lspecies']
            wspecies = self.pf.parameters['wspecies']
            Nrow     = self.pf.parameters['Nrow']
            nstars = lspecies[-1]
            a = self.pf.parameters['aexpn']
            hubble = self.pf.parameters['hubble']
            ud  = self.pf.parameters['r0']*a/hubble #proper Mpc units
            uv  = self.pf.parameters['v0']/(a**1.0)*1.0e5 #proper cm/s
            um  = self.pf.parameters['aM0'] #mass units in solar masses
            um *= 1.989e33 #convert solar masses to grams 
            pbar = get_pbar("Loading Particles   ",5)
            self.pf.particle_position,self.pf.particle_velocity = \
                read_particles(self.pf.file_particle_data,nstars,Nrow)
            pbar.update(1)
            np = lspecies[-1]
            self.pf.particle_position   = self.pf.particle_position[:np]
            self.pf.particle_position  -= 1.0 #fortran indices start with 0
            pbar.update(2)
            self.pf.particle_position  /= self.pf.domain_dimensions #to unitary units (comoving)
            pbar.update(3)
            self.pf.particle_velocity   = self.pf.particle_velocity[:np]
            self.pf.particle_velocity  *= uv #to proper cm/s
            pbar.update(4)
            self.pf.particle_species    = na.zeros(np,dtype='uint8')
            self.pf.particle_mass       = na.zeros(np,dtype='float64')
            
            dist = self.pf['cm']/self.pf.domain_dimensions[0]
            self.pf.conversion_factors['particle_mass'] = um #solar mass in g
            self.pf.conversion_factors['particle_species'] = 1.0
            for ax in 'xyz':
                self.pf.conversion_factors['particle_velocity_%s'%ax] = 1.0
                #already in unitary units
                self.pf.conversion_factors['particle_position_%s'%ax] = 1.0 
            self.pf.conversion_factors['particle_creation_time'] =  31556926.0
            self.pf.conversion_factors['particle_metallicity_fraction']=1.0
            self.pf.conversion_factors['particle_index']=1.0
            
            
            a,b=0,0
            for i,(b,m) in enumerate(zip(lspecies,wspecies)):
                self.pf.particle_species[a:b] = i #particle type
                self.pf.particle_mass[a:b]    = m #mass in solar masses
                a=b
            pbar.finish()
            
            self.pf.particle_star_index = i
            
            if self.pf.file_star_data:
                nstars, mass, imass, tbirth, metallicity1, metallicity2 \
                     = read_stars(self.pf.file_star_data,nstars,Nrow)
                nstars = nstars[0] 
                if nstars > 0 :
                    n=min(1e2,len(tbirth))
                    pbar = get_pbar("Stellar Ages        ",n)
                    self.pf.particle_star_ages  = \
                        b2t(tbirth,n=n,logger=lambda x: pbar.update(x)).astype('float64')
                    self.pf.particle_star_ages *= 1.0e9
                    self.pf.particle_star_ages *= 365*24*3600 #to seconds
                    self.pf.particle_star_ages = self.pf.current_time-self.pf.particle_star_ages
                    pbar.finish()
                    self.pf.particle_star_metallicity1 = metallicity1
                    self.pf.particle_star_metallicity2 = metallicity2
                    self.pf.particle_star_mass_initial = imass*self.pf.parameters['aM0']
                    self.pf.particle_mass[-nstars:] = mass*self.pf.parameters['aM0']
            left = self.pf.particle_position.shape[0]
            pbar = get_pbar("Gridding  Particles ",left)
            pos = self.pf.particle_position.copy()
            #particle indices travel with the particle positions
            pos = na.vstack((na.arange(pos.shape[0]),pos.T)).T 
            max_level = min(self.pf.max_level,self.pf.limit_level)
            grid_particle_count = na.zeros((len(grids),1),dtype='int64')
            
            #grid particles at the finest level, removing them once gridded
            for level in range(max_level,self.pf.min_level-1,-1):
                lidx = self.grid_levels[:,0] == level
                for gi,gidx in enumerate(na.where(lidx)[0]): 
                    g = grids[gidx]
                    assert g is not None
                    le,re = self.grid_left_edge[g._id_offset],self.grid_right_edge[g._id_offset]
                    idx = na.logical_and(na.all(le < pos[:,1:],axis=1),
                                         na.all(re > pos[:,1:],axis=1))
                    np = na.sum(idx)                     
                    g.NumberOfParticles = np
                    grid_particle_count[gidx,0]=np
                    g.hierarchy.grid_particle_count = grid_particle_count
                    if np==0: 
                        g.particle_indices = []
                        #we have no particles in this grid
                    else:
                        fidx = pos[:,0][idx]
                        g.particle_indices = fidx.astype('int64')
                        pos = pos[~idx] #throw out gridded particles from future gridding
                    self.grids[gidx] = g
                    left -= np
                    pbar.update(left)
            pbar.finish()
            
        else:
            
            pbar = get_pbar("Finalizing grids ",len(grids))
            for gi, g in enumerate(grids): 
                self.grids[gi] = g
            pbar.finish()
            

    def _get_grid_parents(self, grid, LE, RE):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(LE, RE)
        mask[grid_ind] = True
        mask = na.logical_and(mask, (self.grid_levels == (grid.Level-1)).flat)
        return self.grids[mask]

    def _populate_grid_objects(self):
        for gi,g in enumerate(self.grids):
            parents = self._get_grid_parents(g,
                            self.grid_left_edge[gi,:],
                            self.grid_right_edge[gi,:])
            if len(parents) > 0:
                g.Parent.extend(parents.tolist())
                for p in parents: p.Children.append(g)
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    # def _populate_grid_objects(self):
    #     mask = na.empty(self.grids.size, dtype='int32')
    #     pb = get_pbar("Populating grids", len(self.grids))
    #     for gi,g in enumerate(self.grids):
    #         pb.update(gi)
    #         amr_utils.get_box_grids_level(self.grid_left_edge[gi,:],
    #                             self.grid_right_edge[gi,:],
    #                             g.Level - 1,
    #                             self.grid_left_edge, self.grid_right_edge,
    #                             self.grid_levels, mask)
    #         parents = self.grids[mask.astype("bool")]
    #         if len(parents) > 0:
    #             g.Parent.extend((p for p in parents.tolist()
    #                     if p.locations[0,0] == g.locations[0,0]))
    #             for p in parents: p.Children.append(g)
    #         # Now we do overlapping siblings; note that one has to "win" with
    #         # siblings, so we assume the lower ID one will "win"
    #         amr_utils.get_box_grids_level(self.grid_left_edge[gi,:],
    #                             self.grid_right_edge[gi,:],
    #                             g.Level,
    #                             self.grid_left_edge, self.grid_right_edge,
    #                             self.grid_levels, mask, gi)
    #         mask[gi] = False
    #         siblings = self.grids[mask.astype("bool")]
    #         if len(siblings) > 0:
    #             g.OverlappingSiblings = siblings.tolist()
    #         g._prepare_grid()
    #         g._setup_dx()
    #     pb.finish()
    #     self.max_level = self.grid_levels.max()

    def _setup_field_list(self):
        if self.parameter_file.use_particles:
            # We know which particle fields will exist -- pending further
            # changes in the future.
            for field in art_particle_field_names:
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
    _handle = None
    
    def __init__(self, filename, data_style='art',
                 storage_filename = None, 
                 file_particle_header=None, 
                 file_particle_data=None,
                 file_star_data=None,
                 discover_particles=False,
                 use_particles=True,
                 limit_level=None):
        import yt.frontends.ramses._ramses_reader as _ramses_reader
        
        
        dirn = os.path.dirname(filename)
        base = os.path.basename(filename)
        aexp = base.split('_')[2].replace('.d','')
        
        self.file_particle_header = file_particle_header
        self.file_particle_data = file_particle_data
        self.file_star_data = file_star_data
        
        if limit_level is None:
            self.limit_level = na.inf
        else:
            mylog.info("Using maximum level: %i",limit_level)
            self.limit_level = limit_level
        
        if discover_particles:
            if file_particle_header is None:
                loc = filename.replace(base,'PMcrd%s.DAT'%aexp)
                if os.path.exists(loc):
                    self.file_particle_header = loc
                    mylog.info("Discovered particle header: %s",os.path.basename(loc))
            if file_particle_data is None:
                loc = filename.replace(base,'PMcrs0%s.DAT'%aexp)
                if os.path.exists(loc):
                    self.file_particle_data = loc
                    mylog.info("Discovered particle data:   %s",os.path.basename(loc))
            if file_star_data is None:
                loc = filename.replace(base,'stars_%s.dat'%aexp)
                if os.path.exists(loc):
                    self.file_star_data = loc
                    mylog.info("Discovered stellar data:    %s",os.path.basename(loc))
        
        self.use_particles = any([self.file_particle_header,
            self.file_star_data, self.file_particle_data])
        StaticOutput.__init__(self, filename, data_style)
        
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'art'
        self.parameters["Time"] = 1. # default unit is 1...
        self.parameters["InitialTime"]=self.current_time
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
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        
        
        z = self.current_redshift
        
        h = self.hubble_constant
        boxcm_cal = self["boxh"]
        boxcm_uncal = boxcm_cal / h
        box_proper = boxcm_uncal/(1+z)
        aexpn = self["aexpn"]
        for unit in mpc_conversion:
            self.units[unit] = mpc_conversion[unit] * box_proper
            self.units[unit+'h'] = mpc_conversion[unit] * box_proper * h
            self.units[unit+'cm'] = mpc_conversion[unit] * boxcm_uncal
            self.units[unit+'hcm'] = mpc_conversion[unit] * boxcm_cal
        # Variable names have been chosen to reflect primary reference
        #Om0 = self["Om0"]
        #boxh = self["boxh"]
        wmu = self["wmu"]
        #ng = self.domain_dimensions[0]
        #r0 = self["cmh"]/ng # comoving cm h^-1
        #t0 = 6.17e17/(self.hubble_constant + na.sqrt(self.omega_matter))
        #v0 = r0 / t0
        #rho0 = 1.8791e-29 * self.hubble_constant**2.0 * self.omega_matter
        #e0 = v0**2.0
        
        wmu = self["wmu"]
        boxh = self["boxh"]
        aexpn = self["aexpn"]
        hubble = self.hubble_constant
        ng = self.domain_dimensions[0]
        self.r0 = boxh/ng
        self.v0 =  self.r0 * 50.0*1.0e5 * na.sqrt(self.omega_matter)  #cm/s
        self.t0 = self.r0/self.v0
        # this is 3H0^2 / (8pi*G) *h*Omega0 with H0=100km/s. 
        # ie, critical density 
        self.rho0 = 1.8791e-29 * hubble**2.0 * self.omega_matter
        self.tr = 2./3. *(3.03e5*self.r0**2.0*wmu*self.omega_matter)*(1.0/(aexpn**2))     
        tr  = self.tr
        
        #factors to multiply the native code units to CGS
        self.conversion_factors['Pressure'] = self.parameters["P0"] #already cgs
        self.conversion_factors['Velocity'] = self.parameters['v0']*1e3 #km/s -> cm/s
        self.conversion_factors["Mass"] = self.parameters["aM0"] * 1.98892e33
        self.conversion_factors["Density"] = self.rho0*(aexpn**-3.0)
        self.conversion_factors["GasEnergy"] = self.rho0*self.v0**2*(aexpn**-5.0)
        self.conversion_factors["Temperature"] = tr
        self.conversion_factors["Potential"] = 1.0
        self.cosmological_simulation = True
        
        # Now our conversion factors
        for ax in 'xyz':
            # Add on the 1e5 to get to cm/s
            self.conversion_factors["%s-velocity" % ax] = self.v0/aexpn
        seconds = self.t0
        self.time_units['Gyr']   = 1.0/(1.0e9*365*3600*24.0)
        self.time_units['Myr']   = 1.0/(1.0e6*365*3600*24.0)
        self.time_units['years'] = 1.0/(365*3600*24.0)
        self.time_units['days']  = 1.0 / (3600*24.0)


        #we were already in seconds, go back in to code units
        #self.current_time /= self.t0 
        #self.current_time = b2t(self.current_time,n=1)
        
    
    def _parse_parameter_file(self):
        # We set our domain to run from 0 .. 1 since we are otherwise
        # unconstrained.
        self.domain_left_edge = na.zeros(3, dtype="float64")
        self.domain_right_edge = na.ones(3, dtype="float64")
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self.parameters = {}

        header_struct = [
            ('>i','pad byte'),
            ('>256s','jname'),
            ('>i','pad byte'),
            
            ('>i','pad byte'),
            ('>i','istep'),
            ('>d','t'),
            ('>d','dt'),
            ('>f','aexpn'),
            ('>f','ainit'),
            ('>i','pad byte'),
            
            ('>i','pad byte'),
            ('>f','boxh'),
            ('>f','Om0'),
            ('>f','Oml0'),
            ('>f','Omb0'),
            ('>f','hubble'),
            ('>i','pad byte'),
            
            ('>i','pad byte'),
            ('>i','nextras'),
            ('>i','pad byte'),

            ('>i','pad byte'),
            ('>f','extra1'),
            ('>f','extra2'),
            ('>i','pad byte'),

            ('>i','pad byte'),
            ('>256s','lextra'),
            ('>256s','lextra'),
            ('>i','pad byte'),
            
            ('>i', 'pad byte'),
            ('>i', 'min_level'),
            ('>i', 'max_level'),
            ('>i', 'pad byte'),
            ]
        
        f = open(self.parameter_filename, "rb")
        header_vals = {}
        for format, name in header_struct:
            size = struct.calcsize(format)
            # We parse single values at a time, so this will
            # always need to be indexed with 0
            output = struct.unpack(format, f.read(size))[0]
            header_vals[name] = output
        self.dimensionality = 3 # We only support three
        self.refine_by = 2 # Octree
        # Update our parameters with the header and with some compile-time
        # constants we will set permanently.
        self.parameters.update(header_vals)
        self.parameters["Y_p"] = 0.245
        self.parameters["wmu"] = 4.0/(8.0-5.0*self.parameters["Y_p"])
        self.parameters["gamma"] = 5./3.
        self.parameters["T_CMB0"] = 2.726  
        self.parameters["T_min"] = 300.0 #T floor in K
        self.parameters["boxh"] = header_vals['boxh']
        self.parameters['ng'] = 128 # of 0 level cells in 1d 
        self.current_redshift = self.parameters["aexpn"]**-1.0 - 1.0
        self.parameters['CosmologyInitialRedshift']=self.current_redshift
        self.data_comment = header_vals['jname']
        self.current_time_raw = header_vals['t']
        self.current_time = header_vals['t']
        self.omega_lambda = header_vals['Oml0']
        self.omega_matter = header_vals['Om0']
        self.hubble_constant = header_vals['hubble']
        self.min_level = header_vals['min_level']
        self.max_level = header_vals['max_level']
        self.nhydro_vars = 10 #this gets updated later, but we'll default to this
        #nchem is nhydrovars-8, so we typically have 2 extra chem species 
        self.hubble_time  = 1.0/(self.hubble_constant*100/3.08568025e19)
        #self.hubble_time /= 3.168876e7 #Gyr in s 
        # def integrand(x,oml=self.omega_lambda,omb=self.omega_matter):
        #     return 1./(x*na.sqrt(oml+omb*x**-3.0))
        # spacings = na.logspace(-5,na.log10(self.parameters['aexpn']),1e5)
        # integrand_arr = integrand(spacings)
        # self.current_time = na.trapz(integrand_arr,dx=na.diff(spacings))
        # self.current_time *= self.hubble_time
        self.current_time = b2t(self.current_time_raw)*1.0e9*365*3600*24         
        for to_skip in ['tl','dtl','tlold','dtlold','iSO']:
            _skip_record(f)

        
        Om0 = self.parameters['Om0']
        hubble = self.parameters['hubble']
        dummy = 100.0 * hubble * na.sqrt(Om0)
        ng = self.parameters['ng']
        wmu = self.parameters["wmu"]
        boxh = header_vals['boxh'] 
        
        #distance unit #boxh is units of h^-1 Mpc
        self.parameters["r0"] = self.parameters["boxh"] / self.parameters['ng']
        r0 = self.parameters["r0"]
        #time, yrs
        self.parameters["t0"] = 2.0 / dummy * 3.0856e19 / 3.15e7
        #velocity velocity units in km/s
        self.parameters["v0"] = 50.0*self.parameters["r0"]*\
                na.sqrt(self.parameters["Om0"])
        #density = 3H0^2 * Om0 / (8*pi*G) - unit of density in Msun/Mpc^3
        self.parameters["rho0"] = 2.776e11 * hubble**2.0 * Om0
        rho0 = self.parameters["rho0"]
        #Pressure = rho0 * v0**2 - unit of pressure in g/cm/s^2
        self.parameters["P0"] = 4.697e-16 * Om0**2.0 * r0**2.0 * hubble**2.0
        #T_0 = unit of temperature in K and in keV)
        #T_0 = 2.61155 * r0**2 * wmu * Om0 ! [keV]
        self.parameters["T_0"] = 3.03e5 * r0**2.0 * wmu * Om0 # [K]
        #S_0 = unit of entropy in keV * cm^2
        self.parameters["S_0"] = 52.077 * wmu**(5.0/3.0) * hubble**(-4.0/3.0)*Om0**(1.0/3.0)*r0**2.0
        
        #mass conversion (Mbox = rho0 * Lbox^3, Mbox_code = Ng^3
        #     for non-cosmological run aM0 must be defined during initialization
        #     [aM0] = [Msun]
        self.parameters["aM0"] = rho0 * (boxh/hubble)**3.0 / ng**3.0
        
        #CGS for everything in the next block
    
        (self.ncell,) = struct.unpack('>l', _read_record(f))
        # Try to figure out the root grid dimensions
        est = int(na.rint(self.ncell**(1.0/3.0)))
        # Note here: this is the number of *cells* on the root grid.
        # This is not the same as the number of Octs.
        self.domain_dimensions = na.ones(3, dtype='int64')*est 

        self.root_grid_mask_offset = f.tell()
        #_skip_record(f) # iOctCh
        root_cells = self.domain_dimensions.prod()
        self.root_iOctCh = _read_frecord(f,'>i')[:root_cells]
        self.root_iOctCh = self.root_iOctCh.reshape(self.domain_dimensions,order='F')
        self.root_grid_offset = f.tell()
        _skip_record(f) # hvar
        _skip_record(f) # var

        self.iOctFree, self.nOct = struct.unpack('>ii', _read_record(f))
        self.child_grid_offset = f.tell()

        f.close()
        
        if self.file_particle_header is not None:
            self._read_particle_header(self.file_particle_header)
        
    def _read_particle_header(self,fn):    
        """ Reads control information, various parameters from the 
            particle data set. Adapted from Daniel Ceverino's 
            Read_Particles_Binary in analysis_ART.F   
        """ 
        header_struct = [
            ('>i','pad'),
            ('45s','header'), 
            ('>f','aexpn'),
            ('>f','aexp0'),
            ('>f','amplt'),
            ('>f','astep'),

            ('>i','istep'),
            ('>f','partw'),
            ('>f','tintg'),

            ('>f','Ekin'),
            ('>f','Ekin1'),
            ('>f','Ekin2'),
            ('>f','au0'),
            ('>f','aeu0'),


            ('>i','Nrow'),
            ('>i','Ngridc'),
            ('>i','Nspecies'),
            ('>i','Nseed'),

            ('>f','Om0'),
            ('>f','Oml0'),
            ('>f','hubble'),
            ('>f','Wp5'),
            ('>f','Ocurv'),
            ('>f','Omb0'),
            ('>%ds'%(396),'extras'),
            ('>f','unknown'),

            ('>i','pad')]
        fh = open(fn,'rb')
        vals = _read_struct(fh,header_struct)
        
        for k,v in vals.iteritems():
            self.parameters[k]=v
        
        seek_extras = 137
        fh.seek(seek_extras)
        n = self.parameters['Nspecies']
        self.parameters['wspecies'] = na.fromfile(fh,dtype='>f',count=10)
        self.parameters['lspecies'] = na.fromfile(fh,dtype='>i',count=10)
        self.parameters['wspecies'] = self.parameters['wspecies'][:n]
        self.parameters['lspecies'] = self.parameters['lspecies'][:n]
        fh.close()
        
        ls_nonzero = [ls for ls in self.parameters['lspecies'] if ls>0 ]
        mylog.debug("Discovered %i species of particles",len(ls_nonzero))
        mylog.debug("Particle populations: "+'%1.1e '*len(ls_nonzero),ls_nonzero)
        

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return False # We make no effort to auto-detect ART data


