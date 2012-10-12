"""
ART-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Christopher Erick Moody <cemoody@ucsc.edu>
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
from yt.utilities.lib import \
    get_box_grids_level
import yt.utilities.lib as amr_utils

import yt.frontends.ramses._ramses_reader as _ramses_reader #do not fail silently;

from yt.utilities.physical_constants import \
    mass_hydrogen_cgs
    
from yt.frontends.art.definitions import art_particle_field_names

from yt.frontends.art.io import _read_child_mask_level
from yt.frontends.art.io import read_particles
from yt.frontends.art.io import read_stars
from yt.frontends.art.io import spread_ages
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

    def __init__(self, id, hierarchy, level, locations,start_index, le,re,gd,
            child_mask=None,np=0):
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
        self.NumberOfParticles=np
        self.particle_type = na.array([])
        self.particle_id= na.array([])
        self.particle_age= na.array([])
        self.particle_position_x = na.array([])
        self.particle_position_y = na.array([])
        self.particle_position_z = na.array([])
        self.particle_velocity_x = na.array([])
        self.particle_velocity_y = na.array([])
        self.particle_velocity_z = na.array([])
        self.particle_mass= na.array([])
        self.star_position_x = na.array([])
        self.star_position_y = na.array([])
        self.star_position_z = na.array([])
        self.star_velocity_x = na.array([])
        self.star_velocity_y = na.array([])
        self.star_velocity_z = na.array([])
        self.star_age = na.array([])
        self.star_metallicity1 = na.array([])
        self.star_metallicity2 = na.array([])
        self.star_mass_initial = na.array([])
        self.star_mass = na.array([])

        self.field_dict = { 'particle_index': self.particle_id,
            'particle_type':self.particle_type,
            'particle_position_x':self.particle_position_x,
            'particle_position_y':self.particle_position_y,
            'particle_position_z':self.particle_position_z,
            'particle_age':self.particle_age,
            'particle_mass':self.particle_mass,
            'particle_mass_initial':self.particle_mass_initial,
            'particle_metallicity':self.particle_metallicity,
            'particle_velocity_x':self.particle_velocity_x,
            'particle_velocity_y':self.particle_velocity_y,
            'particle_velocity_z':self.particle_velocity_z,
            
            #stellar fields
            'star_position_x':self.star_position_x,
            'star_position_y':self.star_position_y,
            'star_position_z':self.star_position_z,
            'star_mass':self.star_mass,
            'star_velocity_x':self.star_velocity_x,
            'star_velocity_y':self.star_velocity_y,
            'star_velocity_z':self.star_velocity_z,
            'star_age':self.star_age,
            'star_metallicity':self.star_metallicity1 + grid.star_metallicity2,
            'star_metallicity1':self.star_metallicity1,
            'star_metallicity2':self.star_metallicity2,
            'star_mass_initial':self.star_mass_initial,
            'star_mass':self.star_mass}
         
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
        self.max_level = pf.max_level
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.float_type = np.float64
        AMRHierarchy.__init__(self,pf,data_style)
        if 'particle_position' in dir(self.pf):
            self._setup_particle_grids()
        self._setup_field_list()
        
    def _setup_particle_grids(self):
        grid_particle_count = na.zeros(len(self.grids),dtype='int64')
        npt = self.pf.particle_position.shape[0]
        if self.pf.do_grid_particles:
            nps = self.pf.star_position.shape[0]
            grid_indices = na.zeros(nps,dtype='int64')
            particle_id= na.arange(nps,dtype='int64')
            pbar = get_pbar("Gridding Particles",len(self.grids))
            grid_indices,grid_particle_count,grids_done = \
                    particle_assignment(self.grids,
                      self.grids[0], 
                      self.pf.star_position,
                      particle_id,
                      grid_indices,
                      grid_particle_count, 
                      self.pf.domain_dimensions,
                      self.pf.max_level,
                      logger=pbar)
            pbar.finish()        
            pbar = get_pbar("Finalizing grids ",len(self.grids))
            for gi, (g,npi) in enumerate(zip(self.grids,grid_particle_count)): 
                star_mask= grid_indices==gi
                if gi==0:
                    #attach all the particles to the root grid
                    g.particle_type = self.pf.particle_type
                    g.particle_id = na.arange(npt)
                    g.particle_mass = self.pf.particle_mass
                    g.particle_mass_initial = self.pf.particle_mass_initial
                    g.particle_age = self.pf.particle_age
                    g.particle_metallicity= self.pf.particle_metallicity
                    g.particle_position_x= self.pf.particle_position[:,0]
                    g.particle_position_y= self.pf.particle_position[:,1]
                    g.particle_position_z= self.pf.particle_position[:,2]
                    g.particle_velocity_x= self.pf.particle_velocity[:,0]
                    g.particle_velocity_y= self.pf.particle_velocity[:,1]
                    g.particle_velocity_z= self.pf.particle_velocity[:,2]
                if star_mask.sum()>0:
                    star_data = self.pf.star_data[star_mask]         
                    (g.star_position_x, \
                        g.star_position_y, \
                        g.star_position_z, \
                        g.star_velocity_x,\
                        g.star_velocity_y,\
                        g.star_velocity_z,\
                        g.star_age,\
                        g.star_metallicity1,\
                        g.star_metallicity2,\
                        g.star_mass_initial,\
                        g.star_mass) = tuple(star_data.T)
                    g.NumberOfParticles = npi        
                self.grids[gi] = g
                pbar.update(gi)
            pbar.finish()
        else:        
            pbar = get_pbar("Finalizing grids ",len(self.grids))
            for gi, g in enumerate(self.grids): 
                if gi==0:
                    #attach all the particles to the root grid
                    g.particle_type = self.pf.particle_type
                    g.particle_id = na.arange(npt)
                    g.particle_mass = self.pf.particle_mass
                    g.particle_mass_initial = self.pf.particle_mass_initial
                    g.particle_age = self.pf.particle_age
                    g.particle_metallicity= self.pf.particle_metallicity
                    g.particle_position_x= self.pf.particle_position[:,0]
                    g.particle_position_y= self.pf.particle_position[:,1]
                    g.particle_position_z= self.pf.particle_position[:,2]
                    g.particle_velocity_x= self.pf.particle_velocity[:,0]
                    g.particle_velocity_y= self.pf.particle_velocity[:,1]
                    g.particle_velocity_z= self.pf.particle_velocity[:,2]
                    if self.pf.do_stars:
                        (g.star_position_x, \
                            g.star_position_y, \
                            g.star_position_z, \
                            g.star_velocity_x,\
                            g.star_velocity_y,\
                            g.star_velocity_z,\
                            g.star_age,\
                            g.star_metallicity1,\
                            g.star_metallicity2,\
                            g.star_mass_initial,\
                            g.star_mass) = tuple(self.pf.star_data.T)
                    g.NumberOfParticles = npt        
                else:
                    g.star_indices = []
                self.grids[gi] = g
            pbar.finish()
            grid_particle_count[0]=npt
        self.grid_particle_count = grid_particle_count

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
        self.pf.level_info = np.array(self.pf.level_info)        
        self.pf.level_offsets = self.pf.level_child_offsets
        self.pf.level_offsets = np.array(self.pf.level_offsets, dtype='int64')
        self.pf.level_offsets[0] = self.pf.root_grid_offset
        
        self.pf.level_art_child_masks = {}
        cm = self.pf.root_iOctCh>0
        cm_shape = (1,)+cm.shape 
        self.pf.level_art_child_masks[0] = cm.reshape(cm_shape).astype('uint8')        
        del cm
        
        root_psg = _ramses_reader.ProtoSubgrid(
                        np.zeros(3, dtype='int64'), # left index of PSG
                        self.pf.domain_dimensions, # dim of PSG
                        np.zeros((1,3), dtype='int64'), # left edges of grids
                        np.zeros((1,6), dtype='int64') # empty
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
                idims = (np.max(dleft_index, axis=0) - initial_left).ravel()+2
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
            eff_mean = np.mean(psg_eff)
            eff_nmin = np.sum([e<=min_eff*tol for e in psg_eff])
            eff_nall = len(psg_eff)
            mylog.info("Average subgrid efficiency %02.1f %%",
                        eff_mean*100.0)
            mylog.info("%02.1f%% (%i/%i) of grids had minimum efficiency",
                        eff_nmin*100.0/eff_nall,eff_nmin,eff_nall)
            
        
            mylog.info("Done with level % 2i; max LE %i", level,na.max(left_index))
            pbar.finish()
            self.proto_grids.append(psgs)
            #print sum(len(psg.grid_file_locations) for psg in psgs)
            #mylog.info("Final grid count: %s", len(self.proto_grids[level]))
            if len(self.proto_grids[level]) == 1: continue
        self.num_grids = sum(len(l) for l in self.proto_grids)
                    
            
            

    num_deep = 0

        
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
                    le = na.zeros(3,dtype='float64')
                    re = na.ones(3,dtype='float64')
                    gd = dd
                self.grid_left_edge[gi,:] = le
                self.grid_right_edge[gi,:] = re
                self.grid_dimensions[gi,:] = gd
                assert na.all(self.grid_left_edge[gi,:]<=1.0)    
                self.grid_levels[gi,:] = level
                child_mask = np.zeros(props[2,:],'uint8')
                amr_utils.fill_child_mask(fl,start_index,
                    self.pf.level_art_child_masks[level],
                    child_mask)
                grids.append(self.grid(gi, self, level, fl, 
                    start_index,le,re,gd))
                gi += 1
        self.grids = np.empty(len(grids), dtype='object')
        

        if self.pf.file_particle_data:
            lspecies = self.pf.parameters['lspecies']
            wspecies = self.pf.parameters['wspecies']
            Nrow     = self.pf.parameters['Nrow']
            nstars = na.diff(lspecies)[-1]
            a = self.pf.parameters['aexpn']
            hubble = self.pf.parameters['hubble']
            ud  = self.pf.parameters['r0']*a/hubble #proper Mpc units
            uv  = self.pf.parameters['v0']/(a**1.0)*1.0e5 #proper cm/s
            um  = self.pf.parameters['aM0'] #mass units in solar masses
            um *= 1.989e33 #convert solar masses to grams 
            pbar = get_pbar("Loading Particles   ",5)
            self.pf.particle_position,self.pf.particle_velocity = \
                read_particles(self.pf.file_particle_data,Nrow)
            pbar.update(1)
            npa,npb=0,0
            npb = lspecies[-1]
            clspecies = np.concatenate(([0,],lspecies))
            if self.pf.only_particle_type is not None:
                npb = lspecies[0]
                if type(self.pf.only_particle_type)==type(5):
                    npa = clspecies[self.pf.only_particle_type]
                    npb = clspecies[self.pf.only_particle_type+1]
            np = npb-npa
            nparticles = np
            #make sure we aren't going to throw out good particles
            if not na.all(self.pf.particle_position[npb:]==0.0):
                print 'WARNING: unused particles discovered from lspecies'
            self.pf.particle_position   = self.pf.particle_position[npa:npb]
            #do NOT correct by an offset of 1.0
            #self.pf.particle_position  -= 1.0 #fortran indices start with 0
            pbar.update(2)
            self.pf.particle_position  /= self.pf.domain_dimensions 
            #to unitary units (comoving)
            pbar.update(3)
            self.pf.particle_velocity   = self.pf.particle_velocity[npa:npb]
            self.pf.particle_velocity  *= uv #to proper cm/s
            pbar.update(4)
            self.pf.particle_type         = np.zeros(np,dtype='uint8')
            self.pf.particle_mass         = np.zeros(np,dtype='float64')
            self.pf.particle_mass_initial = np.zeros(np,dtype='float64')-1
            self.pf.particle_creation_time= np.zeros(np,dtype='float64')-1
            self.pf.particle_metallicity  = np.zeros(np,dtype='float64')-1
            self.pf.particle_metallicity1 = np.zeros(np,dtype='float64')-1
            self.pf.particle_metallicity2 = np.zeros(np,dtype='float64')-1
            self.pf.particle_age          = np.zeros(np,dtype='float64')-1

            dist = self.pf['cm']/self.pf.domain_dimensions[0]
            self.pf.conversion_factors['particle_mass'] = 1.0 #solar mass in g
            self.pf.conversion_factors['particle_mass_initial'] = 1.0
            self.pf.conversion_factors['particle_species'] = 1.0
            for ax in 'xyz':
                self.pf.conversion_factors['particle_velocity_%s'%ax] = 1.0
                #already in unitary units
                self.pf.conversion_factors['particle_position_%s'%ax] = 1.0 
            self.pf.conversion_factors['particle_creation_time'] =  31556926.0
            self.pf.conversion_factors['particle_metallicity']=1.0
            self.pf.conversion_factors['particle_metallicity1']=1.0
            self.pf.conversion_factors['particle_metallicity2']=1.0
            self.pf.conversion_factors['particle_index']=1.0
            self.pf.conversion_factors['particle_type']=1
            self.pf.conversion_factors['particle_age']=1
            self.pf.conversion_factors['Msun'] = 5.027e-34 
            #conversion to solar mass units
            

            a,b=0,0
            self.pf.particle_star_index = len(wspecies)-1
            for i,(b,m) in enumerate(zip(lspecies,wspecies)):
                if type(self.pf.only_particle_type)==type(5):
                    if not i==self.pf.only_particle_type:
                        continue
                    self.pf.particle_type += i
                    self.pf.particle_mass += m*um

                else:
                    self.pf.particle_type[a:b] = i #particle type
                    self.pf.particle_mass[a:b] = m*um #mass in solar masses
                    if m==0.0:
                        self.pf.particle_star_index = i
                a=b
            pbar.finish()

            lparticles = [0,]+list(lspecies)
            for j,np in enumerate(lparticles):
                mylog.debug('found %i of particle type %i'%(j,np))
            
            
            do_stars = (self.pf.only_particle_type is None) or \
                       (self.pf.only_particle_type == -1) or \
                       (self.pf.only_particle_type == len(lspecies))
            self.pf.do_stars = do_stars           
            if self.pf.file_star_data and do_stars: 
                nstars_pa = nstars
                (nstars_rs,), mass, imass, tbirth, metallicity1, metallicity2, \
                        ws_old,ws_oldi,tdum,adum \
                     = read_stars(self.pf.file_star_data)
                self.pf.nstars_rs = nstars_rs     
                self.pf.nstars_pa = nstars_pa
                if not nstars_rs==na.sum(self.pf.particle_type==self.pf.particle_star_index):
                    print 'WARNING!: nstars is inconsistent!'
                if nstars_rs > 0 :
                    n=min(1e2,len(tbirth))
                    pbar = get_pbar("Stellar Ages        ",n)
                    birthtimes= \
                        b2t(tbirth,n=n,logger=lambda x: pbar.update(x)).astype('float64')
                    assert birthtimes.shape == tbirth.shape    
                    birthtimes*= 1.0e9 #from Gyr to yr
                    birthtimes*= 365*24*3600 #to seconds
                    ages = self.pf.current_time-birthtimes
                    spread = self.pf.spread
                    if spread == False:
                        pass
                    elif type(spread)==type(5.5):
                        ages = spread_ages(ages,spread=spread)
                    else:
                        ages = spread_ages(ages)
                    idx = self.pf.particle_type == self.pf.particle_star_index    
                    assert na.sum(idx)==nstars_pa
                    self.pf.star_position = self.pf.particle_position[idx]
                    self.pf.star_velocity = self.pf.particle_velocity[idx]
                    self.pf.particle_age[idx] = ages
                    self.pf.star_age = ages
                    pbar.finish()
                    self.pf.particle_metallicity[idx] = metallicity1+metallicity2
                    self.pf.particle_metallicity1[idx] = metallicity1
                    self.pf.particle_metallicity2[idx] = metallicity2
                    self.pf.particle_mass[idx] = mass*um
                    self.pf.particle_mass_initial[idx] = mass*um
                    self.pf.star_metallicity1 = metallicity1
                    self.pf.star_metallicity2 = metallicity2
                    self.pf.star_mass_initial = imass*um
                    self.pf.star_mass = mass*um
                    self.pf.star_data = na.array([
                        self.pf.star_position[:,0],
                        self.pf.star_position[:,1],
                        self.pf.star_position[:,2],
                        self.pf.star_velocity[:,0],
                        self.pf.star_velocity[:,1],
                        self.pf.star_velocity[:,2],
                        self.pf.star_age,
                        self.pf.star_metallicity1,
                        self.pf.star_metallicity2,
                        self.pf.star_mass_initial,
                        self.pf.star_mass]).T

            done = 0
            init = self.pf.particle_position.shape[0]
            pos = self.pf.particle_position
            #particle indices travel with the particle positions
            #pos = na.vstack((na.arange(pos.shape[0]),pos.T)).T 
        for gi,g in enumerate(grids):    
            self.grids[gi]=g
                    
    def _get_grid_parents(self, grid, LE, RE):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(LE, RE)
        mask[grid_ind] = True
        mask = np.logical_and(mask, (self.grid_levels == (grid.Level-1)).flat)
        return self.grids[mask]

    def _populate_grid_objects(self):
        mask = na.empty(self.grids.size, dtype='int32')
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
        pb.finish()
        self.max_level = self.grid_levels.max()

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
                 discover_particles=True,
                 limit_level=None,
                 only_particle_type = None,
                 do_grid_particles=False,
                 merge_dm_and_stars=False,
                 spread = True,
                 single_particle_mass=False,
                 single_particle_type=0):
        
        #dirn = os.path.dirname(filename)
        base = os.path.basename(filename)
        aexp = base.split('_')[2].replace('.d','')
        if not aexp.startswith('a'):
            aexp = '_'+aexp
        
        self.file_particle_header = file_particle_header
        self.file_particle_data = file_particle_data
        self.file_star_data = file_star_data
        self.only_particle_type = only_particle_type
        self.do_grid_particles = do_grid_particles
        self.single_particle_mass = single_particle_mass
        self.merge_dm_and_stars = merge_dm_and_stars
        self.spread = spread
        
        if limit_level is None:
            self.limit_level = np.inf
        else:
            limit_level = int(limit_level)
            mylog.info("Using maximum level: %i",limit_level)
            self.limit_level = limit_level
        
        def repu(x):
            for i in range(5):
                x=x.replace('__','_')
            return x    
        if discover_particles:
            if file_particle_header is None:
                loc = filename.replace(base,'PMcrd%s.DAT'%aexp)
                loc = repu(loc)
                if os.path.exists(loc):
                    self.file_particle_header = loc
                    mylog.info("Discovered particle header: %s",os.path.basename(loc))
            if file_particle_data is None:
                loc = filename.replace(base,'PMcrs0%s.DAT'%aexp)
                loc = repu(loc)
                if os.path.exists(loc):
                    self.file_particle_data = loc
                    mylog.info("Discovered particle data:   %s",os.path.basename(loc))
            if file_star_data is None:
                loc = filename.replace(base,'stars_%s.dat'%aexp)
                loc = repu(loc)
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
        #t0 = 6.17e17/(self.hubble_constant + np.sqrt(self.omega_matter))
        #v0 = r0 / t0
        #rho0 = 1.8791e-29 * self.hubble_constant**2.0 * self.omega_matter
        #e0 = v0**2.0
        
        wmu = self["wmu"]
        boxh = self["boxh"]
        aexpn = self["aexpn"]
        hubble = self.hubble_constant
        ng = self.domain_dimensions[0]
        self.r0 = boxh/ng
        self.v0 =  self.r0 * 50.0*1.0e5 * np.sqrt(self.omega_matter)  #cm/s
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
        #self.conversion_factors["Temperature"] = tr 
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
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype="float64")
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
        #     return 1./(x*np.sqrt(oml+omb*x**-3.0))
        # spacings = np.logspace(-5,np.log10(self.parameters['aexpn']),1e5)
        # integrand_arr = integrand(spacings)
        # self.current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
        # self.current_time *= self.hubble_time
        self.current_time = b2t(self.current_time_raw)*1.0e9*365*3600*24         
        for to_skip in ['tl','dtl','tlold','dtlold','iSO']:
            _skip_record(f)

        
        Om0 = self.parameters['Om0']
        hubble = self.parameters['hubble']
        dummy = 100.0 * hubble * np.sqrt(Om0)
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
                np.sqrt(self.parameters["Om0"])
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
        est = int(np.rint(self.ncell**(1.0/3.0)))
        # Note here: this is the number of *cells* on the root grid.
        # This is not the same as the number of Octs.
        self.domain_dimensions = np.ones(3, dtype='int64')*est 

        self.root_grid_mask_offset = f.tell()
        #_skip_record(f) # iOctCh
        root_cells = self.domain_dimensions.prod()
        self.root_iOctChfull = _read_frecord(f,'>i')
        self.root_iOctCh = self.root_iOctChfull[:root_cells]
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
        self.parameters['wspeciesf'] = np.fromfile(fh,dtype='>f',count=10)
        self.parameters['lspeciesf'] = np.fromfile(fh,dtype='>i',count=10)
        assert np.all(self.parameters['lspeciesf'][n:]==0.0)
        assert np.all(self.parameters['wspeciesf'][n:]==0.0)
        self.parameters['wspecies'] = self.parameters['wspeciesf'][:n]
        self.parameters['lspecies'] = self.parameters['lspeciesf'][:n]
        fh.close()
        
        ls_nonzero = [ls for ls in self.parameters['lspecies'] if ls>0 ]
        mylog.debug("Discovered %i species of particles",len(ls_nonzero))
        mylog.debug("Particle populations: "+'%1.1e '*len(ls_nonzero),ls_nonzero)
        

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if "10MpcBox" in args[0]:
            return True
        return False

def particle_assignment(grids,this_grid, 
                                  pos,
                                  particle_id,
                                  grid_indices,
                                  grid_particle_count, 
                                  domain_dimensions,
                                  max_level,
                                  subdiv=2,
                                  grids_done=0,
                                  logger=None):
    #for every particle check every child grid to see if it fits inside
    #cast the pos -> cell location index (instead of doing a LE<pos<RE check)
    #find if cell descends into the next mesh
    
    #cast every position into a cell on this grid
    #we may get negative indices or indices outside this grid
    #mask them out
    exp = domain_dimensions*subdiv**this_grid.Level
    lei= na.floor((pos-this_grid.LeftEdge)*exp).astype('int64')

    #now lookup these indices in the child index mask
    #throw out child grids = -1 and particles outside the range
    #default state is to not grid a particle
    child_idx = na.zeros(lei.shape[0],dtype='int64')-1
    #remove particles to the left or right of the grid
    lei_out  = na.any(lei>=this_grid.ActiveDimensions,axis=1)
    lei_out |= na.any(lei<0,axis=1)
    #lookup grids for every particle except the ones to the 
    leio=lei[~lei_out]
    #child_idx[~lei_out]= \
    child_idx[~lei_out]= \
            this_grid.child_index_mask[(leio[:,0],leio[:,1],leio[:,2])]
    mask = (child_idx > -1)
    #only assign the particles if they point to a grid ID that isnt -1
    grid_indices[particle_id[mask]] = child_idx[mask]
    #the number of particles on this grid is equal to those
    #that point to -1
    grid_particle_count[this_grid.id] = na.sum(~mask)
    grids_done +=1
    if logger:
        logger.update(grids_done)

    for child_grid_index in na.unique(this_grid.child_index_mask):
        if child_grid_index == -1: 
            continue
        if grids[child_grid_index].Level == max_level:
            continue
        mask = child_idx == child_grid_index
        if na.sum(mask)==0:continue
        grid_indices,grid_particle_count,grids_done = \
        particle_assignment(grids,grids[child_grid_index],
                pos[mask],particle_id[mask],
                grid_indices,grid_particle_count,
                domain_dimensions,max_level,grids_done=grids_done,
                subdiv=subdiv,logger=logger)
    return grid_indices,grid_particle_count,grids_done

