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

def num_deep_inc(f):
    def wrap(self, *args, **kwargs):
        self.num_deep += 1
        rv = f(self, *args, **kwargs)
        self.num_deep -= 1
        return rv
    return wrap

class ARTGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, hierarchy, level, locations, start_index):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Level = level
        self.Parent = []
        self.Children = []
        self.locations = locations
        self.start_index = start_index.copy()

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

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        # This will need to be generalized to be used elsewhere.
        self.field_list = [ 'Density','TotalEnergy',
                            'x-momentum','y-momentum','z-momentum',
                            'Pressure','Gamma','GasEnergy',
                            'Metal_DensitySNII', 'Metal_DensitySNIa',
                            'Potential_New','Potential_Old']
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        LEVEL_OF_EDGE = 7
        MAX_EDGE = (2 << (LEVEL_OF_EDGE- 1))
        
        min_eff = 0.40
        
        f = open(self.pf.parameter_filename,'rb')
        self.pf.nhydro_vars, self.pf.level_info, self.pf.level_offsetsa = \
                         _count_art_octs(f, 
                          self.pf.child_grid_offset,
                          self.pf.min_level, self.pf.max_level)
        self.pf.level_info[0]=self.pf.ncell
        f.close()
        self.pf.level_info = na.array(self.pf.level_info)
        num_ogrids = sum(self.pf.level_info) + self.pf.iOctFree
        print 'found %i oct grids'%num_ogrids
        num_ogrids *=7
        print 'instantiating... %i grids'%num_ogrids
        ogrid_left_indices = na.zeros((num_ogrids,3), dtype='int64') - 999
        ogrid_levels = na.zeros(num_ogrids, dtype='int64')
        ogrid_file_locations = na.zeros((num_ogrids,6), dtype='int64')
        
        #don't need parents?
        #ogrid_parents = na.zeros(num_ogrids, dtype="int64")
        
        #don't need masks?
        #ochild_masks = na.zeros((num_ogrids, 8), dtype='int64').ravel()
        
        self.pf.level_offsets = amr_utils.read_art_tree(
                                self.pf.parameter_filename, 
                                self.pf.child_grid_offset,
                                self.pf.min_level, self.pf.max_level,
                                ogrid_left_indices, ogrid_levels,
                                ogrid_file_locations)
                                #ochild_masks,
                                #ogrid_parents, 
                                
        self.pf.level_offsets = na.array(self.pf.level_offsets, dtype='int64')
        self.pf.level_offsets[0] = self.pf.root_grid_offset
        #ochild_masks.reshape((num_ogrids, 8), order="F")
        ogrid_levels[ogrid_left_indices[:,0] == -999] = -1
        # This bit of code comes from Chris, and I'm still not sure I have a
        # handle on what it does.
        final_indices =  ogrid_left_indices[na.where(ogrid_levels==self.pf.max_level)[0]]
        divisible=[na.all((final_indices%2**(level))==0) 
            for level in xrange(self.pf.max_level*2)]
        root_level = self.pf.max_level+na.where(na.logical_not(divisible))[0][0] 
        ogrid_dimension = na.zeros(final_indices.shape,dtype='int')+2
        ogrid_left_indices = ogrid_left_indices/2**(root_level - ogrid_levels[:,None] - 1) - 1

        # Now we can rescale
        # root_psg = _ramses_reader.ProtoSubgrid(
        #                 na.zeros(3, dtype='int64'), # left index of PSG
        #                 self.pf.domain_dimensions, # dim of PSG
        #                 na.zeros((1,3), dtype='int64'), # left edges of grids
        #                 self.pf.domain_dimensions[None,:], # right edges of grids
        #                 self.pf.domain_dimensions[None,:], # dims of grids
        #                 na.zeros((1,6), dtype='int64') # empty
        #                 )
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
            ggi = (ogrid_levels == level).ravel()
            nd = self.pf.domain_dimensions * 2**level
            dims = na.ones((ggi.sum(), 3), dtype='int64') * 2
            fl = ogrid_file_locations[ggi,:]
            psgs = []
            
            #refers to the left index for the art octgrid
            left_index = ogrid_left_indices[ggi,:]
            
            # We now calculate the hilbert curve position of every left_index,
            # of the octs, with respect to a lower order hilbert curve.
            left_index_gridpatch = left_index >> LEVEL_OF_EDGE
            order = max(level + 1 - LEVEL_OF_EDGE, 0)
            
            #compute the hilbert indices up to a certain level
            #the indices will associate an oct grid to the nearest
            #hilbert index?
            hilbert_indices = _ramses_reader.get_hilbert_indices(order, left_index_gridpatch)
            
            # Strictly speaking, we don't care about the index of any
            # individual oct at this point.  So we can then split them up.
            unique_indices = na.unique(hilbert_indices)
            mylog.debug("Level % 2i has % 10i unique indices for %0.3e octs",
                        level, unique_indices.size, hilbert_indices.size)
            
            #use the hilbert indices to order oct grids so that consecutive
            #items on a list are spatially near each other
            #this is useful because we will define grid patches over these
            #octs, which are more efficient if the octs are spatially close
            
            #split into list of lists, with domains containing 
            #lists of sub octgrid left indices and an index
            #referring to the domain on which they live
            locs, lefts = _ramses_reader.get_array_indices_lists(
                        hilbert_indices, unique_indices, left_index, fl)
            
            #iterate over the domains    
            step=0
            pbar = get_pbar("Re-gridding  Level %i"%level, len(locs))
            psg_eff = []
            psg_dep = []
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
                self.num_deep = 0
                
                #because grid patches may still be mostly empty, and with octs
                #that only partially fill the grid,it  may be more efficient
                #to split large patches into smaller patches. We split
                #if less than 10% the volume of a patch is covered with octs
                psg_split = _ramses_reader.recursive_patch_splitting(
                    psg, idims, initial_left, 
                    dleft_index, dfl,min_eff=min_eff)
                    
                psgs.extend(psg_split)
                
                tol = 1.00001
                psg_eff  += [x.efficiency for x in psg_split] 
                psg_dep  += [x.num_deep for x in psg_split] 
                
                step+=1
                pbar.update(step)
            eff_mean = na.mean(psg_eff)
            eff_nmin = na.sum([e<=min_eff*tol for e in psg_eff])
            eff_nall = len(psg_eff)
            dep_mean = na.rint(na.mean(psg_dep))
            mylog.info("Average subgrid efficiency %02.1f %% and average depth %i",
                        eff_mean*100.0, dep_mean)
            mylog.info("%02.1f%% (%i/%i) of grids had minimum efficiency",
                        eff_nmin*100.0/eff_nall,eff_nmin,eff_nall)
            mylog.info("Re-gridding level %i: %s octree grids", level, ggi.sum())
            
        
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
            correction=1.0
            if level == 0:
                correction=64.0
            for g in grid_list:
                fl = g.grid_file_locations
                props = g.get_properties()
                dds = ((2**level) * self.pf.domain_dimensions).astype("float64")
                self.grid_left_edge[gi,:] = props[0,:]*correction / dds
                self.grid_right_edge[gi,:] = props[1,:]*correction / dds
                self.grid_dimensions[gi,:] = props[2,:]*correction
                self.grid_levels[gi,:] = level
                grids.append(self.grid(gi, self, level, fl, props[0,:]))
                gi += 1
        self.grids = na.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g

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
                 storage_filename = None):
        import yt.frontends.ramses._ramses_reader as _ramses_reader
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename
        
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'art'
        self.parameters["Time"] = 1. # default unit is 1...
        self.parameters["InitialTime"]=self.current_time
        
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
        self.conversion_factors["Density"] = \
            self.rho0*(aexpn**-3.0)
        self.conversion_factors["GasEnergy"] = \
            self.rho0*self.v0**2*(aexpn**-5.0)
        tr  = self.tr
        self.conversion_factors["Temperature"] = tr
        self.conversion_factors["Metal_Density"] = 1
        self.cosmological_simulation = True
        
        # Now our conversion factors
        for ax in 'xyz':
            # Add on the 1e5 to get to cm/s
            self.conversion_factors["%s-velocity" % ax] = self.v0/aexpn
        seconds = self.t0
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)

        #we were already in seconds, go back in to code units
        self.current_time /= self.t0 
        
        
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
        self.current_redshift = self.parameters["aexpn"]**-1.0 - 1.0
        self.data_comment = header_vals['jname']
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
        def integrand(x,oml=self.omega_lambda,omb=self.omega_matter):
            return 1./(x*na.sqrt(oml+omb*x**-3.0))
        spacings = na.logspace(-5,na.log10(self.parameters['aexpn']),1e5)
        integrand_arr = integrand(spacings)
        self.current_time = na.trapz(integrand_arr,dx=na.diff(spacings))
        self.current_time *= self.hubble_time
                
        for to_skip in ['tl','dtl','tlold','dtlold','iSO']:
            _skip_record(f)

        (self.ncell,) = struct.unpack('>l', _read_record(f))
        # Try to figure out the root grid dimensions
        est = int(na.rint(self.ncell**(1.0/3.0)))
        # Note here: this is the number of *cells* on the root grid.
        # This is not the same as the number of Octs.
        self.domain_dimensions = na.ones(3, dtype='int64')*est 

        self.root_grid_mask_offset = f.tell()
        _skip_record(f) # iOctCh
        self.root_grid_offset = f.tell()
        _skip_record(f) # hvar
        _skip_record(f) # var

        self.iOctFree, self.nOct = struct.unpack('>ii', _read_record(f))
        self.child_grid_offset = f.tell()

        f.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        return False # We make no effort to auto-detect ART data

def _skip_record(f):
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    f.seek(s[0], 1)
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))

def _read_record(f):
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))[0]
    ss = f.read(s)
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    return ss

def _read_record_size(f):
    pos = f.tell()
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    f.seek(pos)
    return s[0]

def _count_art_octs(f, offset, 
                   MinLev, MaxLevelNow):
    level_offsets= []
    f.seek(offset)
    nchild,ntot=8,0
    Level = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iNOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iHOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    for Lev in xrange(MinLev + 1, MaxLevelNow+1):
        level_offsets.append(f.tell())
        
        #Get the info for this level, skip the rest
        #print "Reading oct tree data for level", Lev
        #print 'offset:',f.tell()
        Level[Lev], iNOLL[Lev], iHOLL[Lev] = struct.unpack(
           '>iii', _read_record(f))
        #print 'Level %i : '%Lev, iNOLL
        #print 'offset after level record:',f.tell()
        iOct = iHOLL[Lev] - 1
        nLevel = iNOLL[Lev]
        nLevCells = nLevel * nchild
        ntot = ntot + nLevel

        #Skip all the oct hierarchy data
        ns = _read_record_size(f)
        size = struct.calcsize('>i') + ns + struct.calcsize('>i')
        f.seek(f.tell()+size * nLevel)
        
        #Skip the child vars data
        ns = _read_record_size(f)
        size = struct.calcsize('>i') + ns + struct.calcsize('>i')
        f.seek(f.tell()+size * nLevel*nchild)
        
        #find nhydrovars
        nhydrovars = 8+2
    f.seek(offset)
    return nhydrovars, iNOLL, level_offsets

def _read_art_level(f, level_offsets,level):
    pos = f.tell()
    f.seek(level_offsets[leve])
    #Get the info for this level, skip the rest
    #print "Reading oct tree data for level", Lev
    #print 'offset:',f.tell()
    Level[Lev], iNOLL[Lev], iHOLL[Lev] = struct.unpack(
       '>iii', _read_record(f))
    #print 'Level %i : '%Lev, iNOLL
    #print 'offset after level record:',f.tell()
    iOct = iHOLL[Lev] - 1
    nLevel = iNOLL[Lev]
    nLevCells = nLevel * nchild
    ntot = ntot + nLevel

    #Skip all the oct hierarchy data
    #in the future, break this up into large chunks
    count = nLevel*15
    le  = numpy.zeros((count,3),dtype='int64')
    fl  = numpy.zeros((count,6),dtype='int64')
    idxa,idxb = 0,0
    chunk = 1e9 #this is ~111MB for 15 dimensional 64 bit arrays
    while left > 0 :
        data = na.fromfile(f,dtype='>i',count=chunk*15)
        data.reshape(chunk,15)
        left = count-index
        le[idxa:idxb,:] = data[0:3]
        fl[idxa:idxb,1] = numpy.arange(chunk)
    del data
    f.seek(pos)
    return le,fl

