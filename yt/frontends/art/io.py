"""
ART-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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
import struct

import os
import os.path

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.lib as au

from yt.frontends.art.definitions import art_particle_field_names

class IOHandlerART(BaseIOHandler):
    _data_style = "art"

    def __init__(self, filename, nhydro_vars, level_info, level_offsets,
                 *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.filename = filename
        self.nhydro_vars = nhydro_vars
        self.level_info = level_info
        self.level_offsets = level_offsets
        self.level_data = {}

    def preload_level(self, level,field=None):
        """ Reads in the full ART tree. From the ART source:
            iOctLv :    >0   - level of an oct
            iOctPr :         - parent of an oct
            iOctCh :    >0   - pointer to an oct of children
                        0   - there are no children; the cell is a leaf
            iOctNb :    >0   - pointers to neighbouring cells 
            iOctPs :         - coordinates of Oct centers
            
            iOctLL1:         - doubly linked list of octs
            iOctLL2:         - doubly linked list of octs
            
            tl - current  time moment for level L
            tlold - previous time moment for level L
            dtl - dtime0/2**iTimeBin
            dtlold -  previous time step for level L
            iSO - sweep order
            
            hvar(1,*) - gas density 
            hvar(2,*) - gas energy 
            hvar(3,*) - x-momentum 
            hvar(4,*) - y-momentum
            hvar(5,*) - z-momentum
            hvar(6,*) - pressure
            hvar(7,*) - Gamma
            hvar(8,*) - internal energy 

            var (1,*) - total density 
            var (2,*) - potential (new)
            var (3,*) - potential (old)
            
            
            
        """
        
        if level in self.level_data: return
        if level == 0:
            self.preload_root_level()
            return
        f = open(self.filename, 'rb')
        f.seek(self.level_offsets[level])
        ncells = 8*self.level_info[level]
        nvals = ncells * (self.nhydro_vars + 6) # 2 vars, 2 pads
        arr = na.fromfile(f, dtype='>f', count=nvals)
        arr = arr.reshape((self.nhydro_vars+6, ncells), order="F")
        assert na.all(arr[0,:]==arr[-1,:]) #pads must be equal
        arr = arr[3:-1,:] #skip beginning pad, idc, iOctCh, + ending pad
        if field==None:
            self.level_data[level] = arr.astype('float32')
        else:
            self.level_data[level] = arr.astype('float32')
        del arr

    def preload_root_level(self):
        f = open(self.filename, 'rb')
        f.seek(self.level_offsets[0] + 4) # Ditch the header
        ncells = self.level_info[0]
        nhvals = ncells * (self.nhydro_vars) # 0 vars, 0 pads
        hvar = na.fromfile(f, dtype='>f', count=nhvals).astype("float32")
        hvar = hvar.reshape((self.nhydro_vars, ncells), order="F")
        na.fromfile(f,dtype='>i',count=2) #throw away the pads
        nvars = ncells * (2) # 0 vars, 0 pads
        var = na.fromfile(f, dtype='>f', count=nvars).astype("float32")
        var = var.reshape((2, ncells), order="F")
        arr = na.concatenate((hvar,var))
        self.level_data[0] = arr

    def clear_level(self, level):
        self.level_data.pop(level, None)

    def _read_particle_field(self, grid, field):
        #This will be cleaned up later
        import pdb; pdb.set_trace()
        if field == 'particle_index':
            return grid.particle_id
        if field == 'particle_type':
            return grid.particle_type
        if field == 'particle_position_x':
            return grid.particle_position_x
        if field == 'particle_position_y':
            return grid.particle_position_y
        if field == 'particle_position_z':
            return grid.particle_position_z
        if field == 'particle_age':
            return grid.particle_age
        if field == 'particle_mass':
            return grid.particle_mass
        if field == 'particle_mass_initial':
            return grid.particle_mass_initial
        if field == 'particle_metallicity':
            return grid.particle_metallicity
        if field == 'particle_velocity_x':
            return grid.particle_velocity_x
        if field == 'particle_velocity_y':
            return grid.particle_velocity_y
        if field == 'particle_velocity_z':
            return grid.particle_velocity_z
        
        #stellar fields
        if field == 'star_position_x':
            return grid.star_position_x
        if field == 'star_position_y':
            return grid.star_position_y
        if field == 'star_position_z':
            return grid.star_position_z
        if field == 'star_mass':
            return grid.star_mass
        if field == 'star_velocity_x':
            return grid.star_velocity_x
        if field == 'star_velocity_y':
            return grid.star_velocity_y
        if field == 'star_velocity_z':
            return grid.star_velocity_z
        if field == 'star_age':
            return grid.star_age
        if field == 'star_metallicity':
            return grid.star_metallicity1 +\
                   grid.star_metallicity2
        if field == 'star_metallicity1':
            return grid.star_metallicity1
        if field == 'star_metallicity2':
            return grid.star_metallicity2
        if field == 'star_mass_initial':
            return grid.star_mass_initial
        if field == 'star_mass':
            return grid.star_mass
        
        raise 'Should have matched one of the particle fields...'

        
    def _read_data_set(self, grid, field):
        if field in art_particle_field_names:
            return self._read_particle_field(grid, field)
        pf = grid.pf
        field_id = grid.pf.h.field_list.index(field)
        if grid.Level == 0: # We only have one root grid
            self.preload_level(0)
            tr = self.level_data[0][field_id,:].reshape(
                    pf.domain_dimensions, order="F").copy()
            return tr.swapaxes(0, 2).astype("float64")
        tr = na.zeros(grid.ActiveDimensions, dtype='float32')
        grids = [grid]
        l_delta = 0
        filled = na.zeros(grid.ActiveDimensions, dtype='uint8')
        to_fill = grid.ActiveDimensions.prod()
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                self.preload_level(g.Level,field=field_id)
                #print "Filling %s from %s (%s)" % (grid, g, g.Level)
                to_fill -= au.read_art_grid(field_id, 
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, self.level_data[g.Level],
                        g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        return self._read_data_set(grid, field)[sl]

def _count_art_octs(f, offset, 
                   MinLev, MaxLevelNow):
    level_oct_offsets= [0,]
    level_child_offsets= [0,]
    f.seek(offset)
    nchild,ntot=8,0
    Level = na.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    iNOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    iHOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    for Lev in xrange(MinLev + 1, MaxLevelNow+1):
        level_oct_offsets.append(f.tell())

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

        level_child_offsets.append(f.tell())
        #Skip the child vars data
        ns = _read_record_size(f)
        size = struct.calcsize('>i') + ns + struct.calcsize('>i')
        f.seek(f.tell()+size * nLevel*nchild)

        #find nhydrovars
        nhydrovars = 8+2
    f.seek(offset)
    return nhydrovars, iNOLL, level_oct_offsets, level_child_offsets

def _read_art_level_info(f, level_oct_offsets,level,coarse_grid=128):
    pos = f.tell()
    f.seek(level_oct_offsets[level])
    #Get the info for this level, skip the rest
    junk, nLevel, iOct = struct.unpack(
       '>iii', _read_record(f))
    
    #fortran indices start at 1
    
    #Skip all the oct hierarchy data
    le     = na.zeros((nLevel,3),dtype='int64')
    fl     = na.ones((nLevel,6),dtype='int64')
    iocts  = na.zeros(nLevel+1,dtype='int64')
    idxa,idxb = 0,0
    chunk = long(1e6) #this is ~111MB for 15 dimensional 64 bit arrays
    left = nLevel
    while left > 0 :
        this_chunk = min(chunk,left)
        idxb=idxa+this_chunk
        data = na.fromfile(f,dtype='>i',count=this_chunk*15)
        data=data.reshape(this_chunk,15)
        left-=this_chunk
        le[idxa:idxb,:] = data[:,1:4]
        fl[idxa:idxb,1] = na.arange(idxa,idxb)
        #pad byte is last, LL2, then ioct right before it
        iocts[idxa:idxb] = data[:,-3] 
        idxa=idxa+this_chunk
    del data
    
    #ioct always represents the index of the next variable
    #not the current, so shift forward one index
    #the last index isn't used
    ioctso = iocts.copy()
    iocts[1:]=iocts[:-1] #shift
    iocts = iocts[:nLevel] #chop off the last index
    iocts[0]=iOct #starting value

    #now correct iocts for fortran indices start @ 1
    iocts = iocts-1

    assert na.unique(iocts).shape[0] == nLevel
    
    #ioct tries to access arrays much larger than le & fl
    #just make sure they appear in the right order, skipping
    #the empty space in between
    idx = na.argsort(iocts)
    
    #now rearrange le & fl in order of the ioct
    le = le[idx]
    fl = fl[idx]


    #left edges are expressed as if they were on 
    #level 15, so no matter what level max(le)=2**15 
    #correct to the yt convention
    #le = le/2**(root_level-1-level)-1

    #try to find the root_level first
    root_level=na.floor(na.log2(le.max()*1.0/coarse_grid))
    root_level = root_level.astype('int64')

    #try without the -1
    le = le/2**(root_level+1-level)-1

    #now read the hvars and vars arrays
    #we are looking for iOctCh
    #we record if iOctCh is >0, in which it is subdivided
    iOctCh  = na.zeros((nLevel+1,8),dtype='bool')
    
    
    
    f.seek(pos)
    return le,fl,nLevel,root_level


def read_particles(file,Nrow):
    words = 6 # words (reals) per particle: x,y,z,vx,vy,vz
    real_size = 4 # for file_particle_data; not always true?
    np_per_page = Nrow**2 # defined in ART a_setup.h
    num_pages = os.path.getsize(file)/(real_size*words*np_per_page)

    f = na.fromfile(file, dtype='>f4').astype('float32') # direct access
    pages = na.vsplit(na.reshape(f, (num_pages, words, np_per_page)), num_pages)
    data = na.squeeze(na.dstack(pages)).T # x,y,z,vx,vy,vz
    return data[:,0:3],data[:,3:]

def read_stars(file):
    fh = open(file,'rb')
    tdum,adum   = _read_frecord(fh,'>d')
    nstars      = _read_frecord(fh,'>i')
    ws_old, ws_oldi = _read_frecord(fh,'>d')
    mass    = _read_frecord(fh,'>f') 
    imass   = _read_frecord(fh,'>f') 
    tbirth  = _read_frecord(fh,'>f') 
    if fh.tell() < os.path.getsize(file):
        metallicity1 = _read_frecord(fh,'>f') 
    if fh.tell() < os.path.getsize(file):
        metallicity2 = _read_frecord(fh,'>f')     
    assert fh.tell() == os.path.getsize(file)
    return  nstars, mass, imass, tbirth, metallicity1, metallicity2,\
            ws_old,ws_oldi,tdum,adum

def _read_child_mask_level(f, level_child_offsets,level,nLevel,nhydro_vars):
    f.seek(level_child_offsets[level])
    nvals = nLevel * (nhydro_vars + 6) # 2 vars, 2 pads
    ioctch = na.zeros(nLevel,dtype='uint8')
    idc = na.zeros(nLevel,dtype='int32')
    
    chunk = long(1e6)
    left = nLevel
    width = nhydro_vars+6
    a,b=0,0
    while left > 0:
        chunk = min(chunk,left)
        b += chunk
        arr = na.fromfile(f, dtype='>i', count=chunk*width)
        arr = arr.reshape((width, chunk), order="F")
        assert na.all(arr[0,:]==arr[-1,:]) #pads must be equal
        idc[a:b]    = arr[1,:]-1 #fix fortran indexing
        ioctch[a:b] = arr[2,:]==0 #if it is above zero, then refined info available
        #zero in the mask means there is refinement available
        a=b
        left -= chunk
    assert left==0
    return idc,ioctch
    
nchem=8+2
dtyp = na.dtype(">i4,>i8,>i8"+",>%sf4"%(nchem)+ \
                ",>%sf4"%(2)+",>i4")
def _read_art_child(f, level_child_offsets,level,nLevel,field):
    pos=f.tell()
    f.seek(level_child_offsets[level])
    arr = na.fromfile(f, dtype='>f', count=nLevel * 8)
    arr = arr.reshape((nLevel,16), order="F")
    arr = arr[3:-1,:].astype("float64")
    f.seek(pos)
    return arr[field,:]

def _skip_record(f):
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    f.seek(s[0], 1)
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))

def _read_frecord(f,fmt):
    s1 = struct.unpack('>i', f.read(struct.calcsize('>i')))[0]
    count = s1/na.dtype(fmt).itemsize
    ss = na.fromfile(f,fmt,count=count)
    s2 = struct.unpack('>i', f.read(struct.calcsize('>i')))[0]
    assert s1==s2
    return ss


def _read_record(f,fmt=None):
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))[0]
    ss = f.read(s)
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    if fmt is not None:
        return struct.unpack(ss,fmt)
    return ss

def _read_record_size(f):
    pos = f.tell()
    s = struct.unpack('>i', f.read(struct.calcsize('>i')))
    f.seek(pos)
    return s[0]

def _read_struct(f,structure,verbose=False):
    vals = {}
    for format,name in structure:
        size = struct.calcsize(format)
        (val,) = struct.unpack(format,f.read(size))
        vals[name] = val
        if verbose: print "%s:\t%s\t (%d B)" %(name,val,f.tell())
    return vals



#All of these functions are to convert from hydro time var to 
#proper time
sqrt = na.sqrt
sign = na.sign

def find_root(f,a,b,tol=1e-6):
    c = (a+b)/2.0
    last = -na.inf
    assert(sign(f(a)) != sign(f(b)))  
    while na.abs(f(c)-last) > tol:
        last=f(c)
        if sign(last)==sign(f(b)):
            b=c
        else:
            a=c
        c = (a+b)/2.0
    return c

def quad(fintegrand,xmin,xmax,n=1e4):
    spacings = na.logspace(na.log10(xmin),na.log10(xmax),n)
    integrand_arr = fintegrand(spacings)
    val = na.trapz(integrand_arr,dx=na.diff(spacings))
    return val

def a2b(at,Om0=0.27,Oml0=0.73,h=0.700):
    def f_a2b(x):
        val = 0.5*sqrt(Om0) / x**3.0
        val /= sqrt(Om0/x**3.0 +Oml0 +(1.0 - Om0-Oml0)/x**2.0)
        return val
    #val, err = si.quad(f_a2b,1,at)
    val = quad(f_a2b,1,at)
    return val

def b2a(bt,**kwargs):
    #converts code time into expansion factor 
    #if Om0 ==1and OmL == 0 then b2a is (1 / (1-td))**2
    #if bt < -190.0 or bt > -.10:  raise 'bt outside of range'
    f_b2a = lambda at: a2b(at,**kwargs)-bt
    return find_root(f_b2a,1e-4,1.1)
    #return so.brenth(f_b2a,1e-4,1.1)
    #return brent.brent(f_b2a)

def a2t(at,Om0=0.27,Oml0=0.73,h=0.700):
    integrand = lambda x : 1./(x*sqrt(Oml0+Om0*x**-3.0))
    #current_time,err = si.quad(integrand,0.0,at,epsabs=1e-6,epsrel=1e-6)
    current_time = quad(integrand,1e-4,at)
    #spacings = na.logspace(-5,na.log10(at),1e5)
    #integrand_arr = integrand(spacings)
    #current_time = na.trapz(integrand_arr,dx=na.diff(spacings))
    current_time *= 9.779/h
    return current_time

def b2t(tb,n = 1e2,logger=None,**kwargs):
    tb = na.array(tb)
    if type(tb) == type(1.1): 
        return a2t(b2a(tb))
    if tb.shape == (): 
        return a2t(b2a(tb))
    if len(tb) < n: n= len(tb)
    age_min = a2t(b2a(tb.max(),**kwargs),**kwargs)
    age_max = a2t(b2a(tb.min(),**kwargs),**kwargs)
    tbs  = -1.*na.logspace(na.log10(-tb.min()),
                          na.log10(-tb.max()),n)
    ages = []
    for i,tbi in enumerate(tbs):
        ages += a2t(b2a(tbi)),
        if logger: logger(i)
    ages = na.array(ages)
    fb2t = na.interp(tb,tbs,ages)
    #fb2t = interp1d(tbs,ages)
    return fb2t

def spread_ages(ages,logger=None,spread=1.0e7*365*24*3600):
    #stars are formed in lumps; spread out the ages linearly
    da= na.diff(ages)
    assert na.all(da<=0)
    #ages should always be decreasing, and ordered so
    agesd = na.zeros(ages.shape)
    idx, = na.where(da<0)
    idx+=1 #mark the right edges
    #spread this age evenly out to the next age
    lidx=0
    lage=0
    for i in idx:
        n = i-lidx #n stars affected
        rage = ages[i]
        lage = max(rage-spread,0.0)
        agesd[lidx:i]=na.linspace(lage,rage,n)
        lidx=i
        #lage=rage
        if logger: logger(i)
    #we didn't get the last iter
    i=ages.shape[0]-1
    n = i-lidx #n stars affected
    rage = ages[i]
    lage = max(rage-spread,0.0)
    agesd[lidx:i]=na.linspace(lage,rage,n)
    return agesd
