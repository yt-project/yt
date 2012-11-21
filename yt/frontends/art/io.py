"""
ART-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Chris Moody <cemoody@ucsc.edu>
Affiliation: UCSC
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

import numpy as np
import struct

import os
import os.path

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.lib as au

from yt.frontends.art.definitions import art_particle_field_names

def _count_art_octs(f, offset, 
                   MinLev, MaxLevelNow):
    level_oct_offsets= [0,]
    level_child_offsets= [0,]
    f.seek(offset)
    nchild,ntot=8,0
    Level = np.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iNOLL = np.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iHOLL = np.zeros(MaxLevelNow+1 - MinLev, dtype='i')
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

def _read_art_level_info(f, level_oct_offsets,level,root_level=15):
    pos = f.tell()
    f.seek(level_oct_offsets[level])
    #Get the info for this level, skip the rest
    junk, nLevel, iOct = struct.unpack(
       '>iii', _read_record(f))
    
    #fortran indices start at 1
    
    #Skip all the oct hierarchy data
    le     = np.zeros((nLevel,3),dtype='int64')
    fl     = np.ones((nLevel,6),dtype='int64')
    iocts  = np.zeros(nLevel+1,dtype='int64')
    idxa,idxb = 0,0
    chunk = long(1e6) #this is ~111MB for 15 dimensional 64 bit arrays
    left = nLevel
    while left > 0 :
        this_chunk = min(chunk,left)
        idxb=idxa+this_chunk
        data = np.fromfile(f,dtype='>i',count=this_chunk*15)
        data=data.reshape(this_chunk,15)
        left-=this_chunk
        le[idxa:idxb,:] = data[:,1:4]
        fl[idxa:idxb,1] = np.arange(idxa,idxb)
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

    assert np.unique(iocts).shape[0] == nLevel
    
    #ioct tries to access arrays much larger than le & fl
    #just make sure they appear in the right order, skipping
    #the empty space in between
    idx = np.argsort(iocts)
    
    #now rearrange le & fl in order of the ioct
    le = le[idx]
    fl = fl[idx]

    #left edges are expressed as if they were on 
    #level 15, so no matter what level max(le)=2**15 
    #correct to the yt convention
    #le = le/2**(root_level-1-level)-1

    #try without the -1
    le = le/2**(root_level-2-level)-1

    #now read the hvars and vars arrays
    #we are looking for iOctCh
    #we record if iOctCh is >0, in which it is subdivided
    iOctCh  = np.zeros((nLevel+1,8),dtype='bool')
    
    
    
    f.seek(pos)
    return le,fl,nLevel


def read_particles(file,nstars,Nrow):
    words = 6 # words (reals) per particle: x,y,z,vx,vy,vz
    real_size = 4 # for file_particle_data; not always true?
    np = nstars # number of particles including stars, should come from lspecies[-1]
    np_per_page = Nrow**2 # defined in ART a_setup.h
    num_pages = os.path.getsize(file)/(real_size*words*np_per_page)

    f = np.fromfile(file, dtype='>f4').astype('float32') # direct access
    pages = np.vsplit(np.reshape(f, (num_pages, words, np_per_page)), num_pages)
    data = np.squeeze(np.dstack(pages)).T # x,y,z,vx,vy,vz
    return data[:,0:3],data[:,3:]

def read_stars(file,nstars,Nrow):
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
    return nstars, mass, imass, tbirth, metallicity1, metallicity2

def _read_child_mask_level(f, level_child_offsets,level,nLevel,nhydro_vars):
    f.seek(level_child_offsets[level])
    nvals = nLevel * (nhydro_vars + 6) # 2 vars, 2 pads
    ioctch = np.zeros(nLevel,dtype='uint8')
    idc = np.zeros(nLevel,dtype='int32')
    
    chunk = long(1e6)
    left = nLevel
    width = nhydro_vars+6
    a,b=0,0
    while left > 0:
        chunk = min(chunk,left)
        b += chunk
        arr = np.fromfile(f, dtype='>i', count=chunk*width)
        arr = arr.reshape((width, chunk), order="F")
        assert np.all(arr[0,:]==arr[-1,:]) #pads must be equal
        idc[a:b]    = arr[1,:]-1 #fix fortran indexing
        ioctch[a:b] = arr[2,:]==0 #if it is above zero, then refined info available
        #zero in the mask means there is refinement available
        a=b
        left -= chunk
    assert left==0
    return idc,ioctch
    
nchem=8+2
dtyp = np.dtype(">i4,>i8,>i8"+",>%sf4"%(nchem)+ \
                ",>%sf4"%(2)+",>i4")
def _read_art_child(f, level_child_offsets,level,nLevel,field):
    pos=f.tell()
    f.seek(level_child_offsets[level])
    arr = np.fromfile(f, dtype='>f', count=nLevel * 8)
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
    count = s1/np.dtype(fmt).itemsize
    ss = np.fromfile(f,fmt,count=count)
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
sqrt = np.sqrt
sign = np.sign

def find_root(f,a,b,tol=1e-6):
    c = (a+b)/2.0
    last = -np.inf
    assert(sign(f(a)) != sign(f(b)))  
    while np.abs(f(c)-last) > tol:
        last=f(c)
        if sign(last)==sign(f(b)):
            b=c
        else:
            a=c
        c = (a+b)/2.0
    return c

def quad(fintegrand,xmin,xmax,n=1e4):
    spacings = np.logspace(np.log10(xmin),np.log10(xmax),n)
    integrand_arr = fintegrand(spacings)
    val = np.trapz(integrand_arr,dx=np.diff(spacings))
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
    #spacings = np.logspace(-5,np.log10(at),1e5)
    #integrand_arr = integrand(spacings)
    #current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
    current_time *= 9.779/h
    return current_time

def b2t(tb,n = 1e2,logger=None,**kwargs):
    tb = np.array(tb)
    if type(tb) == type(1.1): 
        return a2t(b2a(tb))
    if tb.shape == (): 
        return a2t(b2a(tb))
    if len(tb) < n: n= len(tb)
    age_min = a2t(b2a(tb.max(),**kwargs),**kwargs)
    age_max = a2t(b2a(tb.min(),**kwargs),**kwargs)
    tbs  = -1.*np.logspace(np.log10(-tb.min()),
                          np.log10(-tb.max()),n)
    ages = []
    for i,tbi in enumerate(tbs):
        ages += a2t(b2a(tbi)),
        if logger: logger(i)
    ages = np.array(ages)
    fb2t = np.interp(tb,tbs,ages)
    #fb2t = interp1d(tbs,ages)
    return fb2t

