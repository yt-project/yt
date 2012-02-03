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
import pdb

from yt.utilities.io_handler import \
    BaseIOHandler
import numpy as na

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.amr_utils as au

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

    def preload_level(self, level):
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
        arr = arr[3:-1,:].astype("float64")
        self.level_data[level] = arr

    def preload_root_level(self):
        f = open(self.filename, 'rb')
        f.seek(self.level_offsets[0] + 4) # Ditch the header
        ncells = self.level_info[0]
        #pdb.set_trace()
        nhvals = ncells * (self.nhydro_vars) # 0 vars, 0 pads
        hvar = na.fromfile(f, dtype='>f', count=nhvals).astype("float64")
        hvar = hvar.reshape((self.nhydro_vars, ncells), order="F")
        na.fromfile(f,dtype='>i',count=2) #throw away the pads
        nvars = ncells * (2) # 0 vars, 0 pads
        var = na.fromfile(f, dtype='>f', count=nvars).astype("float64")
        var = var.reshape((2, ncells), order="F")
        arr = na.concatenate((hvar,var))
        self.level_data[0] = arr

    def clear_level(self, level):
        self.level_data.pop(level, None)
        
    def _read_data_set(self, grid, field):
        pf = grid.pf
        field_id = grid.pf.h.field_list.index(field)
        if grid.Level == 0: # We only have one root grid
            self.preload_level(0)
            tr = self.level_data[0][field_id,:].reshape(
                    pf.domain_dimensions, order="F").copy()
            return tr.swapaxes(0, 2)
        tr = na.zeros(grid.ActiveDimensions, dtype='float64')
        filled = na.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                self.preload_level(g.Level)
                #print "Filling %s from %s (%s)" % (grid, g, g.Level)
                to_fill -= au.read_art_grid(field_id, 
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, self.level_data[g.Level],
                        g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

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
    Level = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iNOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
    iHOLL = na.zeros(MaxLevelNow+1 - MinLev, dtype='i')
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
    #in the future, break this up into large chunks
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
    le = le/2**(root_level-1-level)-1
    f.seek(pos)
    return le,fl,iocts,nLevel

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
