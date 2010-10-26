"""
ART-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerART(BaseIOHandler):
    _data_style = "art"
    #at which position in the child record does the field occur
    field_dict = {'Density':0} 
    
    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)


    def _read_data_set(self,grid,field):
        fn = grid.pf.storage_filename
        # because of reuse promote this handler to pf class?
        min_level = grid.pf.min_level
        max_level = grid.pf.max_level
        nhydro_vars = grid.pf.nhydro_vars 
        grid_level  = grid.Level 
        #grid_idc = (grid.id-1)*8+grid.pf.ncell+1 (verbatim analysis_ART)
        grid_idc = (grid.id)*8+grid.pf.ncell #this gives only the first cell of an oct
        header_offset = grid.pf.child_grid_offset

        record_size = 2+1+nhydro_vars+3 #two pads, idc, hydro vars, 3 vars
        #go past the previous levels
        # skip the header
        offset  = (header_offset 
            #per level:
            # 4 bytes per integer, float or pad
            # first section has three integers + 2 pads
            + 4*5*(grid_level)
            # second section has 13 floats +2 pads per oct 
            # (level_info) is the number of octs per level 
            + 4*15*sum(level_info[:grid_level])
            # after the oct section is the child section.
            # there are 2 pads, 1 integer child ID (idc)
            # then #nhydro_vars of floats + 3 vars 
            #there 8 times as many children as octs
            + 4*(2+1+nhydro_vars+3)*8*sum(level_info[:grid_level-1]))
        pdb.set_trace()
        fh = open(fn,'rb')
        fh.seek(offset)
        dtype='>i,'+'>i,'+'>%df,'%(nhydro_vars+3)+'>i'
        first_idc = na.fromfile(fh,dtype=dtype,count=1)[0][1]
        #this is the first idc of our section. we can guess
        # how many records away our idc is since idc increments 
        # by 1 for almost every record. because idc can increase 
        # by more than 1 sometimes, we will always overestimate
        # the number of records we should skip ahead. so work
        # backwards from our guess
        seek_guess = (grid_idc - first_idc)*record_size*4
        fh.seek(offset+seek_guess)
        while True: 
            record = na.fromfile(fh,dtype=dtype,count=1)[0]
            # make sure that the pad bytes line up
            # otherwise we may have a phase offset or have stepped out  
            # our section
            assert record[0]==record[-1]
            # we better have overestimated the idc
            assert record[1]>= grid_idc 
            if record[1] == grid_idc:
                fh.close()
                return record[2][field_dict[field]] 
            else:
                # in the next iteration we'll read the previous record
                # so rewind the last read, then rewind one further record
                fh.seek(fh.tell()-2*record_size) 
        fh.close()
    
    
