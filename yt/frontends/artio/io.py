"""
ARTIO-specific IO

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
from collections import defaultdict
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from .definitions import yt_to_art
import yt.utilities.fortran_utils as fpu
import cStringIO

class IOHandlerARTIO(BaseIOHandler):
    _data_style = "artio"

    def _read_fluid_selection(self, chunks, selector, fields ):
    	print 'reading in fluid data'
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        cp = 0
        for onechunk in chunks:
            for artchunk in onechunk.objs :
                rv = artchunk.fill(fields)  
                for f in fields:
                    tr[f].resize(cp+artchunk.data_size)
                    tr[f][cp:cp+artchunk.data_size] = rv.pop(f)
                cp += artchunk.data_size 
        return tr

    def _read_particle_selection(self, chunks, selector, fields):
        # http://yt-project.org/doc/analyzing/particles.html
        print "reading particle data"
        
        # get particle types
        for f in fields :
            if particles_yt_to_art[f[0]] not in self.particle_types :
                assert (particles_yt_to_art.has_key(f[0])) #particle types must exist in ART
                self.particle_types.append(yt_to_art[f[0]])
        accessed_species = self.particle_types #duplicate fields added in caller

    	print 'io.py particle fields ',fields
        sys.exit(1)
        cp = 0
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        for onechunk in chunks:
            rv = onechunk.fill_particles( accessed_species, fields)

            # make sure size of all fieldnames within a fieldtype is the same :
            prev_fieldtype = None
            prev_fieldname_size  = 0
            for fieldtype, fieldname in fields:
                fieldname_size = len(rv[fieldtype,fieldname])
                if prev_fieldtype == field_type and fieldname_size != prev_fieldname_size :
                    print 'size varies between fields! exiting'
                    sys.exit(1)
                prev_fieldtype = fieldtype
                prev_fieldname_size = len(rv[fieldtype,fieldname])
                
            for f in fields:
                cp_newsize = cp[f]+len(rv[f])
                tr[f].resize(cp_newsize) 
                tr[f][cp[f]:cp_newsize] = rv.pop(f)
                cp[f] = cp_newsize
        return tr
