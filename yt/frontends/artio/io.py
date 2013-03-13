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
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from .definitions import yt_to_art

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
        print 'reading in particle data'
        # TODO: determine proper datatype for fields
        tr = dict((ftuple, np.empty(0, dtype='float32')) for ftuple in fields)
        for onechunk in chunks:
            for artchunk in onechunk.objs :
                artchunk.fill_particles(tr, fields)
        return tr
