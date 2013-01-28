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
import yt.utilities.fortran_utils as fpu
import cStringIO

class IOHandlerARTIO(BaseIOHandler):
    _data_style = "artio"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        tr = dict((ftuple, np.empty(size, dtype='float64')) for ftuple in fields)
        cp = 0
        for onechunk in chunks:
            for subset in onechunk.objs:
                print 'reading from ',fields, subset.domain.grid_fn
                rv = subset.fill(fields) 
                for fieldtype, fieldname in fields:
                    mylog.debug("Filling %s with %s (%0.3e %0.3e) (%s:%s)",
                        fieldname, subset.masked_cell_count, rv[fieldname].min(), 
                        rv[fieldname].max(), cp, cp+subset.masked_cell_count)
                    tr[(fieldtype, fieldname)][cp:cp+subset.masked_cell_count] = rv.pop(fieldname)
                cp += subset.masked_cell_count
        return tr

    def _read_particle_selection(self, chunks, selector, fields):
        # First pass to generate mask

        # FIX need an input for particle type (in fields?)
        # http://yt-project.org/doc/analyzing/particles.html
        # ->creation_time >0 used to indicate star particles
#        accessed_species = ['N-BODY']#,'STAR']
#
#        totsize = 0
#        onesize = 0
#        sizes = {}
#        masks = {}
#        (fieldnames, fieldtypes) = fields
#        print 'fieldnames in io.py', fieldnames
#        print 'quitting from io.py '
#        sys.exit(1)
#        for onechunk in chunks:
#            for subset in onechunk.objs:
#                print 'getting mask from ', subset.domain.part_fn
#                # list of all x positions all y positions and all z positions
#                selection = subset.get_particle_pos(accessed_species, fieldnames) 
#                mask = selector.select_points(selection['x'],
#                            selection['y'], selection['z'])
#                if mask is None: continue
#                onesize = mask.sum()
#                totsize += onesize
#                sizes[id(subset)] = onesize
#                masks[id(subset)] = mask
#
#        # Second pass fills particles where masked
#        tr = dict((f, np.empty(size, dtype="float64")) for f in fields)
#        cp = 0
#        for onechunk in chunks:
#            for subset in onechunk.objs:
#                print 'reading values from', subset.domain.part_fn
#                mask = masks.pop(id(subset), None)
#                if mask is None: continue
#                rv = subset.fill_particles(fields, accessed_species, mask, sizes[id(subset)],mask) 
#                for fieldtype, fieldname in fields:
#                    tr[(fieldtype,fieldname)][cp:cp+sizes[id(subset)]] = rv.pop(fieldname)
#                cp += sizes[id(subset)]

        raise NotImplementedError 
        return tr

