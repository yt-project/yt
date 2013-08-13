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


class IOHandlerARTIO(BaseIOHandler):
    _data_style = "artio"

    def _read_fluid_selection(self, chunks, selector, fields):
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        cp = 0
        for onechunk in chunks:
            for artchunk in onechunk.objs:
                rv = artchunk.fill(fields, selector)
                for f in fields:
                    tr[f].resize(cp+artchunk.data_size)
                    tr[f][cp:cp+artchunk.data_size] = rv.pop(f)
                cp += artchunk.data_size
        return tr

    def _read_particle_selection(self, chunks, selector, fields):
        fd = dict((ftuple, []) for ftuple in fields)
        ftypes = set(ftype for (ftype, fname) in fields)
        for ftype in ftypes:
            for ax in 'xyz':
                ftuple = (ftype, "particle_position_%s" % ax)
                if ftuple not in fd:
                    fd[ftuple] = []
        for onechunk in chunks:
            # We're going to read here
            for artchunk in onechunk.objs:
                rv = artchunk.fill_particles(fd.keys())
                # Now we count and also cut
                for ftype in rv:
                    mask = selector.select_points(
                        rv[ftype]["particle_position_x"],
                        rv[ftype]["particle_position_y"],
                        rv[ftype]["particle_position_z"])
                    if mask is None: continue
                    for fname in rv[ftype]:
                        if (ftype, fname) not in fields: continue
                        fd[ftype, fname].append(rv[ftype][fname][mask])
        # This needs to be pre-initialized in case we are later concatenated.
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        for f in fd.keys():
            v = fd.pop(f)
            if f not in fields: continue
            if len(v) == 0: continue
            tr[f] = np.concatenate(v)
        return tr
