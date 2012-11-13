"""
Definitions for ARTIO files

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

# These functions are ARTIO-specific

def artio_header(hvals):
    header = ( ('ncpu', 1, 'i'),
               ('ndim', 1, 'i'),
               ('nx', 3, 'i'),
               ('nlevelmax', 1, 'i'),
               ('ngridmax', 1, 'i'),
               ('nboundary', 1, 'i'),
               ('ngrid_current', 1, 'i'),
               ('boxlen', 1, 'd'),
               ('nout', 3, 'I')
              )
    yield header
    noutput, iout, ifout = hvals['nout']
    next_set = ( ('tout', noutput, 'd'),
                 ('aout', noutput, 'd'),
                 ('t', 1, 'd'),
                 ('dtold', hvals['nlevelmax'], 'd'),
                 ('dtnew', hvals['nlevelmax'], 'd'),
                 ('nstep',  2, 'i'),
                 ('stat', 3, 'd'),
                 ('cosm', 7, 'd'),
                 ('timing', 5, 'd'),
                 ('mass_sph', 1, 'd') )
    yield next_set
    tree_header = ( ('headl', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                    ('taill', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                    ('numbl', hvals['nlevelmax'] * hvals['ncpu'], 'i'),
                  )
    yield tree_header
