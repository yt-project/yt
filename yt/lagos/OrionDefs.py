"""
Various definitions for various other modules and routines

@author: U{JS Oishi<http://www.jsoishi.org/>}
@organization: U{UC Berkeley<http://www.astro.berkeley.edu/>}
@contact: U{jsoishi@astro.berkeley.edu<mailto:jsoishi@astro.berkeley.edu>}

@todo: Move into yt.Defs, along with enki.EnkiDefs
@license:
  Copyright (C) 2008 J.S. Oishi.  All Rights Reserved.

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

### this assumes EnzoDefs.py has *already* been imported

# converts the Orion inputs file name to the Enzo/yt name expected
# throughout the code. key is Orion name, value is Enzo/yt equivalent
orion2enzoDict = {"amr.n_cell": "TopGridRank",
                  }

orion2ytFieldsDict = {"x-velocity": "xvel",
                      "y-velocity": "yvel",
                      "z-velocity": "zvel",
                      "Density": "density",
                      "Total_Energy": "eden",
                     }
