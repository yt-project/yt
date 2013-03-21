"""
Data structures helper for Enzo

Author: Stephen Skory <s@skory.us>
Affiliation: Univ of Colorado
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2013 Matthew Turk.  All Rights Reserved.

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

from yt.utilities import hdf5_light_reader

def get_field_names_helper(filename, id, results):
    try:
        names= hdf5_light_reader.ReadListOfDatasets(
                    filename, "/Grid%08i" % id)
        results.put((names, "Grid %s has: %s" % (grid.id, names)))
    except:
        results.put((None, "Grid %s is a bit funky?" % id))
