"""
A base-class representing an astrophysical object

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

from .astrophysical_object import \
    AstrophysicalObject, identification_method, correlation_method
    
class ClumpedRegion(AstrophysicalObject):
    _type_name = "clumped_region"
    def __init__(self, data_source):
        AstrophysicalObject.__init__(self, data_source)

@identification_method("clumped_region", "level_set")
def clumps(obj, field, min_val):
    ds = obj.data_source
    mi, ma = ds.quantities["Extrema"](field)[0]
    cls = obj.data_source.extract_connected_sets(field, 1, min_val, ma)
    return [ClumpedRegion(o) for o in cls[1][0]]
