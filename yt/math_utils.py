"""
Commonly used mathematical functions.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD Physics/CASS
Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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
import math

def periodic_dist(a, b, period):
    """
    Find the Euclidian periodic distance between two points.
    *a* and *b* are array-like vectors, and *period* is a float or
    array-like value for the periodic size of the volume.
    """
    a = na.array(a)
    b = na.array(b)
    if a.size != b.size: RunTimeError("Arrays must be the same shape.")
    c = na.empty((2, a.size), dtype="float64")
    c[0,:] = abs(a - b)
    c[1,:] = period - abs(a - b)
    d = na.amin(c, axis=0)**2
    return math.sqrt(d.sum())

    