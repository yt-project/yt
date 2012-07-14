"""
Create a Catmull-Rom spline.

Author: John Wise <jwise@physics.gatech.edu>
Affiliation: Georgia Tech
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 John Wise.  All Rights Reserved.

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

def create_spline(old_x, old_y, new_x, tension=0.5, sorted=False):
    """
Inputs:
  old_x: array of floats
    Original x-data to be fit with a Catmull-Rom spline
  old_y: array of floats
    Original y-data to be fit with a Catmull-Rom spline
  new_x: array of floats
    interpolate to these x-coordinates
  tension: float, optional
    controls the tension at the specified coordinates
  sorted: boolean, optional
    If True, then the old_x and old_y arrays are sorted, and then this routine
    does not try to sort the coordinates
Outputs:
  result: array of floats
    interpolated y-coordinates
    """
    ndata = len(old_x)
    N = len(new_x)
    result = na.zeros(N)
    if not sorted:
        isort = na.argsort(old_x)
        old_x = old_x[isort]
        old_y = old_y[isort]
    # Floor/ceiling of values outside of the original data
    new_x = na.minimum(new_x, old_x[-1])
    new_x = na.maximum(new_x, old_x[0])
    ind = na.searchsorted(old_x, new_x)
    im2 = na.maximum(ind-2, 0)
    im1 = na.maximum(ind-1, 0)
    ip1 = na.minimum(ind+1, ndata-1)
    for i in range(N):
        if ind[i] != im1[i]:
            u = (new_x[i] - old_x[im1[i]]) / (old_x[ind[i]] - old_x[im1[i]])
        elif ind[i] == im1[i]:
            u = 0
        else:
            print "Bad index during interpolation?"
            sys.exit()
        b0 = -tension * u + 2*tension * u**2 - tension * u**3
        b1 = 1.0 + (tension-3) * u**2 + (2-tension) * u**3
        b2 = tension * u + (3 - 2*tension) * u**2 + (tension-2) * u**3
        b3 = -tension * u**2 + tension * u**3
        result[i] = b0 * old_y[im2[i]] + b1 * old_y[im1[i]] + \
                    b2 * old_y[ind[i]] + b3 * old_y[ip1[i]]
    return result
    
