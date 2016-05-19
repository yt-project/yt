"""
Create a Catmull-Rom spline.



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import sys

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
    result = np.zeros(N)
    if not sorted:
        isort = np.argsort(old_x)
        old_x = old_x[isort]
        old_y = old_y[isort]
    # Floor/ceiling of values outside of the original data
    new_x = np.minimum(new_x, old_x[-1])
    new_x = np.maximum(new_x, old_x[0])
    ind = np.searchsorted(old_x, new_x)
    im2 = np.maximum(ind-2, 0)
    im1 = np.maximum(ind-1, 0)
    ip1 = np.minimum(ind+1, ndata-1)
    for i in range(N):
        if ind[i] != im1[i]:
            u = (new_x[i] - old_x[im1[i]]) / (old_x[ind[i]] - old_x[im1[i]])
        elif ind[i] == im1[i]:
            u = 0
        else:
            print ("Bad index during interpolation?")
            sys.exit()
        b0 = -tension * u + 2*tension * u**2 - tension * u**3
        b1 = 1.0 + (tension-3) * u**2 + (2-tension) * u**3
        b2 = tension * u + (3 - 2*tension) * u**2 + (tension-2) * u**3
        b3 = -tension * u**2 + tension * u**3
        result[i] = b0 * old_y[im2[i]] + b1 * old_y[im1[i]] + \
                    b2 * old_y[ind[i]] + b3 * old_y[ip1[i]]
    return result
    
