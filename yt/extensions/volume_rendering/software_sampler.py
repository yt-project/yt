"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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
from yt.extensions.volume_rendering import *
from yt.funcs import *

def direct_ray_cast(pf, L, center, W, Nvec, tf, 
                    partitioned_grids = None):
    center = na.array(center, dtype='float64')

    # This just helps us keep track of stuff, and it's cheap
    cp = pf.h.cutting(L, center)
    back_center = center - cp._norm_vec * W
    front_center = center + cp._norm_vec * W
    cylinder = pf.h.disk(back_center, L, na.sqrt(2)*W, 2*W)

    if partitioned_grids == None:
        partitioned_grids = partition_all_grids(cylinder._grids,
                                    eval_func = lambda g: na.any(cylinder._get_point_indices(g)))
    #partitioned_grids = partition_all_grids(pf.h.grids)

    LE = (na.array([grid.LeftEdge for grid in partitioned_grids]) - back_center) * cp._norm_vec
    RE = (na.array([grid.RightEdge for grid in partitioned_grids]) - back_center) * cp._norm_vec
    DL = na.sum(LE, axis=1); del LE
    DR = na.sum(RE, axis=1); del RE
    dist = na.minimum(DL, DR)
    ind = na.argsort(dist)
    
    image = na.zeros((Nvec,Nvec,4), dtype='float64', order='F')
    image[:,:,3] = 1.0

    # Now we need to generate regular x,y,z values in regular space for our vector
    # starting places.
    px, py = na.mgrid[-W:W:Nvec*1j,-W:W:Nvec*1j]
    xv = cp._inv_mat[0,0]*px + cp._inv_mat[0,1]*py + cp.center[0]
    yv = cp._inv_mat[1,0]*px + cp._inv_mat[1,1]*py + cp.center[1]
    zv = cp._inv_mat[2,0]*px + cp._inv_mat[2,1]*py + cp.center[2]
    vectors = na.array([xv, yv, zv], dtype='float64').transpose()
    vectors = vectors.copy('F')
    xp0, xp1 = px.min(), px.max()
    yp0, yp1 = py.min(), py.max()

    ng = partitioned_grids.size
    norm_vec = cp._norm_vec
    norm_vec = cp._norm_vec * (2.0*W)
    hit = 0
    tnow = time.time()
    every = na.ceil(len(partitioned_grids) / 100.0)

    vp = VectorPlane(vectors, norm_vec, back_center,
                     (xp0, xp1, yp0, yp1), image, cp._x_vec, cp._y_vec)

    tfp = TransferFunctionProxy(tf)

    pbar = get_pbar("Ray casting ", len(partitioned_grids))
    for i,g in enumerate(partitioned_grids[ind]):
        if (i % every) == 0: 
            pbar.update(i)
        pos = g.cast_plane(tfp, vp)
    pbar.finish()

    return partitioned_grids, image, vectors, norm_vec, pos
