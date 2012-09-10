"""
Transformations between coordinate systems

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: JS Oishi <jsoishi@astro.berkeley.edu>
Organization: UC Berkeley
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk, J. S. Oishi.  All Rights Reserved.

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
from yt.funcs import *

from yt.utilities.linear_interpolators import \
    TrilinearFieldInterpolator

def spherical_regrid(pf, nr, ntheta, nphi, rmax, fields,
                     center=None, smoothed=True):
    """
    This function takes a parameter file (*pf*) along with the *nr*, *ntheta*
    and *nphi* points to generate out to *rmax*, and it grids *fields* onto
    those points and returns a dict.  *center* if supplied will be the center,
    otherwise the most dense point will be chosen.  *smoothed* governs whether
    regular covering grids or smoothed covering grids will be used.
    """
    mylog.warning("This code may produce some artifacts of interpolation")
    mylog.warning("See yt/extensions/coordinate_transforms.py for plotting information")
    if center is None: center = pf.h.find_max("Density")[1]
    fields = ensure_list(fields)
    r,theta,phi = np.mgrid[0:rmax:nr*1j,
                           0:np.pi:ntheta*1j,
                           0:2*np.pi:nphi*1j]
    new_grid = dict(r=r, theta=theta, phi=phi)
    new_grid['x'] = r*np.sin(theta)*np.cos(phi) + center[0]
    new_grid['y'] = r*np.sin(theta)*np.sin(phi) + center[1]
    new_grid['z'] = r*np.cos(theta)             + center[2]
    sphere = pf.h.sphere(center, rmax)
    return arbitrary_regrid(new_grid, sphere, fields, smoothed)

def arbitrary_regrid(new_grid, data_source, fields, smoothed=True):
    """
    This function accepts a dict of points 'x', 'y' and 'z' and a data source
    from which to interpolate new points, along with a list of fields it needs
    to regrid onto those xyz points.  It then returns interpolated points.
    This has not been well-tested other than for regular spherical regridding.
    """
    fields = ensure_list(fields)
    new_grid['handled'] = np.zeros(new_grid['x'].shape, dtype='bool')
    for field in fields:
        new_grid[field] = np.zeros(new_grid['x'].shape, dtype='float64')
    grid_order = np.argsort(data_source.gridLevels)
    ng = len(data_source._grids)

    for i,grid in enumerate(data_source._grids[grid_order][::-1]):
        mylog.info("Regridding grid % 4i / % 4i (%s - %s)", i, ng, grid.id, grid.Level)
        cg = grid.retrieve_ghost_zones(1, fields, smoothed=smoothed)

        # makes x0,x1,y0,y1,z0,z1
        bounds = np.concatenate(zip(cg.left_edge, cg.right_edge)) 

        
        # Now we figure out which of our points are inside this grid
        # Note that we're only looking at the grid, not the grid-with-ghost-zones
        point_ind = np.ones(new_grid['handled'].shape, dtype='bool') # everything at first
        for i,ax in enumerate('xyz'): # i = 0,1,2 ; ax = x, y, z
            # &= does a logical_and on the array
            point_ind &= ( ( grid.LeftEdge[i] <= new_grid[ax]      )
                         & ( new_grid[ax]     <= grid.RightEdge[i] ) )
        point_ind &= (new_grid['handled'] == False) # only want unhandled points

        # If we don't have any, we can just leave
        if point_ind.sum() == 0: continue

        # because of the funky way the interpolator takes points, we have to make a
        # new dict of just the points inside this grid
        point_grid = {'x' : new_grid['x'][point_ind],
                      'y' : new_grid['y'][point_ind],
                      'z' : new_grid['z'][point_ind]}

        # Now we know which of the points in new_grid are inside this grid
        for field in fields:
            interpolator = TrilinearFieldInterpolator(
                cg[field],bounds,['x','y','z'])
            new_grid[field][point_ind] = interpolator(point_grid)

        new_grid['handled'][point_ind] = True

    mylog.info("Finished with %s dangling points",
        new_grid['handled'].size - new_grid['handled'].sum())

    return new_grid

"""
# The following will work to plot through different slices:

import pylab
for i in range(n_theta):
    print "Doing % 3i / % 3i" % (i, n_theta)
    pylab.clf()
    ax=pylab.subplot(1,1,1, projection="polar", aspect=1.)
    ax.pcolormesh(phi[:,i,:], r[:,i,:],
                  np.log10(sph_grid[field][:,i,:]))
    pylab.savefig("polar/latitude_%03i.png" % i)

for i in range(n_phi):
    print "Doing % 3i / % 3i" % (i, n_phi)
    pylab.clf()
    ax=pylab.subplot(1,1,1, projection="polar", aspect=1.)
    ax.pcolormesh(theta[:,:,i], r[:,:,i],
                  np.log10(sph_grid[field][:,:,i]))
    pylab.savefig("polar/longitude_%03i.png" % i)
"""
