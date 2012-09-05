"""Ray tracing algorithms for rays in non-cartesian coordinate systems.

Author: Anthony Scopatz <scopatz@gmail.com>
Affiliation: University of Chicago
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Anthony Scopatz.  All Rights Reserved.

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
cimport numpy as np
cimport cython
cimport libc.math as math

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _cyl2cart(np.ndarray[np.float64_t, ndim=2] x):
    """Converts points in cylindrical coordinates to cartesian points."""
    # NOTE this should be removed once the coord interface comes online
    cdef int i, I
    cdef np.ndarray[np.float64_t, ndim=2] xcart 
    I = x.shape[0]
    xcart = np.empty((I, x.shape[1]), dtype='float64')
    for i in range(I):
        xcart[i,0] = x[i,0] * math.cos(x[i,2])
        xcart[i,1] = x[i,0] * math.sin(x[i,2])
        xcart[i,2] = x[i,1]
    return xcart

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _cart_intersect(np.ndarray[np.float64_t, ndim=1] a,
                    np.ndarray[np.float64_t, ndim=1] b,
                    np.ndarray[np.float64_t, ndim=2] c,
                    np.ndarray[np.float64_t, ndim=2] d):
    """Finds the times and locations of the lines defined by a->b and 
    c->d in cartesian space.  a and b must be 1d points.  c and d must 
    be 2d arrays of equal shape whose second dim is the same length as
    a and b.
    """
    cdef int i, I
    cdef np.ndarray[np.float64_t, ndim=1] r, s, t
    cdef np.ndarray[np.float64_t, ndim=2] loc
    I = c.shape[0]
    shape = (I, c.shape[1])
    t = np.empty(I, dtype='float64')
    loc = np.empty(shape, dtype='float64')

    r = b - a
    for i in range(I):
        s = d[i] - c[i]
        t[i] = (np.cross(c[i] - a, s).sum()) / (np.cross(r, s).sum())
        loc[i] = a + t[i]*r
    return t, loc


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def clyindrical_ray_trace(np.ndarray[np.float64_t, ndim=1] p1, 
                          np.ndarray[np.float64_t, ndim=1] p2, 
                          np.ndarray[np.float64_t, ndim=2] left_edges, 
                          np.ndarray[np.float64_t, ndim=2] right_edges)

    """Computes straight (cartesian) rays being traced through a 
    cylindrical geometry.

    Parameters
    ----------
    p1 : length 3 float ndarray
        start point for ray
    p2 : length 3 float ndarray
        stop point for ray
    left_edges : 2d ndarray
        left edges of grid cells
    right_edges : 2d ndarray
        right edges of grid cells

    Returns
    -------
    t : 1d float ndarray
        ray parametric time on range [0,1]
    s : 1d float ndarray
        ray parametric distance on range [0,len(ray)]
    rztheta : 2d float ndarray
        ray grid cell intersections in cylidrical coordinates
    inds : 1d int ndarray
        indexes into the grid cells which the ray crosses in order.

    """
    int i, I
    cdef np.float64_t a, b,b2, twoa
    cdef np.ndarray[np.float64_t, ndim=1] dp, p1cart, p2cart, dpcart, t, s, 
                                          rleft, rright, zleft, zright, 
                                          cleft, cright, thetaleft, thetaright,
                                          tmleft, tpleft, tmright, tpright
    cdef np.ndarray[np.int32_t, ndim=1] inds
    cdef np.ndarray[np.float64_t, ndim=2] rztheta, ptemp, b1, b2

    # set up  points
    dp = p2 - p1
    ptemp = np.array([p1, p2])
    ptemp = _cyl2cart(ptemp)
    p1cart = ptemp[0]
    p2cart = ptemp[1]
    dpcart = p2cart - p1cart

    # set up components
    rleft = left_edge[:,0]
    rright = right_edge[:,0]
    zleft = left_edge[:,1]
    zright = right_edge[:,1]

    a = (dpcart[0]**2) + (dpcart[1]**2)
    b = (2*dpcart[0]*p1cart[0]) + (2*dpcart[1]*p1cart[1])
    cleft = ((p1cart[0]**2) + (p1cart[1]**2)) - rleft**2
    cright = ((p1cart[0]**2) + (p1cart[1]**2)) - rright**2
    twoa = 2*a
    b2 = b**2

    # Compute positive and negative times and associated masks
    I = left_edge.shape[0]
    tmleft = np.empty(I, dtype='float64')
    tpleft = np.empty(I, dtype='float64')
    tmright = np.empty(I, dtype='float64')
    tpright = np.empty(I, dtype='float64')
    for i in range(I):
        tmleft[i] = (-b - math.sqrt(b2 - 4*a*cleft[i])) / twoa
        tpleft[i] = (-b + math.sqrt(b2 - 4*a*cleft[i])) / twoa  
        tmright[i] = (-b - math.sqrt(b2 - 4*a*cright[i])) / twoa
        tpright[i] = (-b + math.sqrt(b2 - 4*a*cright[i])) / twoa

    tmmright = np.logical_and(~np.isnan(tmright), rright <= p1[0])
    tpmright = np.logical_and(~np.isnan(tpright), rright <= p2[0])

    tmmleft = np.logical_and(~np.isnan(tmleft), rleft <= p1[0])
    tpmleft = np.logical_and(~np.isnan(tpleft), rleft <= p2[0])

    # compute first cut of indexes and thetas, which 
    # have been filtered by those values for which intersection
    # times are impossible (see above masks). Note that this is
    # still independnent of z.
    inds = np.unique(np.concatenate([np.argwhere(tmmleft).flat, 
                                     np.argwhere(tpmleft).flat, 
                                     np.argwhere(tmmright).flat, 
                                     np.argwhere(tpmright).flat,]))
    if 0 == inds.shape[0]:
        inds = np.arange(I)
        thetaleft = np.empty(I)
        thetaleft.fill(p1[2])
        thetaright = np.empty(I)
        thetaleft.fill(p2[2])
    else:
        thetaleft = np.arctan2((p1cart[1] + tmleft[inds]*dpcart[1]), 
                               (p1cart[0] + tmleft[inds]*dpcart[0]))
        nans = np.isnan(thetaleft)
        thetaleft[nans] = np.arctan2((p1cart[1] + tpleft[inds[nans]]*dpcart[1]), 
                                     (p1cart[0] + tpleft[inds[nans]]*dpcart[0]))

        thetaright = np.arctan2((p1cart[1] + tmright[inds]*dpcart[1]), 
                                (p1cart[0] + tmright[inds]*dpcart[0]))
        nans = np.isnan(thetaright)
        thetaright[nans] = np.arctan2((p1cart[1] + tpright[inds[nans]]*dpcart[1]), 
                                      (p1cart[0] + tpright[inds[nans]]*dpcart[0]))

    # Set up the cell boundary arrays    
    b1 = np.concatenate()              np.array([rleft[ind], zright[ind], thetaleft]).T
c = np.append(c, np.array([rleft[ind], zleft[ind], thetaleft]).T, axis=0)
c = np.append(c, np.array([rleft[ind], zleft[ind], thetaleft]).T, axis=0)
c = np.append(c, np.array([rright[ind], zleft[ind], thetaright]).T, axis=0)
c = np.append(c, np.array([rleft[ind], zleft[ind], thetaleft]).T, axis=0)
c = np.append(c, np.array([rright[ind], zright[ind], thetaleft]).T, axis=0)
c = np.append(c, np.array([rleft[ind], zright[ind], thetaleft]).T, axis=0)
c = np.append(c, np.array([rright[ind], zleft[ind], thetaleft]).T, axis=0)

d =              np.array([rright[ind], zright[ind], thetaright]).T
d = np.append(d, np.array([rright[ind], zleft[ind], thetaright]).T, axis=0)
d = np.append(d, np.array([rleft[ind], zright[ind], thetaleft]).T, axis=0)
d = np.append(d, np.array([rright[ind], zright[ind], thetaright]).T, axis=0)
d = np.append(d, np.array([rleft[ind], zleft[ind], thetaright]).T, axis=0)
d = np.append(d, np.array([rright[ind], zright[ind], thetaright]).T, axis=0)
d = np.append(d, np.array([rleft[ind], zright[ind], thetaright]).T, axis=0)
d = np.append(d, np.array([rright[ind], zleft[ind], thetaright]).T, axis=0)

origind = ind
ind = np.append(ind, origind)
ind = np.append(ind, origind)
ind = np.append(ind, origind)
ind = np.append(ind, origind)
ind = np.append(ind, origind)
ind = np.append(ind, origind)
ind = np.append(ind, origind)

"""\ 
tsec, intsec = intersect(a, b, c, d)
print tsec
print np.isnan(tsec).all(), len(tsec)
tmask = np.logical_and(0.0<=tsec, tsec<=1.0)
tsec, utmask = np.unique(tsec[tmask], return_index=True)
print tsec.max(), tsec.ptp(), len(utmask)
ind = ind[tmask][utmask]
xyz = intsec[tmask][utmask]
s = np.sqrt(((xyz - Ecart)**2).sum(axis=1))
s, smask = np.unique(s, return_index=True)
ind = ind[smask]
xyz = xyz[smask]
t = s/np.sqrt((Dcart**2).sum())
si = s.argsort()
s = s[si]
t = t[si]
ind = ind[si]
xyz = xyz[si]
rztheta = np.array([np.sqrt(xyz[:,0]**2 + xyz[:,1]**2), xyz[:,2], np.arctan2(xyz[:,1], xyz[:,0])]).T
rztheta[:,2] = 0.0 + (rztheta[:,2] - np.pi*3/2)%(2*np.pi)




"""
