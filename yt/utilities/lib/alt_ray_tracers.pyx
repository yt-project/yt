"""




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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
def cylindrical_ray_trace(np.ndarray[np.float64_t, ndim=1] p1, 
                          np.ndarray[np.float64_t, ndim=1] p2, 
                          np.ndarray[np.float64_t, ndim=2] left_edges, 
                          np.ndarray[np.float64_t, ndim=2] right_edges):
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
    cdef int i, I
    cdef np.float64_t a, b, bsqrd, twoa
    cdef np.ndarray[np.float64_t, ndim=1] p1cart, p2cart, dpcart, t, s, \
                                          rleft, rright, zleft, zright, \
                                          cleft, cright, thetaleft, thetaright, \
                                          tmleft, tpleft, tmright, tpright, tsect
    cdef np.ndarray[np.int64_t, ndim=1, cast=True] inds, tinds, sinds
    cdef np.ndarray[np.float64_t, ndim=2] xyz, rztheta, ptemp, b1, b2, dsect

    # set up  points
    ptemp = np.array([p1, p2])
    ptemp = _cyl2cart(ptemp)
    p1cart = ptemp[0]
    p2cart = ptemp[1]
    dpcart = p2cart - p1cart

    # set up components
    rleft = left_edges[:,0]
    rright = right_edges[:,0]
    zleft = left_edges[:,1]
    zright = right_edges[:,1]

    a = (dpcart[0]**2) + (dpcart[1]**2)
    b = (2*dpcart[0]*p1cart[0]) + (2*dpcart[1]*p1cart[1])
    cleft = ((p1cart[0]**2) + (p1cart[1]**2)) - rleft**2
    cright = ((p1cart[0]**2) + (p1cart[1]**2)) - rright**2
    twoa = 2*a
    bsqrd = b**2

    # Compute positive and negative times and associated masks
    I = np.intp(left_edges.shape[0])
    tmleft = np.empty(I, dtype='float64')
    tpleft = np.empty(I, dtype='float64')
    tmright = np.empty(I, dtype='float64')
    tpright = np.empty(I, dtype='float64')
    for i in range(I):
        tmleft[i] = (-b - math.sqrt(bsqrd - 4*a*cleft[i])) / twoa
        tpleft[i] = (-b + math.sqrt(bsqrd - 4*a*cleft[i])) / twoa  
        tmright[i] = (-b - math.sqrt(bsqrd - 4*a*cright[i])) / twoa
        tpright[i] = (-b + math.sqrt(bsqrd - 4*a*cright[i])) / twoa

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
        inds = np.arange(np.intp(I))
        thetaleft = np.empty(I)
        thetaleft.fill(p1[2])
        thetaright = np.empty(I)
        thetaleft.fill(p2[2])
    else:
        rleft = rleft[inds]
        rright = rright[inds]

        zleft = zleft[inds]
        zright = zright[inds]

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
    b1 = np.concatenate([[rleft,  zright, thetaleft],
                         [rleft,  zleft,  thetaleft],
                         [rleft,  zleft,  thetaleft],
                         [rright, zleft,  thetaright],
                         [rleft,  zleft,  thetaleft],
                         [rright, zright, thetaleft],
                         [rleft,  zright, thetaleft],
                         [rright, zleft,  thetaleft],
                         ], axis=1).T

    b2 = np.concatenate([[rright, zright, thetaright],
                         [rright, zleft,  thetaright], 
                         [rleft,  zright, thetaleft], 
                         [rright, zright, thetaright], 
                         [rleft,  zleft,  thetaright], 
                         [rright, zright, thetaright], 
                         [rleft,  zright, thetaright], 
                         [rright, zleft,  thetaright],
                         ], axis=1).T

    inds = np.concatenate([inds, inds, inds, inds, inds, inds, inds, inds])

    # find intersections and compute return values
    tsect, dsect = _cart_intersect(p1cart, p2cart, _cyl2cart(b1), _cyl2cart(b2))
    tmask = np.logical_and(0.0<=tsect, tsect<=1.0)
    tsect, tinds = np.unique(tsect[tmask], return_index=True)
    inds = inds[tmask][tinds]
    xyz = dsect[tmask][tinds]
    s = np.sqrt(((xyz - p1cart)**2).sum(axis=1))
    s, sinds = np.unique(s, return_index=True)
    inds = inds[sinds]
    xyz = xyz[sinds]
    t = s/np.sqrt((dpcart**2).sum())
    sinds = s.argsort()
    s = s[sinds]
    t = t[sinds]
    inds = inds[sinds]
    xyz = xyz[sinds]
    rztheta = np.concatenate([np.sqrt(xyz[:,0]**2 + xyz[:,1]**2)[:,np.newaxis], 
                              xyz[:,2:3],
                              np.arctan2(xyz[:,1], xyz[:,0])[:,np.newaxis]], axis=1)
    return t, s, rztheta, inds
    #rztheta[:,2] = 0.0 + (rztheta[:,2] - np.pi*3/2)%(2*np.pi)
