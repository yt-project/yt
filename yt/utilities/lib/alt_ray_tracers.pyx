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


"""\ 
D = F - E

Ecart = na.array((E[0]*na.cos(E[2]), E[0]*na.sin(E[2]), E[1]))
Fcart = na.array((F[0]*na.cos(F[2]), F[0]*na.sin(F[2]), F[1]))
Dcart = Fcart - Ecart

# <codecell>

rleft = pf.h.grid_left_edge[:,0]
rright = pf.h.grid_right_edge[:,0]
zleft = pf.h.grid_left_edge[:,1]
zright = pf.h.grid_right_edge[:,1]

a = (Dcart**2)[:2].sum()
b = (2*Dcart*Ecart)[:2].sum()
cleft = (Ecart**2)[:2].sum() - rleft**2
cright = (Ecart**2)[:2].sum() - rright**2

tmleft = (-b - na.sqrt(b**2 - 4*a*cleft)) / (2*a)
tpleft = (-b + na.sqrt(b**2 - 4*a*cleft)) / (2*a)  
tmright = (-b - na.sqrt(b**2 - 4*a*cright)) / (2*a)
tpright = (-b + na.sqrt(b**2 - 4*a*cright)) / (2*a)  

tmmright = np.logical_and(~np.isnan(tmright), rright <= E[0])
tpmright = np.logical_and(~np.isnan(tpright), rright <= F[0])

tmmleft = np.logical_and(~np.isnan(tmleft), rleft <= E[0])
tpmleft = np.logical_and(~np.isnan(tpleft), rleft <= F[0])

# <codecell>

ind = np.unique(np.concatenate([np.argwhere(tmmleft).flat, np.argwhere(tpmleft).flat, np.argwhere(tmmright).flat, np.argwhere(tpmright).flat,]))

thetaleft = np.arctan2((Ecart[1] + tmleft[ind]*Dcart[1]), (Ecart[0] + tmleft[ind]*Dcart[0]))
nans = np.isnan(thetaleft)
thetaleft[nans] = np.arctan2((Ecart[1] + tpleft[ind[nans]]*Dcart[1]), (Ecart[0] + tpleft[ind[nans]]*Dcart[0]))

thetaright = np.arctan2((Ecart[1] + tmright[ind]*Dcart[1]), (Ecart[0] + tmright[ind]*Dcart[0]))
nans = np.isnan(thetaleft)
thetaright[nans] = np.arctan2((Ecart[1] + tpright[ind[nans]]*Dcart[1]), (Ecart[0] + tpright[ind[nans]]*Dcart[0]))

#thetaleft += np.pi*3/2
#thetaright += np.pi*3/2

if 0 == len(ind):
    print "Ind len zero"
    I = len(rleft)
    ind = np.arange(I)
    thetaleft = np.empty(I)
    thetaleft.fill(E[2])
    thetaright = np.empty(I)
    thetaleft.fill(F[2])

# <codecell>

a = E
b = F

#c = np.array([rleft[ind], zleft[ind], thetaleft]).T
#d = np.array([rleft[ind], zright[ind], thetaleft]).T

#c = np.array([rright[ind], zleft[ind], thetaright]).T
#d = np.array([rright[ind], zright[ind], thetaright]).T

c =              np.array([rleft[ind], zright[ind], thetaleft]).T
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
