"""
Quantities that can be derived from Enzo data that may also required additional
arguments.  (Standard arguments -- such as the center of a distribution of
points -- are excluded here, and left to the EnzoDerivedFields.)

Is this needlessly complicated?

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *

quantityInfo = {}

def addQuantity(func):
    n = func.func_name
    quantityInfo[n] = tuple([a() for a in func()])

def Phi():
    def GetLog():
        return False
    def GetUnits():
        return "radians"
    def GetFunc():
        def getPhi(data, vector, r_vec = None):
            # Get the declination here
            if r_vec == None:
                r_vec = (data.coords - na.reshape(data.center,(3,1)))
            vec = vector.reshape((3,1))
            nv = vec[0,0] * vec[0,0] + \
                 vec[1,0] * vec[1,0] + \
                 vec[2,0] * vec[2,0]
            nv = na.sqrt(nv)
            nr = r_vec[0,:] * r_vec[0,:] + \
                 r_vec[1,:] * r_vec[1,:] + \
                 r_vec[2,:] * r_vec[2,:]
            nr = na.sqrt(nr)
            dp = r_vec[0,:] * vec[0,0] + \
                 r_vec[1,:] * vec[1,0] + \
                 r_vec[2,:] * vec[2,0]
            phi = na.arccos(dp / (nr * nv))
            return phi
        return getPhi
    def GetHelp():
        return "This returns an array of angles-of-declination."
    return GetLog, GetUnits, GetFunc, GetHelp
addQuantity(Phi)

def AngularMomentumVector():
    def GetLog():
        return False
    def GetUnits():
        return ""
    def GetFunc():
        def getVector(data, weight="CellMassCode"):
            # We want the mass-weighted angular momentum vector for the entire
            # system
            # So we first want to get the angular momentum vectors for every
            # cell:  L = mass * r x v
            #          = cellmass * angular_velocity
            # Then we mass-weight a sum
            #        L' = sum(mass * L) / sum(mass)
            v = na.sum( \
                data[weight]* data["CellMassCode"] * data["AngularVelocity"], axis=1) / \
                na.sum(data[weight])
            vec = v/na.sqrt(na.sum(v*v))
            return v, vec
        return getVector
    def GetHelp():
        return "Returns a single mass-weighted angular momentum vector"
    return GetLog, GetUnits, GetFunc, GetHelp
addQuantity(AngularMomentumVector)

def LDeclination()
    def GetLog()
        return False
    def GetUnits():
        return ""
    def GetFunc():
        def getDeclination(data, weight="CellMassCode"):
