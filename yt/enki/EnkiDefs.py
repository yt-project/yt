"""
Definitions for Enki

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
# We're going to be setting default vals here

MAX_SPHERE = 10             # Taken from Grid.h
MAX_DEPTH_OF_HIERARCHY = 50 # Taken from macros_and_parameters.h
DEFAULT_GHOST_ZONES    = 3  # Taken from macros_and_parameters.h

TRUE = 1
FALSE = 0

ENZO_SUCCESS = 1
FAIL = 0

FieldTypes = {}

from yt.enki import *

if has_SWIG:
    FieldTypes["Density"] = EnzoInterface.Density
    FieldTypes["TotalEnergy"] = EnzoInterface.TotalEnergy
    FieldTypes["GasEnergy"] = EnzoInterface.InternalEnergy
    FieldTypes["x-velocity"] = EnzoInterface.Velocity1
    FieldTypes["y-velocity"] = EnzoInterface.Velocity2
    FieldTypes["z-velocity"] = EnzoInterface.Velocity3
    FieldTypes["colour"] = EnzoInterface.Metallicity
    FieldTypes["Electron_Density"] = EnzoInterface.ElectronDensity
    FieldTypes["HI_Density"] = EnzoInterface.HIDensity
    FieldTypes["HII_Density"] = EnzoInterface.HIIDensity
    FieldTypes["HeI_Density"] = EnzoInterface.HeIDensity
    FieldTypes["HeII_Density"] = EnzoInterface.HeIIDensity
    FieldTypes["HeIII_Density"] = EnzoInterface.HeIIIDensity
    FieldTypes["HM_Density"] = EnzoInterface.HMDensity
    FieldTypes["H2I_Density"] = EnzoInterface.H2IDensity
    FieldTypes["H2II_Density"] = EnzoInterface.H2IIDensity
    FieldTypes["DI_Density"] = EnzoInterface.DIDensity
    FieldTypes["DII_Density"] = EnzoInterface.DIIDensity
    FieldTypes["HDI_Density"] = EnzoInterface.HDIDensity
    FieldTypes["Metal_Density"] = EnzoInterface.Metallicity

#def LevelHierarchyGen(LevelHierarchy):
#    lh = LevelHierarchy
#    while lh:
#        yield lh
#        lh = lh.NextGridThisLevel

class LevelHierarchyEntryIter:
    # This class is exclusively an iterator!
    def __init__(self,lh):
        self.lh = lh
    def next(self):
        if self.lh == None:
            raise StopIteration
        (oldlh, self.lh) = (self.lh, self.lh.NextGridThisLevel)
        return oldlh
    def __iter__(self):
        return self

class LevelArrayWrapper:
    # this is perhaps not necessary
    def __init__(self):
        self.la = EnzoInterface.new_LevelHierarchyEntryArray(MAX_DEPTH_OF_HIERARCHY)
        for i in range(MAX_DEPTH_OF_HIERARCHY):
            EnzoInterface.LevelHierarchyEntryArray_setitem(self.la, i, None)
    def __getitem__(self, item):
        return LevelHierarchyEntryIter(self.GetLevel(item))
    def __setitem__(self, item, val):
        EnzoInterface.LevelHierarchyEntryArray_setitem(self.la, item, val)
    def GetLevel(self, level):
        return EnzoInterface.LevelHierarchyEntryArray_getitem(self.la, level)
