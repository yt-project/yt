"""
Various definitions for various other modules and routines

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/

@todo: Move into yt.Defs, along with enki.EnkiDefs
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

# The number of levels we expect to have at most
MAXLEVEL=48

# Number of levels to toss back to check for maximum
NUMTOCHECK=2

axis_labels = [('y','z'),('x','z'),('x','y')]
axis_names = {0: 'x', 1: 'y', 2: 'z', 4:''}
inv_axis_names = {'x':0,'y':1,'z':2}

vm_axis_names = {0:'x', 1:'y', 2:'z', 3:'dx', 4:'dy'}

# The appropriate axes for which way we are slicing
x_dict = [1,0,0]
y_dict = [2,2,1]

mh = 1.67e-24
mu = 1.22

# All the parameters we read from the parameter file, along with how to convert
# them from a string
parameterDict = {"CosmologyCurrentRedshift": float,
                 "CosmologyComovingBoxSize": float,
                 "CosmologyOmegaMatterNow": float,
                 "CosmologyHubbleConstantNow": float,
                 "CosmologyInitialRedshift": float,
                 "DualEnergyFormalismEta1": float,
                 "DualEnergyFormalismEta2": float,
                 "MetaDataString": str,
                 "HydroMethod": int,
                 "DualEnergyFormalism": int,
                 "InitialTime": float,
                 "ComovingCoordinates": int,
                 "DensityUnits": float,
                 "LengthUnits": float,
                 "TemperatureUnits": float,
                 "TimeUnits": float,
                 "GravitationalConstant": float,
                 "Gamma": float,
                 "MultiSpecies": int,
                 "CompilerPrecision": str,
                 "CurrentTimeIdentifier": int,
                 "RefineBy": int,
                 "BoundaryConditionName": str,
                 "TopGridRank": int,
                 "TopGridDimensions": int,
                }

mpc_conversion = {'mpc'   : 1e0,
                  'kpc'   : 1e3,
                  'pc'    : 1e6,
                  'au'    : 2.063e11,
                  'rsun'  : 2.2167e13,
                  'cm'    : 3.0857e24,
                  'miles' : 1.917e19}

axis_labels = [('y','z'),('x','z'),('x','y')]
