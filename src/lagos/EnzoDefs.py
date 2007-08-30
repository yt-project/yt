"""
Various definitions for various other modules and routines

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Move into yt.Defs, along with enki.EnkiDefs
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

# The number of levels we expect to have at most
MAXLEVEL=48

# Number of levels to toss back to check for maximum
NUMTOCHECK=2

axis_labels = [('y','z'),('x','z'),('x','y')]
axis_names = {0: 'x', 1: 'y', 2: 'z'}

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
                 "BoundaryConditionName": str
                }

"""
@todo: Add miles!
"""
unitList = {'mpc'   : 1e0,
            'kpc'   : 1e3,
            'pc'    : 1e6,
            'au'    : 2.063e11,
            'rsun'  : 2.2167e13,
            'cm'    : 3.0857e24 }

axis_labels = [('y','z'),('x','z'),('x','y')]

rates_out_key = \
      [ "tgas", \
        "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11",  \
        "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k21", "k22",\
        "k23", "k50", "k51", "k52", "k53", "k54", "k55", "k56", \
        "k13_1", "k13_2", "k13_3", "k13_4", "k13_5", "k13_6", "k13_7"         ]

cool_out_key = \
      [ "tgas", \
        "ceHI", "ceHeI", "ceHeII", "ciHI", "ciHeI", "ciHeIS", "ciHeII",\
        "reHII", "reHeII1", "reHeII2", "reHeIII", "brem", "comp", \
        "gphdl", "gpldl", "cieco", "vibh", "hyd01k", "rotl", "roth", \
        "h2k01", "hdlte", "hdlow", "hdc_1", "hdc_2", "hdc_3", "hdc_4", \
        "hdc_5"]

# Note that we default to CellMass as the weight

field_weights = \
    {"CellMass" : -999}
