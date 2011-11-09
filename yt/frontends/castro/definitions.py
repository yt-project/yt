"""
Various definitions for various other modules and routines

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2010 J.S. Oishi.  All Rights Reserved.

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
from yt.funcs import *

def boxlib_bool_to_int(v):
    try:
        return int(v)
    except ValueError:
        pass
    v = v.upper().strip()
    if v[0] == 'T':
        return 1
    elif v[0] == 'F':
        return 0

# TODO: get rid of enzo parameters we do not need
parameterDict = {"CosmologyCurrentRedshift": float,
                 "CosmologyComovingBoxSize": float,
                 "CosmologyOmegaMatterNow": float,
                 "CosmologyOmegaLambdaNow": float,
                 "CosmologyHubbleConstantNow": float,
                 "CosmologyInitialRedshift": float,
                 "DualEnergyFormalismEta1": float,
                 "DualEnergyFormalismEta2": float,
                 "MetaDataString": str,
                 "HydroMethod": int,
                 "DualEnergyFormalism": int,
                 "InitialTime": float,
                 "ComovingCoordinates": boxlib_bool_to_int,
                 "DensityUnits": float,
                 "LengthUnits": float,
                 "LengthUnit": float,
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
                 "EOSSoundSpeed": float,
                 "EOSType": int,
                 "NumberOfParticleAttributes": int,
                }

# converts the Castro inputs file name to the Enzo/yt name expected
# throughout the code. key is Castro name, value is Enzo/yt equivalent
castro2enzoDict = {"amr.n_cell": "TopGridDimensions",
                  "materials.gamma": "Gamma",
                  "amr.ref_ratio": "RefineBy",
                  "castro.use_comoving": "ComovingCoordinates",
                  "castro.redshift_in": "CosmologyInitialRedshift",
                  "comoving_OmL": "CosmologyOmegaLambdaNow",
                  "comoving_OmM": "CosmologyOmegaMatterNow",
                  "comoving_h": "CosmologyHubbleConstantNow"
                  }

yt2castroFieldsDict = {}
castro2ytFieldsDict = {}

castro_FAB_header_pattern = r"^FAB \(\((\d+), \([0-9 ]+\)\),\(\d+, \(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n"

castro_particle_field_names = \
    ['particle_position_%s' % ax for ax in 'xyz'] + \
    ['particle_mass'] +  \
    ['particle_velocity_%s' % ax for ax in 'xyz']
