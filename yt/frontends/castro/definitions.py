"""
Various definitions for various other modules and routines


Authors:
 * J. S. Oishi 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
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
