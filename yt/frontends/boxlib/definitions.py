"""
Various definitions for various other modules and routines



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


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
                 "ComovingCoordinates": int,
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


# converts the Orion inputs file name to the Enzo/yt name expected
# throughout the code. key is Orion name, value is Enzo/yt equivalent
orion2enzoDict = {"amr.n_cell": "TopGridDimensions",
                  "materials.gamma": "Gamma",
                  "amr.ref_ratio": "RefineBy",
                  "castro.use_comoving": "ComovingCoordinates",
                  "castro.redshift_in": "CosmologyInitialRedshift",
                  "comoving_OmL": "CosmologyOmegaLambdaNow",
                  "comoving_OmM": "CosmologyOmegaMatterNow",
                  "comoving_h": "CosmologyHubbleConstantNow"
                  }

