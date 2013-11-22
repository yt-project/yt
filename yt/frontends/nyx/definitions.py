"""
Definitions specific to Nyx



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from utils import boxlib_bool_to_int

# This gives the type of each parameter we want, used to cast with a ``map()``
# @todo: get rid of enzo parameters we do not need
parameter_type_dict = {
    "CosmologyCurrentRedshift": float,
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

# Provides translation between parameters in the nyx `inputs` file names to the
# enzo/yt name expected throughout the code. The key is nyx name, value is
# enzo/yt equivalent.
nyx_to_enzo_dict = {
    "amr.n_cell": "TopGridDimensions",
    "amr.ref_ratio": "RefineBy",
    "materials.gamma": "Gamma",
    "castro.use_comoving": "ComovingCoordinates",
    "castro.redshift_in": "CosmologyInitialRedshift",
    "comoving_OmL": "CosmologyOmegaLambdaNow",
    "comoving_OmM": "CosmologyOmegaMatterNow",
    "comoving_h": "CosmologyHubbleConstantNow"
}

# is this the same as nyx_to_enzo.
yt_to_nyx_fields_dict = {}
nyx_to_yt_fields_dict = {}

fab_header_pattern = r"^FAB \(\((\d+), \([0-9 ]+\)\),\(\d+, \(([0-9 ]+)\)\)\)\(\((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\) \((\d+,\d+,\d+)\)\) (\d+)\n"

# need to specify units eventually
nyx_particle_field_names = ['particle_position_%s' % ax for ax in 'xyz'] + \
                           ['particle_mass'] +  \
                           ['particle_velocity_%s' % ax for ax in 'xyz']

