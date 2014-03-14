"""
ART-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.units.yt_array import \
    YTArray
from yt.frontends.art.definitions import *

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
# NOTE: ARTIO uses momentum density.
mom_units = "code_mass / (code_length**2 * code_time)"
en_units = "code_mass*code_velocity**2/code_length**3"

class ARTFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Density", (rho_units, ["density"], None)),
        ("TotalEnergy", (en_units, ["total_energy"], None)),
        ("XMomentumDensity", (mom_units, ["momentum_x"], None)),
        ("YMomentumDensity", (mom_units, ["momentum_y"], None)),
        ("ZMomentumDensity", (mom_units, ["momentum_z"], None)),
        ("Pressure", ("", ["pressure"], None)), # Unused
        ("Gamma", ("", "gamma", None)),
        ("GasEnergy", (en_units, ["thermal_energy"], None)),
        ("MetalDensitySNII", (rho_units, ["metal_ii_density"], None)),
        ("MetalDensitySNIa", (rho_units, ["metal_ia_density"], None)),
        ("PotentialNew", ("", ["potential"], None)),
        ("PotentialOld", ("", ["gas_potential"], None)),
    )

    known_particle_fields = (
        ("POSITION_X", ("code_length", ["particle_position_x"], None)),
        ("POSITION_Y", ("code_length", ["particle_position_y"], None)),
        ("POSITION_Z", ("code_length", ["particle_position_z"], None)),
        ("VELOCITY_X", (vel_units, ["particle_velocity_x"], None)),
        ("VELOCITY_Y", (vel_units, ["particle_velocity_y"], None)),
        ("VELOCITY_Z", (vel_units, ["particle_velocity_z"], None)),
        ("MASS", ("code_mass", ["particle_mass"], None)),
        ("PID", ("", ["particle_index"], None)),
        ("SPECIES", ("", ["particle_type"], None)),
        ("BIRTH_TIME", ("code_time", ["creation_time"], None)),
        ("INITIAL_MASS", ("code_mass", ["initial_mass"], None)),
        ("METALLICITY_SNIa", ("", ["metallicity_snia"], None)),
        ("METALLICITY_SNII", ("", ["metallicity_snii"], None)),
    )
