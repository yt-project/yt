"""
Fields specific to Streaming data



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
import yt.fields.universal_fields
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions

class StreamFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("code_mass/code_length**3", [])),
        ("number_density", ("1/code_length**3", [])),
        ("pressure", ("dyne/code_length**2", [])),
        ("thermal_energy", ("erg / g", [])),
        ("temperature", ("K", [])),
        ("x-velocity", ("code_length/code_time", [])),
        ("y-velocity", ("code_length/code_time", [])),
        ("z-velocity", ("code_length/code_time", [])),
        ("magnetic_field_x", ("gauss", [])),
        ("magnetic_field_y", ("gauss", [])),
        ("magnetic_field_z", ("gauss", [])),
        ("radiation_acceleration_x", ("code_length/code_time**2", [])),
        ("radiation_acceleration_y", ("code_length/code_time**2", [])),
        ("radiation_acceleration_z", ("code_length/code_time**2", [])),
    )

    known_particle_fields = (
        ("particle_position_x", ("code_length", [])),
        ("particle_position_y", ("code_length", [])),
        ("particle_position_z", ("code_length", [])),
        ("particle_velocity_x", ("code_length/code_time", [])),
        ("particle_velocity_y", ("code_length/code_time", [])),
        ("particle_velocity_z", ("code_length/code_time", [])),
        ("particle_index", ("", [])),
        ("particle_gas_density", ("code_mass/code_length**3", [])),
        ("particle_gas_temperature", ("K", [])),
        ("particle_mass", ("code_mass", [])),
        ("particle_position_x", ("code_length", [])),
        ("particle_position_y", ("code_length", [])),
        ("particle_position_z", ("code_length", [])),
        ("particle_index", ("", [])),
        ("particle_gas_density", ("code_mass/code_length**3", [])),
        ("particle_gas_temperature", ("K", [])),
        ("particle_mass", ("code_mass", [])),
    )
