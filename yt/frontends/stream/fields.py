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

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.fields.universal_fields
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions

StreamFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo, "SFI")
add_field = StreamFieldInfo.add_field

KnownStreamFields = FieldInfoContainer()
add_stream_field = KnownStreamFields.add_field

add_stream_field("density", function = NullFunc, units='code_mass/code_length**3')
add_stream_field("number_density", function = NullFunc, units='1/code_length**3')
add_stream_field("pressure", function = NullFunc, units='dyne/code_length**2')
add_stream_field("temperature", function = NullFunc, units='K')
add_stream_field("x-velocity", function = NullFunc, units='code_length/code_time')
add_stream_field("y-velocity", function = NullFunc, units='code_length/code_time')
add_stream_field("z-velocity", function = NullFunc, units='code_length/code_time')
add_stream_field("magnetic_field_x", function = NullFunc, units='gauss')
add_stream_field("magnetic_field_y", function = NullFunc, units='gauss')
add_stream_field("magnetic_field_z", function = NullFunc, units='gauss')
add_stream_field("radiation_acceleration_x", function = NullFunc, units='code_length/code_time**2')
add_stream_field("radiation_acceleration_y", function = NullFunc, units='code_length/code_time**2')
add_stream_field("radiation_acceleration_z", function = NullFunc, units='code_length/code_time**2')

add_stream_field(("all", "particle_position_x"), function = NullFunc, particle_type=True,
          units='code_length')
add_stream_field(("all", "particle_position_y"), function = NullFunc, particle_type=True,
          units='code_length')
add_stream_field(("all", "particle_position_z"), function = NullFunc, particle_type=True,
          units='code_length')
add_stream_field(("all", "particle_velocity_x"), function = NullFunc, particle_type=True,
          units='code_length/code_time')
add_stream_field(("all", "particle_velocity_y"), function = NullFunc, particle_type=True,
          units='code_length/code_time')
add_stream_field(("all", "particle_velocity_z"), function = NullFunc, particle_type=True,
          units='code_length/code_time')
add_stream_field(("all", "particle_index"), function = NullFunc, particle_type=True,
          units='')
add_stream_field(("all", "particle_gas_density"), function = NullFunc, particle_type=True,
          units='code_mass/code_length**3')
add_stream_field(("all", "particle_gas_temperature"), function = NullFunc, particle_type=True,
          units='K')
add_stream_field(("all", "particle_mass"), function = NullFunc, particle_type=True, units='code_mass')
add_stream_field(("all", "particle_position_x"), function = NullFunc,
          particle_type=True, units='code_length')
add_stream_field(("all", "particle_position_y"), function = NullFunc,
          particle_type=True, units='code_length')
add_stream_field(("all", "particle_position_z"), function = NullFunc,
          particle_type=True, units='code_length')
add_stream_field(("all", "particle_index"), function = NullFunc, particle_type=True,
          units='')
add_stream_field(("all", "particle_gas_density"), function = NullFunc,
          particle_type=True, units='code_mass/code_length**3')
add_stream_field(("all", "particle_gas_temperature"), function = NullFunc,
          particle_type=True, units='K')
add_stream_field(("all", "particle_mass"), function = NullFunc, particle_type=True,
          units='code_mass')

add_stream_field("dark_matter_density", function = NullFunc, units='code_mass/code_length**3')
add_stream_field("star_density", function = NullFunc, units='code_mass/code_length**3')

particle_vector_functions("all",
        ["particle_position_%s" % ax for ax in 'xyz'],
        ["particle_velocity_%s" % ax for ax in 'xyz'],
        StreamFieldInfo)
particle_deposition_functions("all", "Coordinates",
    "particle_mass", StreamFieldInfo)
