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

KnownStreamFields = FieldInfoContainer()
add_stream_field = KnownStreamFields.add_field

StreamFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = StreamFieldInfo.add_field

add_stream_field("density", function = NullFunc)
add_stream_field("x-velocity", function = NullFunc)
add_stream_field("y-velocity", function = NullFunc)
add_stream_field("z-velocity", function = NullFunc)

add_field("Density", function = TranslationFunc("density"))

add_stream_field("particle_position_x", function = NullFunc, particle_type=True)
add_stream_field("particle_position_y", function = NullFunc, particle_type=True)
add_stream_field("particle_position_z", function = NullFunc, particle_type=True)
add_stream_field("particle_index", function = NullFunc, particle_type=True)
add_stream_field("particle_gas_density", function = NullFunc, particle_type=True)
add_stream_field("particle_gas_temperature", function = NullFunc, particle_type=True)
add_stream_field("particle_mass", function = NullFunc, particle_type=True)

add_field("ParticleMass", function = TranslationFunc("particle_mass"),
          particle_type=True)

add_stream_field(("all", "particle_position_x"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_position_y"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_position_z"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_index"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_gas_density"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_gas_temperature"), function = NullFunc, particle_type=True)
add_stream_field(("all", "particle_mass"), function = NullFunc, particle_type=True)

add_field(("all", "ParticleMass"), function = TranslationFunc("particle_mass"),
          particle_type=True)

particle_vector_functions("all", ["particle_position_%s" % ax for ax in 'xyz'],
                                 ["particle_velocity_%s" % ax for ax in 'xyz'],
                          StreamFieldInfo)
particle_deposition_functions("all", "Coordinates", "ParticleMass",
                               StreamFieldInfo)
