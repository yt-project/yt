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

def _setup_particle_fields(registry, ptype):
    particle_vector_functions(ptype,
        ["particle_position_%s" % ax for ax in 'xyz'],
        ["particle_velocity_%s" % ax for ax in 'xyz'],
        registry)
    particle_deposition_functions(ptype,
        "Coordinates", "particle_mass", registry)
    for fn in ["creation_time", "dynamical_time", "metallicity_fraction"] + \
              ["particle_type", "particle_index", "ParticleMass"] + \
              ["particle_position_%s" % ax for ax in 'xyz']:
        registry.add_field((ptype, fn), function=NullFunc, particle_type=True)
