"""
Charm-specific fields



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
    FieldInfo, \
    NullFunc, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
import numpy as np

KnownCharmFields = FieldInfoContainer()
add_charm_field = KnownCharmFields.add_field

CharmFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = CharmFieldInfo.add_field

add_charm_field("potential", function=NullFunc, take_log=False,
                validators = [ValidateDataField("potential")],
                units=r"")

def particle_func(p_field, dtype='float64'):
    def _Particles(field, data):
        io = data.hierarchy.io
        if not data.NumberOfParticles > 0:
            return np.array([], dtype=dtype)
        else:
            return io._read_particles(data, p_field).astype(dtype)
        
    return _Particles

_particle_field_list = ["mass",
                        "position_x",
                        "position_y",
                        "position_z",
                        "velocity_x",
                        "velocity_y",
                        "velocity_z",
                        "particle_id"]

for pf in _particle_field_list:
    pfunc = particle_func("particle_%s" % (pf))
    add_field("particle_%s" % pf, function=pfunc,
              validators = [ValidateSpatial(0)],
              particle_type=True)

def _ParticleMass(field, data):
    particles = data["particle_mass"].astype('float64')
    return particles

def _ParticleMassMsun(field, data):
    particles = data["particle_mass"].astype('float64')
    return particles/1.989e33

add_field("ParticleMass",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True)
add_field("ParticleMassMsun",
          function=_ParticleMassMsun, validators=[ValidateSpatial(0)],
          particle_type=True)
