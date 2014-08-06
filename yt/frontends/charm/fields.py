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

CharmFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = CharmFieldInfo.add_field

KnownCharmFields = FieldInfoContainer()
add_charm_field = KnownCharmFields.add_field

add_charm_field("potential", function=NullFunc, take_log=False,
                validators = [ValidateDataField("potential")],
                units=r"")

add_charm_field("density", function=NullFunc, take_log=False,
                validators = [ValidateDataField("density")],
                units=r"")

add_charm_field("gravitational_field_x", function=NullFunc, take_log=False,
                validators = [ValidateDataField("gravitational_field_x")],
                units=r"")

add_charm_field("gravitational_field_y", function=NullFunc, take_log=False,
                validators = [ValidateDataField("gravitational_field_y")],
                units=r"")

add_charm_field("gravitational_field_z", function=NullFunc, take_log=False,
                validators = [ValidateDataField("gravitational_field_z")],
                units=r"")

def _Density(field, data):
    return data["density"]
add_field("Density",function=_Density, take_log=True,
          units=r'\rm{g}/\rm{cm^3}')

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
                        "acceleration_x",
                        "acceleration_y",
                        "acceleration_z"]

for pf in _particle_field_list:
    pfunc = particle_func("%s" % (pf))
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

#do overrides for 2D

Charm2DFieldInfo = FieldInfoContainer.create_with_fallback(CharmFieldInfo)
add_charm_2d_field = Charm2DFieldInfo.add_field

def _gravitational_field_z(field, data):
    return np.zeros(data['gravitational_field_x'].shape,
                    dtype='float64')
add_charm_2d_field("gravitational_field_z", function=_gravitational_field_z)

def _particle_position_z(field, data):
    return np.zeros(data['particle_position_x'].shape, dtype='float64')
add_charm_2d_field("particle_position_z", function=_particle_position_z)

def _particle_velocity_z(field, data):
    return np.zeros(data['particle_velocity_x'].shape, dtype='float64')
add_charm_2d_field("particle_velocity_z", function=_particle_velocity_z)

def _particle_acceleration_z(field, data):
    return np.zeros(data['particle_acceleration_x'].shape, dtype='float64')
add_charm_2d_field("particle_acceleration_z", function=_particle_acceleration_z)

#do overrides for 1D

Charm1DFieldInfo = FieldInfoContainer.create_with_fallback(CharmFieldInfo)
add_charm_1d_field = Charm1DFieldInfo.add_field

def _gravitational_field_y(field, data):
    return np.zeros(data['gravitational_field_y'].shape,
                    dtype='float64')

def _particle_position_y(field, data):
    return np.zeros(data['particle_position_x'].shape, dtype='float64')

def _particle_velocity_y(field, data):
    return np.zeros(data['particle_velocity_x'].shape, dtype='float64')

def _particle_acceleration_y(field, data):
    return np.zeros(data['particle_acceleration_x'].shape, dtype='float64')

add_charm_1d_field("gravitational_field_z", function=_gravitational_field_z)
add_charm_1d_field("gravitational_field_y", function=_gravitational_field_y)

add_charm_1d_field("particle_position_z", function=_particle_position_z)
add_charm_1d_field("particle_velocity_z", function=_particle_velocity_z)
add_charm_1d_field("particle_acceleration_z", function=_particle_acceleration_z)

add_charm_1d_field("particle_position_y", function=_particle_position_y)
add_charm_1d_field("particle_velocity_y", function=_particle_velocity_y)
add_charm_1d_field("particle_acceleration_y", function=_particle_acceleration_y)