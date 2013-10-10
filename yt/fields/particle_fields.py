"""
These are common particle deposition fields.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import *
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    NullFunc, \
    TranslationFunc
from yt.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_sun_cgs, \
    mh
from .universal_fields import \
    get_radius, \
    _get_conv

def _field_concat(fname):
    def _AllFields(field, data):
        v = []
        for ptype in data.pf.particle_types:
            data.pf._last_freq = (ptype, None)
            if ptype == "all" or \
                ptype in data.pf.known_filters:
                  continue
            v.append(data[ptype, fname].copy())
        rv = np.concatenate(v, axis=0)
        return rv
    return _AllFields

def _field_concat_slice(fname, axi):
    def _AllFields(field, data):
        v = []
        for ptype in data.pf.particle_types:
            data.pf._last_freq = (ptype, None)
            if ptype == "all" or \
                ptype in data.pf.known_filters:
                  continue
            v.append(data[ptype, fname][:,axi])
        rv = np.concatenate(v, axis=0)
        return rv
    return _AllFields

def particle_deposition_functions(ptype, coord_name, mass_name, registry):
    orig = set(registry.keys())
    def particle_count(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, method = "count")
        return d

    registry.add_field(("deposit", "%s_count" % ptype),
             function = particle_count,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Count}" % ptype,
             projection_conversion = '1')

    def particle_mass(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        return d

    registry.add_field(("deposit", "%s_mass" % ptype),
             function = particle_mass,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Mass}" % ptype,
             units = r"\mathrm{g}",
             projected_units = r"\mathrm{g}\/\mathrm{cm}",
             projection_conversion = 'cm')

    def particle_density(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        d /= data["gas","CellVolume"]
        return d

    registry.add_field(("deposit", "%s_density" % ptype),
             function = particle_density,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    def particle_cic(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
        d /= data["gas","CellVolume"]
        return d

    registry.add_field(("deposit", "%s_cic" % ptype),
             function = particle_cic,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s CIC Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    # Now some translation functions.

    def particle_ones(field, data):
        return np.ones(data[ptype, mass_name].shape, dtype="float64")

    registry.add_field((ptype, "particle_ones"),
                       function = particle_ones,
                       particle_type = True,
                       units = "")

    registry.add_field((ptype, "ParticleMass"),
            function = TranslationFunc((ptype, mass_name)),
            particle_type = True,
            units = r"\mathrm{g}")

    def _ParticleMassMsun(field, data):
        return data[ptype, mass_name].copy()
    def _conv_Msun(data):
        return 1.0/mass_sun_cgs

    registry.add_field((ptype, "ParticleMassMsun"),
            function = _ParticleMassMsun,
            convert_function = _conv_Msun,
            particle_type = True,
            units = r"\mathrm{M}_\odot")

    def particle_mesh_ids(field, data):
        pos = data[ptype, coord_name]
        ids = np.zeros(pos.shape[0], dtype="float64") - 1
        # This is float64 in name only.  It will be properly cast inside the
        # deposit operation.
        #_ids = ids.view("float64")
        data.deposit(pos, [ids], method = "mesh_id")
        return ids
    registry.add_field((ptype, "mesh_id"),
            function = particle_mesh_ids,
            validators = [ValidateSpatial()],
            particle_type = True)

    return list(set(registry.keys()).difference(orig))

def particle_scalar_functions(ptype, coord_name, vel_name, registry):

    # Now we have to set up the various velocity and coordinate things.  In the
    # future, we'll actually invert this and use the 3-component items
    # elsewhere, and stop using these.
    
    # Note that we pass in _ptype here so that it's defined inside the closure.
    orig = set(registry.keys())

    def _get_coord_funcs(axi, _ptype):
        def _particle_velocity(field, data):
            return data[_ptype, vel_name][:,axi]
        def _particle_position(field, data):
            return data[_ptype, coord_name][:,axi]
        return _particle_velocity, _particle_position
    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs(axi, ptype)
        registry.add_field((ptype, "particle_velocity_%s" % ax),
            particle_type = True, function = v)
        registry.add_field((ptype, "particle_position_%s" % ax),
            particle_type = True, function = p)

    return list(set(registry.keys()).difference(orig))

def particle_vector_functions(ptype, coord_names, vel_names, registry):

    # This will column_stack a set of scalars to create vector fields.
    orig = set(registry.keys())

    def _get_vec_func(_ptype, names):
        def particle_vectors(field, data):
            return np.column_stack([data[_ptype, name] for name in names])
        return particle_vectors
    registry.add_field((ptype, "Coordinates"),
                       function=_get_vec_func(ptype, coord_names),
                       particle_type=True)
    registry.add_field((ptype, "Velocities"),
                       function=_get_vec_func(ptype, vel_names),
                       particle_type=True)
    return list(set(registry.keys()).difference(orig))

def standard_particle_fields(registry, ptype,
                             spos = "particle_position_%s",
                             svel = "particle_velocity_%s"):
    # This function will set things up based on the scalar fields and standard
    # yt field names.
    
    def _ParticleVelocityMagnitude(field, data):
        bulk_velocity = data.get_field_parameter("bulk_velocity")
        if bulk_velocity == None:
            bulk_velocity = np.zeros(3)
        return ( (data[ptype, svel % 'x']-bulk_velocity[0])**2.0 + \
                 (data[ptype, svel % 'y']-bulk_velocity[1])**2.0 + \
                 (data[ptype, svel % 'z']-bulk_velocity[2])**2.0 )**(1.0/2.0)
    registry.add_field((ptype, "ParticleVelocityMagnitude"),
              function=_ParticleVelocityMagnitude,
              particle_type=True, 
              take_log=False, units=r"\rm{cm}/\rm{s}")

    # Angular momentum
    def _ParticleSpecificAngularMomentumX(field, data):
        if data.has_field_parameter("bulk_velocity"):
            bv = data.get_field_parameter("bulk_velocity")
        else: bv = np.zeros(3, dtype='float64')
        center = data.get_field_parameter('center')
        y = data[ptype, spos % 'y'] - center[1]
        z = data[ptype, spos % 'z'] - center[2]
        yv = data[ptype, svel % 'y'] - bv[1]
        zv = data[ptype, svel % 'z'] - bv[2]
        return yv*z - zv*y
    registry.add_field((ptype, "ParticleSpecificAngularMomentumX"), 
              function=_ParticleSpecificAngularMomentumX,
              particle_type=True,
              convert_function=_get_conv("cm"),
              units=r"\rm{cm}^2/\rm{s}",
              validators=[ValidateParameter("center")])
    def _ParticleSpecificAngularMomentumY(field, data):
        if data.has_field_parameter("bulk_velocity"):
            bv = data.get_field_parameter("bulk_velocity")
        else: bv = np.zeros(3, dtype='float64')
        center = data.get_field_parameter('center')
        x = data[ptype, spos % 'x'] - center[0]
        z = data[ptype, spos % 'z'] - center[2]
        xv = data[ptype, svel % 'x'] - bv[0]
        zv = data[ptype, svel % 'z'] - bv[2]
        return -(xv*z - zv*x)
    registry.add_field((ptype, "ParticleSpecificAngularMomentumY"), 
              function=_ParticleSpecificAngularMomentumY,
              particle_type=True,
              convert_function=_get_conv("cm"),
              units=r"\rm{cm}^2/\rm{s}",
              validators=[ValidateParameter("center")])
    def _ParticleSpecificAngularMomentumZ(field, data):
        if data.has_field_parameter("bulk_velocity"):
            bv = data.get_field_parameter("bulk_velocity")
        else: bv = np.zeros(3, dtype='float64')
        center = data.get_field_parameter('center')
        x = data[ptype, spos % 'x'] - center[0]
        y = data[ptype, spos % 'y'] - center[1]
        xv = data[ptype, svel % 'x'] - bv[0]
        yv = data[ptype, svel % 'y'] - bv[1]
        return xv*y - yv*x
    registry.add_field((ptype, "ParticleSpecificAngularMomentumZ"), 
              function=_ParticleSpecificAngularMomentumZ,
              particle_type=True,
              convert_function=_get_conv("cm"),
              units=r"\rm{cm}^2/\rm{s}",
              validators=[ValidateParameter("center")])

    def _ParticleAngularMomentumX(field, data):
        return data[ptype, "particle_mass"] * \
            data[ptype, "ParticleSpecificAngularMomentumX"]
    registry.add_field((ptype, "ParticleAngularMomentumX"),
             function=_ParticleAngularMomentumX,
             units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
             validators=[ValidateParameter('center')])
    def _ParticleAngularMomentumY(field, data):
        return data[ptype, "particle_mass"] * \
            data[ptype, "ParticleSpecificAngularMomentumY"]
    registry.add_field((ptype, "ParticleAngularMomentumY"),
             function=_ParticleAngularMomentumY,
             units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
             validators=[ValidateParameter('center')])
    def _ParticleAngularMomentumZ(field, data):
        return data[ptype, "particle_mass"] * \
            data[ptype, "ParticleSpecificAngularMomentumZ"]
    registry.add_field((ptype, "ParticleAngularMomentumZ"),
             function=_ParticleAngularMomentumZ,
             units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
             validators=[ValidateParameter('center')])

    def _ParticleRadius(field, data):
        return get_radius(data, "particle_position_")

    registry.add_field((ptype, "ParticleRadius"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("cm"), units=r"\rm{cm}",
              particle_type = True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiusMpc"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("mpc"), units=r"\rm{Mpc}",
              particle_type=True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiuskpc"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("kpc"), units=r"\rm{kpc}",
              particle_type=True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiuskpch"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("kpch"), units=r"\rm{kpc}/\rm{h}",
              particle_type=True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiuspc"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("pc"), units=r"\rm{pc}",
              particle_type=True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiusAU"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              convert_function = _get_conv("au"), units=r"\rm{AU}",
              particle_type=True,
              display_name = "Particle Radius")
    registry.add_field((ptype, "ParticleRadiusCode"),
              function=_ParticleRadius,
              validators=[ValidateParameter("center")],
              particle_type=True,
              display_name = "Particle Radius (code)")
    def _ParticleRadiusSpherical(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = spos
        pos = np.array([data[ptype, pos % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        sphr = get_sph_r_component(pos, theta, phi, normal)
        return sphr

    registry.add_field((ptype, "ParticleRadiusSpherical"),
              function=_ParticleRadiusSpherical,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

    def _ParticleThetaSpherical(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = spos
        pos = np.array([data[ptype, pos % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        spht = get_sph_theta_component(pos, theta, phi, normal)
        return spht

    registry.add_field((ptype, "ParticleThetaSpherical"),
              function=_ParticleThetaSpherical,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

    def _ParticlePhiSpherical(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = spos
        pos = np.array([data[ptype, pos % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        vel = vel - np.reshape(bv, (3, 1))
        sphp = get_sph_phi_component(pos, theta, phi, normal)
        return sphp

    registry.add_field((ptype, "ParticlePhiSpherical"),
              function=_ParticleThetaSpherical,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

    def _ParticleRadialVelocity(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = spos
        pos = np.array([data[ptype, pos % ax] for ax in "xyz"])
        vel = svel
        vel = np.array([data[ptype, vel % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        vel = vel - np.reshape(bv, (3, 1))
        sphr = get_sph_r_component(vel, theta, phi, normal)
        return sphr

    registry.add_field((ptype, "ParticleRadialVelocity"),
              function=_ParticleRadialVelocity,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

    def _ParticleThetaVelocity(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = spos
        pos = np.array([data[ptype, pos % ax] for ax in "xyz"])
        vel = svel
        vel = np.array([data[ptype, vel % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        vel = vel - np.reshape(bv, (3, 1))
        spht = get_sph_theta_component(vel, theta, phi, normal)
        return spht

    registry.add_field((ptype, "ParticleThetaVelocity"),
              function=_ParticleThetaVelocity,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

    def _ParticlePhiVelocity(field, data):
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        pos = np.array([data[ptype, spos % ax] for ax in "xyz"])
        vel = np.array([data[ptype, svel % ax] for ax in "xyz"])
        theta = get_sph_theta(pos, center)
        phi = get_sph_phi(pos, center)
        pos = pos - np.reshape(center, (3, 1))
        vel = vel - np.reshape(bv, (3, 1))
        sphp = get_sph_phi_component(vel, phi, normal)
        return sphp

    registry.add_field((ptype, "ParticlePhiVelocity"),
              function=_ParticleThetaVelocity,
              particle_type=True, units=r"\rm{cm}/\rm{s}",
              validators=[ValidateParameter("normal"), 
                          ValidateParameter("center")])

