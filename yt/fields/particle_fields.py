"""
These are common particle fields.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.fields.derived_field import \
    ValidateParameter, \
    ValidateSpatial

from yt.fields.field_detector import \
    FieldDetector

from yt.funcs import \
    issue_deprecation_warning

from yt.units.yt_array import \
    uconcatenate, \
    ucross

from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component, \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, \
    get_sph_theta, get_sph_phi, \
    modify_reference_frame

from .vector_operations import \
    create_magnitude_field

from .field_functions import \
    get_radius

sph_whitelist_fields = (
    'density',
    'temperature',
    'metallicity',
    'thermal_energy',
    'smoothing_length',
    'H_fraction',
    'He_fraction',
    'C_fraction',
    'N_fraction',
    'O_fraction',
    'Ne_fraction',
    'Mg_fraction',
    'Si_fraction',
    'Fe_fraction',
    'C_density',
    'O_density',
    'Si_density',
    'Fe_density'
)


def _field_concat(fname):
    def _AllFields(field, data):
        v = []
        for ptype in data.ds.particle_types:
            data.ds._last_freq = (ptype, None)
            if ptype == "all" or \
                ptype in data.ds.known_filters:
                  continue
            v.append(data[ptype, fname].copy())
        rv = uconcatenate(v, axis=0)
        return rv
    return _AllFields

def _field_concat_slice(fname, axi):
    def _AllFields(field, data):
        v = []
        for ptype in data.ds.particle_types:
            data.ds._last_freq = (ptype, None)
            if ptype == "all" or \
                ptype in data.ds.known_filters:
                  continue
            v.append(data[ptype, fname][:,axi])
        rv = uconcatenate(v, axis=0)
        return rv
    return _AllFields

def particle_deposition_functions(ptype, coord_name, mass_name, registry):
    unit_system = registry.ds.unit_system
    orig = set(registry.keys())
    ptype_dn = ptype.replace("_"," ").title()
    def particle_count(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, method = "count")
        d = data.ds.arr(d, input_units = "cm**-3")
        return data.apply_units(d, field.units)

    registry.add_field(("deposit", "%s_count" % ptype),
                       sampling_type="cell",
                       function = particle_count,
                       validators = [ValidateSpatial()],
                       units = '',
                       display_name = r"\mathrm{%s Count}" % ptype_dn)

    def particle_mass(field, data):
        pos = data[ptype, coord_name]
        pmass = data[ptype, mass_name]
        pmass.convert_to_units(field.units)
        d = data.deposit(pos, [pmass], method = "sum")
        return data.apply_units(d, field.units)

    registry.add_field(("deposit", "%s_mass" % ptype),
                       sampling_type="cell",
                       function = particle_mass,
                       validators = [ValidateSpatial()],
                       display_name = r"\mathrm{%s Mass}" % ptype_dn,
                       units = unit_system["mass"])

    def particle_density(field, data):
        pos = data[ptype, coord_name].convert_to_units("code_length")
        mass = data[ptype, mass_name].convert_to_units("code_mass")
        d = data.deposit(pos, [mass], method = "sum")
        d = data.ds.arr(d, "code_mass")
        d /= data["index", "cell_volume"]
        return d

    registry.add_field(("deposit", "%s_density" % ptype),
                       sampling_type="cell",
                       function = particle_density,
                       validators = [ValidateSpatial()],
                       display_name = r"\mathrm{%s Density}" % ptype_dn,
                       units = unit_system["density"])

    def particle_cic(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
        d = data.apply_units(d, data[ptype, mass_name].units)
        d /= data["index", "cell_volume"]
        return d

    registry.add_field(("deposit", "%s_cic" % ptype),
                       sampling_type="cell",
                       function = particle_cic,
                       validators = [ValidateSpatial()],
                       display_name = r"\mathrm{%s CIC Density}" % ptype_dn,
                       units = unit_system["density"])

    def _get_density_weighted_deposit_field(fname, units, method):
        def _deposit_field(field, data):
            """
            Create a grid field for particle quantities weighted by particle
            mass, using cloud-in-cell deposit.
            """
            pos = data[ptype, "particle_position"]
            # Get back into density
            pden = data[ptype, 'particle_mass']
            top = data.deposit(pos, [data[(ptype, fname)]*pden], method=method)
            bottom = data.deposit(pos, [pden], method=method)
            top[bottom == 0] = 0.0
            bnz = bottom.nonzero()
            top[bnz] /= bottom[bnz]
            d = data.ds.arr(top, input_units=units)
            return d
        return _deposit_field

    for ax in 'xyz':
        for method, name in zip(("cic", "sum"), ("cic", "nn")):
            function = _get_density_weighted_deposit_field(
                "particle_velocity_%s" % ax, "code_velocity", method)
            registry.add_field(
                ("deposit", ("%s_"+name+"_velocity_%s") % (ptype, ax)),
                sampling_type="cell",
                function=function,
                units=unit_system["velocity"],
                take_log=False,
                validators=[ValidateSpatial(0)])

    for method, name in zip(("cic", "sum"), ("cic", "nn")):
        function = _get_density_weighted_deposit_field(
            "age", "code_time", method)
        registry.add_field(
            ("deposit", ("%s_"+name+"_age") % (ptype)),
            sampling_type="cell",
            function=function,
            units=unit_system["time"],
            take_log=False,
            validators=[ValidateSpatial(0)])

    # Now some translation functions.

    def particle_ones(field, data):
        v = np.ones(data[ptype, coord_name].shape[0], dtype="float64")
        return data.apply_units(v, field.units)

    registry.add_field((ptype, "particle_ones"),
                       sampling_type="particle",
                       function = particle_ones,
                       units = "",
                       display_name = r"Particle Count")

    def particle_mesh_ids(field, data):
        pos = data[ptype, coord_name]
        ids = np.zeros(pos.shape[0], dtype="float64") - 1
        # This is float64 in name only.  It will be properly cast inside the
        # deposit operation.
        #_ids = ids.view("float64")
        data.deposit(pos, [ids], method = "mesh_id")
        return data.apply_units(ids, "")
    registry.add_field((ptype, "mesh_id"),
                       sampling_type="particle",
                       function = particle_mesh_ids,
                       validators = [ValidateSpatial()],
                       units = '')

    return list(set(registry.keys()).difference(orig))

def particle_scalar_functions(ptype, coord_name, vel_name, registry):

    # Now we have to set up the various velocity and coordinate things.  In the
    # future, we'll actually invert this and use the 3-component items
    # elsewhere, and stop using these.

    # Note that we pass in _ptype here so that it's defined inside the closure.

    def _get_coord_funcs(axi, _ptype):
        def _particle_velocity(field, data):
            return data[_ptype, vel_name][:,axi]
        def _particle_position(field, data):
            return data[_ptype, coord_name][:, axi]
        return _particle_velocity, _particle_position
    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs(axi, ptype)
        registry.add_field((ptype, "particle_velocity_%s" % ax),
                           sampling_type="particle",
                           function = v,
                           units = "code_velocity")
        registry.add_field((ptype, "particle_position_%s" % ax),
                           sampling_type="particle",
                           function = p,
                           units = "code_length")

def particle_vector_functions(ptype, coord_names, vel_names, registry):

    unit_system = registry.ds.unit_system

    # This will column_stack a set of scalars to create vector fields.

    def _get_vec_func(_ptype, names):
        def particle_vectors(field, data):
            v = [data[_ptype, name].in_units(field.units)
                  for name in names]
            c = np.column_stack(v)
            return data.apply_units(c, field.units)
        return particle_vectors

    registry.add_field((ptype, "particle_position"),
                       sampling_type="particle",
                       function=_get_vec_func(ptype, coord_names),
                       units = "code_length")

    registry.add_field((ptype, "particle_velocity"),
                       sampling_type="particle",
                       function=_get_vec_func(ptype, vel_names),
                       units = unit_system["velocity"])

def get_angular_momentum_components(ptype, data, spos, svel):
    if data.has_field_parameter("normal"):
        normal = data.get_field_parameter("normal")
    else:
        normal = data.ds.arr([0.0,0.0,1.0],"code_length") # default to simulation axis
    pos = data.ds.arr([data[ptype, spos % ax] for ax in "xyz"]).T
    vel = data.ds.arr(
        [data[ptype, 'relative_%s' % (svel % ax)] for ax in "xyz"]).T
    return pos, vel, normal

def standard_particle_fields(registry, ptype,
                             spos = "particle_position_%s",
                             svel = "particle_velocity_%s"):
    unit_system = registry.ds.unit_system

    def _particle_velocity_magnitude(field, data):
        """ M{|v|} """
        return np.sqrt(data[ptype, 'relative_%s' % (svel % 'x')]**2 +
                       data[ptype, 'relative_%s' % (svel % 'y')]**2 +
                       data[ptype, 'relative_%s' % (svel % 'z')]**2)

    registry.add_field((ptype, "particle_velocity_magnitude"),
                       sampling_type="particle",
                       function=_particle_velocity_magnitude,
                       take_log=False,
                       units=unit_system["velocity"])

    def _particle_specific_angular_momentum(field, data):
        """Calculate the angular of a particle velocity.

        Returns a vector for each particle.
        """
        center = data.get_field_parameter('center')
        pos, vel, normal = get_angular_momentum_components(
            ptype, data, spos, svel)
        L, r_vec, v_vec = modify_reference_frame(center, normal, P=pos, V=vel)
        # adding in the unit registry allows us to have a reference to the
        # dataset and thus we will always get the correct units after applying
        # the cross product.
        return -ucross(r_vec, v_vec, registry=data.ds.unit_registry)


    registry.add_field((ptype, "particle_specific_angular_momentum"),
                       sampling_type="particle",
                       function=_particle_specific_angular_momentum,
                       units=unit_system["specific_angular_momentum"],
                       validators=[ValidateParameter("center")])

    def _get_spec_ang_mom_comp(axi, ax, _ptype):
        def _particle_specific_angular_momentum_component(field, data):
            return data[_ptype, "particle_specific_angular_momentum"][:, axi]
        def _particle_angular_momentum_component(field, data):
            return data[_ptype, "particle_mass"] * \
                data[ptype, "particle_specific_angular_momentum_%s" % ax]
        return _particle_specific_angular_momentum_component, \
            _particle_angular_momentum_component
    for axi, ax in enumerate("xyz"):
        f, v = _get_spec_ang_mom_comp(axi, ax, ptype)
        registry.add_field(
            (ptype, "particle_specific_angular_momentum_%s" % ax),
            sampling_type="particle",
            function=f,
            units=unit_system["specific_angular_momentum"],
            validators=[ValidateParameter("center")]
        )
        registry.add_field(
            (ptype, "particle_angular_momentum_%s" % ax),
            sampling_type="particle",
            function=v, units=unit_system["angular_momentum"],
            validators=[ValidateParameter('center')])

    def _particle_angular_momentum(field, data):
        am = data[ptype, "particle_mass"] * data[ptype, "particle_specific_angular_momentum"].T
        return am.T

    registry.add_field((ptype, "particle_angular_momentum"),
                       sampling_type="particle",
                       function=_particle_angular_momentum,
                       units=unit_system["angular_momentum"],
                       validators=[ValidateParameter("center")])

    create_magnitude_field(registry, "particle_angular_momentum",
                           unit_system["angular_momentum"],
                           sampling_type='particle',
                           ftype=ptype)

    def _particle_radius(field, data):
        """The spherical radius component of the particle positions

        Relative to the coordinate system defined by the *normal* vector,
        and *center* field parameters.
        """
        return get_radius(data, "particle_position_", field.name[0])

    registry.add_field((ptype, "particle_radius"),
                       sampling_type="particle",
                       function=_particle_radius,
                       units=unit_system["length"],
                       validators=[ValidateParameter("center")])

    def _relative_particle_position(field, data):
        """The cartesian particle positions in a rotated reference frame

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.

        Note that the orientation of the x and y axes are arbitrary.
        """
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        pos = data.ds.arr([data[ptype, spos % ax] for ax in "xyz"]).T
        L, pos = modify_reference_frame(center, normal, P=pos)
        return pos

    def _particle_position_relative(field, data):
        if not isinstance(data, FieldDetector):
            issue_deprecation_warning(
                "The 'particle_position_relative' field has been deprecated in " +
                "favor of 'relative_particle_position'.")
        if isinstance(field.name, tuple):
            return data[field.name[0], 'relative_particle_position']
        else:
            return data['relative_particle_position']

    for name, func in zip(
            ["particle_position_relative", "relative_particle_position"],
            [_particle_position_relative,  _relative_particle_position]):
        registry.add_field(
            (ptype, name),
            sampling_type="particle",
            function=func,
            units=unit_system["length"],
            validators=[ValidateParameter("normal"),
                        ValidateParameter("center")])

    def _relative_particle_velocity(field, data):
        """The vector particle velocities in an arbitrary coordinate system

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.

        Note that the orientation of the x and y axes are arbitrary.
        """
        normal = data.get_field_parameter('normal')
        center = data.get_field_parameter('center')
        bv = data.get_field_parameter("bulk_velocity")
        vel = data.ds.arr(
            [data[ptype, svel % ax] - bv[iax]
             for iax, ax in enumerate("xyz")]).T
        L, vel = modify_reference_frame(center, normal, V=vel)
        return vel

    def _particle_velocity_relative(field, data):
        if not isinstance(data, FieldDetector):
            issue_deprecation_warning(
                "The 'particle_velocity_relative' field has been deprecated in " +
                "favor of 'relative_particle_velocity'.")
        if isinstance(field.name, tuple):
            return data[field.name[0], 'relative_particle_velocity']
        else:
            return data['relative_particle_velocity']

    for name, func in zip(
            ["particle_velocity_relative", "relative_particle_velocity"],
            [_particle_velocity_relative,  _relative_particle_velocity]):
        registry.add_field(
            (ptype, name),
            sampling_type="particle",
            function=func,
            units=unit_system["velocity"],
            validators=[ValidateParameter("normal"),
                        ValidateParameter("center")])

    def _get_coord_funcs_relative(axi, _ptype):
        def _particle_pos_rel(field, data):
            return data[_ptype, "relative_particle_position"][:, axi]
        def _particle_vel_rel(field, data):
            return data[_ptype, "relative_particle_velocity"][:, axi]
        return _particle_vel_rel, _particle_pos_rel
    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs_relative(axi, ptype)
        registry.add_field((ptype, "particle_velocity_relative_%s" % ax),
                           sampling_type="particle",
                           function=v,
                           units="code_velocity")
        registry.add_field((ptype, "particle_position_relative_%s" % ax),
            sampling_type="particle", function = p, units = "code_length")
        registry.add_field((ptype, "relative_particle_velocity_%s" % ax),
            sampling_type="particle", function = v, units = "code_velocity")
        registry.add_field((ptype, "relative_particle_position_%s" % ax),
            sampling_type="particle", function = p, units = "code_length")

    # this is just particle radius but we add it with an alias for the sake of
    # consistent naming
    registry.add_field((ptype, "particle_position_spherical_radius"),
                       sampling_type="particle",
                       function=_particle_radius,
                       units=unit_system["length"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_spherical_position_radius(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_position_spherical_radius']

    registry.add_field((ptype, "particle_spherical_position_radius"),
                       sampling_type="particle",
                       function=_particle_spherical_position_radius,
                       units=unit_system["length"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_position_spherical_theta(field, data):
        """The spherical theta coordinate of the particle positions.

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        return data.ds.arr(get_sph_theta(pos, normal), "")

    registry.add_field((ptype, "particle_position_spherical_theta"),
                       sampling_type="particle",
                       function=_particle_position_spherical_theta,
                       units="",
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")])

    def _particle_spherical_position_theta(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_position_spherical_theta']

    registry.add_field((ptype, "particle_spherical_position_theta"),
                       sampling_type="particle",
                       function=_particle_spherical_position_theta,
                       units="",
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_position_spherical_phi(field, data):
        """The spherical phi component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        return data.ds.arr(get_sph_phi(pos, normal), "")

    registry.add_field((ptype, "particle_position_spherical_phi"),
                       sampling_type="particle",
                       function=_particle_position_spherical_phi,
                       units="",
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_spherical_position_phi(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_position_spherical_phi']

    registry.add_field((ptype, "particle_spherical_position_phi"),
                       sampling_type="particle",
                       function=_particle_spherical_position_phi,
                       units="",
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")])

    def _particle_velocity_spherical_radius(field, data):
        """The spherical radius component of the particle velocities in an
         arbitrary coordinate system

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        vel = data['relative_particle_velocity'].T
        theta = get_sph_theta(pos, normal)
        phi = get_sph_phi(pos, normal)
        sphr = get_sph_r_component(vel, theta, phi, normal)
        return sphr

    registry.add_field((ptype, "particle_velocity_spherical_radius"),
                       sampling_type="particle",
                       function=_particle_velocity_spherical_radius,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_spherical_velocity_radius(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_velocity_spherical_radius']

    registry.add_field((ptype, "particle_spherical_velocity_radius"),
                       sampling_type="particle",
                       function=_particle_spherical_velocity_radius,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    # particel_velocity_spherical_radius is simply aliased to
    # "particle_radial_velocity" for convenience
    registry.add_field((ptype, "particle_radial_velocity"),
                       sampling_type="particle",
                       function=_particle_spherical_velocity_radius,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_velocity_spherical_theta(field, data):
        """The spherical theta component of the particle velocities in an
         arbitrary coordinate system

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        vel = data['relative_particle_velocity'].T
        theta = get_sph_theta(pos, normal)
        phi = get_sph_phi(pos, normal)
        spht = get_sph_theta_component(vel, theta, phi, normal)
        return spht

    registry.add_field((ptype, "particle_velocity_spherical_theta"),
                       sampling_type="particle",
                       function=_particle_velocity_spherical_theta,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_spherical_velocity_theta(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_velocity_spherical_theta']

    registry.add_field((ptype, "particle_spherical_velocity_theta"),
                       sampling_type="particle",
                       function=_particle_spherical_velocity_theta,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_velocity_spherical_phi(field, data):
        """The spherical phi component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        vel = data['relative_particle_velocity'].T
        phi = get_sph_phi(pos, normal)
        sphp = get_sph_phi_component(vel, phi, normal)
        return sphp

    registry.add_field((ptype, "particle_velocity_spherical_phi"),
                       sampling_type="particle",
                       function=_particle_velocity_spherical_phi,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_spherical_velocity_phi(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_spherical_velocity_theta']

    registry.add_field((ptype, "particle_spherical_velocity_phi"),
                       sampling_type="particle",
                       function=_particle_spherical_velocity_phi,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_position_cylindrical_radius(field, data):
        """The cylindrical radius component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        return data.ds.arr(get_cyl_r(pos, normal), 'code_length')

    registry.add_field((ptype, "particle_position_cylindrical_radius"),
                       sampling_type="particle",
                       function=_particle_position_cylindrical_radius,
                       units=unit_system["length"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_position_cylindrical_theta(field,data):
        """The cylindrical theta component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        return data.ds.arr(get_cyl_theta(pos, normal), "")

    registry.add_field((ptype, "particle_position_cylindrical_theta"),
                       sampling_type="particle",
                       function=_particle_position_cylindrical_theta,
                       units="",
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")])

    def _particle_position_cylindrical_z(field,data):
        """The cylindrical z component of the particle positions

        Relative to the coordinate system defined by the *normal* vector
        and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        return data.ds.arr(get_cyl_z(pos, normal), 'code_length')

    registry.add_field((ptype, "particle_position_cylindrical_z"),
                       sampling_type="particle",
                       function=_particle_position_cylindrical_z,
                       units=unit_system["length"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_velocity_cylindrical_radius(field, data):
        """The cylindrical radius component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        vel = data['relative_particle_velocity'].T
        theta = get_cyl_theta(pos, normal)
        cylr = get_cyl_r_component(vel, theta, normal)
        return cylr

    registry.add_field((ptype, "particle_velocity_cylindrical_radius"),
                       sampling_type="particle",
                       function=_particle_velocity_cylindrical_radius,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_velocity_cylindrical_theta(field, data):
        """The cylindrical theta component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        pos = data['relative_particle_position'].T
        vel = data['relative_particle_velocity'].T
        theta = get_cyl_theta(pos, normal)
        cylt = get_cyl_theta_component(vel, theta, normal)
        return cylt

    registry.add_field((ptype, "particle_velocity_cylindrical_theta"),
                       sampling_type="particle",
                       function=_particle_velocity_cylindrical_theta,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_cylindrical_velocity_theta(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, 'particle_velocity_cylindrical_theta']

    registry.add_field((ptype, "particle_cylindrical_velocity_theta"),
                       sampling_type="particle",
                       function=_particle_cylindrical_velocity_theta,
                       units="cm/s",
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_velocity_cylindrical_z(field, data):
        """The cylindrical z component of the particle velocities

        Relative to the coordinate system defined by the *normal* vector,
        *bulk_velocity* vector and *center* field parameters.
        """
        normal = data.get_field_parameter('normal')
        vel = data['relative_particle_velocity'].T
        cylz = get_cyl_z_component(vel, normal)
        return cylz

    registry.add_field((ptype, "particle_velocity_cylindrical_z"),
                       sampling_type="particle",
                       function=_particle_velocity_cylindrical_z,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    def _particle_cylindrical_velocity_z(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data[ptype, "particle_velocity_cylindrical_z"]

    registry.add_field((ptype, "particle_cylindrical_velocity_z"),
                       sampling_type="particle",
                       function=_particle_cylindrical_velocity_z,
                       units=unit_system["velocity"],
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])


def add_particle_average(registry, ptype, field_name,
                         weight = "particle_mass",
                         density = True):
    field_units = registry[ptype, field_name].units
    def _pfunc_avg(field, data):
        pos = data[ptype, "particle_position"]
        f = data[ptype, field_name]
        wf = data[ptype, weight]
        f *= wf
        v = data.deposit(pos, [f], method = "sum")
        w = data.deposit(pos, [wf], method = "sum")
        v /= w
        if density: v /= data["index", "cell_volume"]
        v[np.isnan(v)] = 0.0
        return v
    fn = ("deposit", "%s_avg_%s" % (ptype, field_name))
    registry.add_field(fn,
                       sampling_type="cell",
                       function=_pfunc_avg,
                       validators = [ValidateSpatial(0)],
                       units = field_units)
    return fn

def add_volume_weighted_smoothed_field(ptype, coord_name, mass_name,
        smoothing_length_name, density_name, smoothed_field, registry,
        nneighbors = 64, kernel_name = 'cubic'):
    unit_system = registry.ds.unit_system
    if kernel_name == 'cubic':
        field_name = ("deposit", "%s_smoothed_%s" % (ptype, smoothed_field))
    else:
        field_name = ("deposit", "%s_%s_smoothed_%s" % (ptype, kernel_name,
                      smoothed_field))
    field_units = registry[ptype, smoothed_field].units
    def _vol_weight(field, data):
        pos = data[ptype, coord_name]
        pos = pos.convert_to_units("code_length")
        mass = data[ptype, mass_name].in_base(unit_system.name)
        dens = data[ptype, density_name].in_base(unit_system.name)
        quan = data[ptype, smoothed_field]
        if hasattr(quan, "units"):
            quan = quan.convert_to_units(field_units)

        if smoothing_length_name is None:
            hsml = np.zeros(quan.shape, dtype='float64') - 1
            hsml = data.apply_units(hsml, "code_length")
        else:
            hsml = data[ptype, smoothing_length_name]
            hsml.convert_to_units("code_length")
        # This is for applying cutoffs, similar to in the SPLASH paper.
        smooth_cutoff = data["index","cell_volume"]**(1./3)
        smooth_cutoff.convert_to_units("code_length")
        # volume_weighted smooth operations return lists of length 1.
        rv = data.smooth(pos, [mass, hsml, dens, quan],
                         index_fields=[smooth_cutoff],
                         method="volume_weighted",
                         create_octree=True,
                         nneighbors=nneighbors,
                         kernel_name=kernel_name)[0]
        rv[np.isnan(rv)] = 0.0
        # Now some quick unit conversions.
        # This should be used when seeking a non-normalized value:
        rv /= hsml.uq**3 / hsml.uq.in_base(unit_system.name).uq**3
        rv = data.apply_units(rv, field_units)
        return rv
    registry.add_field(field_name,
                       sampling_type="cell",
                       function = _vol_weight,
                       validators = [ValidateSpatial(0)],
                       units = field_units)
    registry.find_dependencies((field_name,))
    return [field_name]

def add_nearest_neighbor_field(ptype, coord_name, registry, nneighbors = 64):
    field_name = (ptype, "nearest_neighbor_distance_%s" % (nneighbors))
    def _nth_neighbor(field, data):
        pos = data[ptype, coord_name]
        pos.convert_to_units("code_length")
        distances = 0.0 * pos[:,0]
        data.particle_operation(pos, [distances],
                         method="nth_neighbor",
                         nneighbors = nneighbors)
        # Now some quick unit conversions.
        return distances
    registry.add_field(field_name,
                       sampling_type="particle",
                       function = _nth_neighbor,
                       validators = [ValidateSpatial(0)],
                       units = "code_length")
    return [field_name]

def add_union_field(registry, ptype, field_name, units):
    """
    Create a field that is the concatenation of multiple particle types.
    This allows us to create fields for particle unions using alias names.
    """

    def _cat_field(field, data):
        return uconcatenate([data[dep_type, field_name]
                             for dep_type in data.ds.particle_types_raw])

    registry.add_field((ptype, field_name),
                       sampling_type="particle",
                       function=_cat_field,
                       units=units)
