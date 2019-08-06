"""
This is where we define a handful of vector operations for fields.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from .derived_field import \
    ValidateParameter, \
    ValidateSpatial

from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component

from yt.funcs import \
    just_one, \
    handle_mks_cgs

from yt.geometry.geometry_handler import \
    is_curvilinear

from yt.utilities.lib.misc_utilities import obtain_relative_velocity_vector

def get_bulk(data, basename, unit):
    if data.has_field_parameter("bulk_%s" % basename):
        bulk = data.get_field_parameter("bulk_%s" % basename)
    else:
        bulk = [0, 0, 0]*unit
    return bulk


def create_magnitude_field(registry, basename, field_units,
                           ftype="gas", slice_info=None,
                           validators=None, particle_type=False):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in axis_order]

    def _magnitude(field, data):
        fn = field_components[0]
        if data.has_field_parameter('bulk_%s' % basename):
            fn = (fn[0], 'relative_%s' % fn[1])
        d = data[fn]
        mag = (d)**2
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            if data.has_field_parameter('bulk_%s' % basename):
                fn = (fn[0], 'relative_%s' % fn[1])
            mag += (data[fn])**2
        return np.sqrt(mag)

    if particle_type is True:
        sampling_type = 'particle'
    else:
        sampling_type = 'cell'

    registry.add_field((ftype, "%s_magnitude" % basename),
                       sampling_type=sampling_type,
                       function=_magnitude,
                       units=field_units,
                       validators=validators)

def create_relative_field(registry, basename, field_units, ftype='gas',
                          slice_info=None, validators=None,
                          particle_type=False):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in axis_order]

    def relative_vector(ax):
        def _relative_vector(field, data):
            iax = axis_order.index(ax)
            d = handle_mks_cgs(data[field_components[iax]], field.units)
            bulk = handle_mks_cgs(get_bulk(data, basename, d.unit_quantity),
                                  field.units)
            return d - bulk[iax]
        return _relative_vector

    if particle_type is True:
        sampling_type = 'particle'
    else:
        sampling_type = 'cell'

    for d in axis_order:
        registry.add_field((ftype, "relative_%s_%s" % (basename, d)),
                           sampling_type=sampling_type,
                           function=relative_vector(d),
                           units=field_units,
                           validators=validators)


def create_squared_field(registry, basename, field_units,
                         ftype="gas", slice_info=None,
                         validators=None, particle_type=False):

    axis_order = registry.ds.coordinates.axis_order

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in axis_order]

    def _squared(field, data):
        fn = field_components[0]
        if data.has_field_parameter('bulk_%s' % basename):
            fn = (fn[0], 'relative_%s' % fn[1])
        squared = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            squared += data[fn] * data[fn]
        return squared

    registry.add_field((ftype, "%s_squared" % basename), sampling_type="cell",
                       function=_squared, units=field_units,
                       validators=validators, particle_type=particle_type)

def create_vector_fields(registry, basename, field_units,
                         ftype="gas", slice_info=None):
    from yt.units.unit_object import Unit
    # slice_info would be the left, the right, and the factor.
    # For example, with the old Enzo-ZEUS fields, this would be:
    # slice(None, -2, None)
    # slice(1, -1, None)
    # 1.0
    # Otherwise, we default to a centered difference.
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info

    xn, yn, zn = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    # Is this safe?
    if registry.ds.dimensionality < 3:
        zn = ("index", "zeros")
    if registry.ds.dimensionality < 2:
        yn = ("index", "zeros")

    create_relative_field(
        registry, basename, field_units, ftype=ftype, slice_info=slice_info,
        validators=[ValidateParameter('bulk_%s' % basename)])

    create_magnitude_field(
        registry, basename, field_units, ftype=ftype, slice_info=slice_info,
        validators=[ValidateParameter('bulk_%s' % basename)])

    if not is_curvilinear(registry.ds.geometry):

        # The following fields are invalid for curvilinear geometries
        def _spherical_radius_component(field, data):
            """The spherical radius component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            theta = data['index', 'spherical_theta']
            phi = data['index', 'spherical_phi']
            rv = get_sph_r_component(vectors, theta, phi, normal)
            # Now, anywhere that radius is in fact zero, we want to zero out our
            # return values.
            rv[np.isnan(theta)] = 0.0
            return rv

        registry.add_field((ftype, "%s_spherical_radius" % basename),
                           sampling_type="cell",
                           function=_spherical_radius_component,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center"),
                                       ValidateParameter("bulk_%s" % basename)])

        def _radial(field, data):
            return data[ftype, "%s_spherical_radius" % basename]

        def _radial_absolute(field, data):
            return np.abs(data[ftype, "%s_spherical_radius" % basename])

        def _tangential(field, data):
            return np.sqrt(data[ftype, "%s_spherical_theta" % basename]**2.0 +
                           data[ftype, "%s_spherical_phi" % basename]**2.0)

        registry.add_field((ftype, "radial_%s" % basename),
                           sampling_type="cell",
                           function=_radial,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center")])

        registry.add_field((ftype, "radial_%s_absolute" % basename),
                           sampling_type="cell",
                           function=_radial_absolute,
                           units=field_units)

        registry.add_field((ftype, "tangential_%s" % basename),
                           sampling_type="cell",
                           function=_tangential,
                           units=field_units)

        def _spherical_theta_component(field, data):
            """The spherical theta component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            theta = data["index", "spherical_theta"]
            phi = data["index", "spherical_phi"]
            return get_sph_theta_component(vectors, theta, phi, normal)

        registry.add_field((ftype, "%s_spherical_theta" % basename),
                           sampling_type="cell",
                           function=_spherical_theta_component,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center"),
                                       ValidateParameter("bulk_%s" % basename)])

        def _spherical_phi_component(field, data):
            """The spherical phi component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            phi = data["index", "spherical_phi"]
            return get_sph_phi_component(vectors, phi, normal)

        registry.add_field((ftype, "%s_spherical_phi" % basename),
                           sampling_type="cell",
                           function=_spherical_phi_component,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center"),
                                       ValidateParameter("bulk_%s" % basename)])

        def _cp_vectors(ax):
            def _cp_val(field, data):
                vec = data.get_field_parameter("cp_%s_vec" % (ax))
                tr  = data[xn[0], 'relative_%s' % xn[1]] * vec.d[0]
                tr += data[yn[0], 'relative_%s' % yn[1]] * vec.d[1]
                tr += data[zn[0], 'relative_%s' % zn[1]] * vec.d[2]
                return tr
            return _cp_val

        for ax in 'xyz':
            registry.add_field((ftype, "cutting_plane_%s_%s" % (basename, ax)),
                               sampling_type="cell",
                               function=_cp_vectors(ax),
                               units=field_units)

        def _divergence(field, data):
            ds = div_fac * just_one(data["index", "dx"])
            f  = data[xn[0], 'relative_%s' % xn[1]][sl_right,1:-1,1:-1]/ds
            f -= data[xn[0], 'relative_%s' % xn[1]][sl_left ,1:-1,1:-1]/ds
            ds = div_fac * just_one(data["index", "dy"])
            f += data[yn[0], 'relative_%s' % yn[1]][1:-1,sl_right,1:-1]/ds
            f -= data[yn[0], 'relative_%s' % yn[1]][1:-1,sl_left ,1:-1]/ds
            ds = div_fac * just_one(data["index", "dz"])
            f += data[zn[0], 'relative_%s' % zn[1]][1:-1,1:-1,sl_right]/ds
            f -= data[zn[0], 'relative_%s' % zn[1]][1:-1,1:-1,sl_left ]/ds
            new_field = data.ds.arr(
                np.zeros(data[xn].shape, dtype=np.float64), f.units)
            new_field[1:-1,1:-1,1:-1] = f
            return new_field

        def _divergence_abs(field, data):
            return np.abs(data[ftype, "%s_divergence" % basename])

        field_units = Unit(field_units, registry=registry.ds.unit_registry)
        div_units = field_units / registry.ds.unit_system["length"]

        registry.add_field((ftype, "%s_divergence" % basename),
                           sampling_type="cell",
                           function=_divergence,
                           units=div_units,
                           validators=[ValidateSpatial(1),
                                       ValidateParameter('bulk_%s' % basename)])

        registry.add_field((ftype, "%s_divergence_absolute" % basename),
                           sampling_type="cell",
                           function=_divergence_abs,
                           units=div_units)

        def _tangential_over_magnitude(field, data):
            tr = (data[ftype, "tangential_%s" % basename] /
                  data[ftype, '%s_magnitude' % basename])
            return np.abs(tr)
        registry.add_field((ftype, "tangential_over_%s_magnitude" % basename),
                           sampling_type="cell",
                           function=_tangential_over_magnitude,
                           take_log=False)

        def _cylindrical_radius_component(field, data):
            """The cylindrical radius component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            theta = data["index", 'cylindrical_theta']
            return get_cyl_r_component(vectors, theta, normal)

        registry.add_field((ftype, "%s_cylindrical_radius" % basename),
                           sampling_type="cell",
                           function=_cylindrical_radius_component,
                           units=field_units,
                           validators=[ValidateParameter("normal")])

        def _cylindrical_radial(field, data):
            """This field is deprecated and will be removed in a future version"""
            return data[ftype, '%s_cylindrical_radius' % basename]

        registry.add_field((ftype, "cylindrical_radial_%s" % basename),
                           sampling_type="cell",
                           function=_cylindrical_radial,
                           units=field_units)

        def _cylindrical_radial_absolute(field, data):
            """This field is deprecated and will be removed in a future version"""
            return np.abs(data[ftype, '%s_cylindrical_radius' % basename])

        registry.add_field((ftype, "cylindrical_radial_%s_absolute" % basename),
                           sampling_type="cell",
                           function=_cylindrical_radial_absolute,
                           units=field_units,
                           validators=[ValidateParameter("normal")])

        def _cylindrical_theta_component(field, data):
            """The cylindrical theta component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            theta = data["index", 'cylindrical_theta'].copy()
            theta = np.tile(theta, (3,) + (1,)*len(theta.shape))
            return get_cyl_theta_component(vectors, theta, normal)

        registry.add_field((ftype, "%s_cylindrical_theta" % basename),
                           sampling_type="cell",
                           function=_cylindrical_theta_component,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center"),
                                       ValidateParameter("bulk_%s" % basename)])

        def _cylindrical_tangential(field, data):
            """This field is deprecated and will be removed in a future release"""
            return data["%s_cylindrical_theta" % basename]

        def _cylindrical_tangential_absolute(field, data):
            """This field is deprecated and will be removed in a future release"""
            return np.abs(data['cylindrical_tangential_%s' % basename])

        registry.add_field((ftype, "cylindrical_tangential_%s" % basename),
                           sampling_type="cell",
                           function=_cylindrical_tangential,
                           units=field_units)

        registry.add_field(
            (ftype, "cylindrical_tangential_%s_absolute" % basename),
            sampling_type="cell",
            function=_cylindrical_tangential_absolute,
            units=field_units)

        def _cylindrical_z_component(field, data):
            """The cylindrical z component of the vector field

            Relative to the coordinate system defined by the *normal* vector,
            *center*, and *bulk_* field parameters.
            """
            normal = data.get_field_parameter("normal")
            vectors = obtain_relative_velocity_vector(
                data, (xn, yn, zn), "bulk_%s" % basename)
            return get_cyl_z_component(vectors, normal)

        registry.add_field((ftype, "%s_cylindrical_z" % basename),
                           sampling_type="cell",
                           function=_cylindrical_z_component,
                           units=field_units,
                           validators=[ValidateParameter("normal"),
                                       ValidateParameter("center"),
                                       ValidateParameter("bulk_%s" % basename)])

    else: # Create Cartesian fields for curvilinear coordinates

        def _cartesian_x(field,data):

          if registry.ds.geometry is "polar":

            return data["%s_r" % basename] * np.cos(data["theta"])

          elif registry.ds.geometry == "cylindrical":

            if data.ds.dimensionality==2:
              return data["%s_r" % basename]
            elif data.ds.dimensionality==3:
              return (data["%s_r" % basename] * np.cos(data["theta"]) -
                      data["%s_theta" % basename] * np.sin(data["theta"]))

          elif registry.ds.geometry == "spherical":

            if data.ds.dimensionality==2:
              return (data["%s_r" % basename] * np.sin(data["theta"]) +
                      data["%s_theta" % basename] * np.cos(data["theta"]))
            elif data.ds.dimensionality==3:
              return (data["%s_r" % basename] * np.sin(data["theta"]) * np.cos(data["phi"]) +
                      data["%s_theta" % basename] * np.cos(data["theta"]) * np.cos(["phi"]) -
                      data["%s_phi" % basename] * np.sin(data["phi"]))

        # it's redundant to define a cartesian x field for 1D data
        if registry.ds.dimensionality > 1:
          registry.add_field((ftype, "%s_cartesian_x" % basename), sampling_type="cell", 
                             function=_cartesian_x,
                             units=field_units,
                             display_field=True)

        def _cartesian_y(field,data):

          if registry.ds.geometry=="polar":

            return data["%s_r" % basename] * np.sin(data["theta"])

          elif registry.ds.geometry == "cylindrical":

            if data.ds.dimensionality==2:
              return data["%s_z" % basename]
            elif data.ds.dimensionality==3:
              return (data["%s_r" % basename] * np.sin(data["theta"]) +
                      data["%s_theta" % basename] * np.cos(data["theta"]))

          elif registry.ds.geometry == "spherical":

            if data.ds.dimensionality==2:
              return (data["%s_r" % basename] * np.cos(data["theta"]) -
                     data["%s_theta" % basename] * np.sin(data["theta"]))
            elif data.ds.dimensionality==3:
              return (data["%s_r" % basename] * np.sin(data["theta"]) * np.sin(data["phi"]) +
                      data["%s_theta" % basename] * np.cos(data["theta"]) * np.sin(["phi"]) +
                      data["%s_phi" % basename] * np.cos(data["phi"]))

        if registry.ds.dimensionality >= 2:
          registry.add_field((ftype, "%s_cartesian_y" % basename), sampling_type="cell", 
                             function=_cartesian_y,
                             units=field_units,
                             display_field=True)

        def _cartesian_z(field,data):

          if registry.ds.geometry == "cylindrical":
            return data["%s_z" % basename]
          elif registry.ds.geometry == "spherical":
            return (data["%s_r" % basename] * np.cos(data["theta"]) -
                    data["%s_theta" % basename] * np.sin(data["theta"]))

        if registry.ds.dimensionality == 3:
          registry.add_field((ftype, "%s_cartesian_z" % basename), sampling_type="cell", 
                             function=_cartesian_z,
                             units=field_units,
                             display_field=True)


def create_averaged_field(registry, basename, field_units, ftype="gas",
                          slice_info=None, validators=None,
                          weight="cell_mass"):

    if validators is None:
        validators = []
    validators += [ValidateSpatial(1, [(ftype, basename)])]

    def _averaged_field(field, data):
        nx, ny, nz = data[(ftype, basename)].shape
        new_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2), dtype=np.float64),
                                (just_one(data[(ftype, basename)]) *
                                 just_one(data[(ftype, weight)])).units)
        weight_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2),
                                            dtype=np.float64),
                                   data[(ftype, weight)].units)
        i_i, j_i, k_i = np.mgrid[0:3, 0:3, 0:3]

        for i, j, k in zip(i_i.ravel(), j_i.ravel(), k_i.ravel()):
            sl = (slice(i, nx-(2-i)), slice(j, ny-(2-j)), slice(k, nz-(2-k)))
            new_field += data[(ftype, basename)][sl] * data[(ftype, weight)][sl]
            weight_field += data[(ftype, weight)][sl]

        # Now some fancy footwork
        new_field2 = data.ds.arr(np.zeros((nx, ny, nz)), 
                                 data[(ftype, basename)].units)
        new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
        return new_field2

    registry.add_field((ftype, "averaged_%s" % basename),
                       sampling_type="cell",
                       function=_averaged_field,
                       units=field_units,
                       validators=validators)
