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
    get_cyl_theta_component, \
    resize_vector

from yt.funcs import just_one

from yt.utilities.lib.misc_utilities import obtain_rv_vec

def create_magnitude_field(registry, basename, field_units,
                           ftype="gas", slice_info=None,
                           validators=None, particle_type=False):

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    def _magnitude(field, data):
        fn = field_components[0]
        mag = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            mag += data[fn] * data[fn]
        return np.sqrt(mag)

    registry.add_field((ftype, "%s_magnitude" % basename),
                       function=_magnitude, units=field_units,
                       validators=validators, particle_type=particle_type)

def create_squared_field(registry, basename, field_units,
                         ftype="gas", slice_info=None,
                         validators=None, particle_type=False):

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    def _squared(field, data):
        fn = field_components[0]
        squared = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            squared += data[fn] * data[fn]
        return squared

    registry.add_field((ftype, "%s_squared" % basename),
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

    create_magnitude_field(registry, basename, field_units,
                           ftype=ftype, slice_info=slice_info)

    def _spherical_radius_component(field, data):
        """The spherical radius component of the vector field

        Relative to the coordinate system defined by the *normal* vector,
        *center*, and *bulk_* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data['index', 'spherical_theta']
        phi = data['index', 'spherical_phi']
        rv = get_sph_r_component(vectors, theta, phi, normal)
        # Now, anywhere that radius is in fact zero, we want to zero out our
        # return values.
        rv[np.isnan(theta)] = 0.0
        return rv

    registry.add_field((ftype, "%s_spherical_radius" % basename),
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
        return np.sqrt(data[ftype, "%s_magnitude" % basename]**2.0 -
                       data[ftype, "%s_spherical_radius" % basename]**2.0)

    registry.add_field((ftype, "radial_%s" % basename),
                       function=_radial,
                       units=field_units,
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])

    registry.add_field((ftype, "radial_%s_absolute" % basename),
                       function=_radial_absolute,
                       units=field_units)

    registry.add_field((ftype, "tangential_%s" % basename),
                       function=_tangential,
                       units=field_units)

    def _spherical_theta_component(field, data):
        """The spherical theta component of the vector field

        Relative to the coordinate system defined by the *normal* vector,
        *center*, and *bulk_* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = resize_vector(data["index", "spherical_theta"], vectors)
        phi = resize_vector(data["index", "spherical_phi"], vectors)
        return get_sph_theta_component(vectors, theta, phi, normal)

    registry.add_field((ftype, "%s_spherical_theta" % basename),
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
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        phi = resize_vector(data["index", "spherical_phi"], vectors)
        return get_sph_phi_component(vectors, phi, normal)

    registry.add_field((ftype, "%s_spherical_phi" % basename),
                       function=_spherical_phi_component,
                       units=field_units,
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center"),
                                   ValidateParameter("bulk_%s" % basename)])

    def _cp_vectors(ax):
        def _cp_val(field, data):
            vec = data.get_field_parameter("cp_%s_vec" % (ax))
            bv = data.get_field_parameter("bulk_%s" % basename, None)
            if bv is None:
                bv = data.ds.arr(np.zeros(3), data[xn].units)
            tr  = (data[xn] - bv[0]) * vec.d[0]
            tr += (data[yn] - bv[1]) * vec.d[1]
            tr += (data[zn] - bv[2]) * vec.d[2]
            return tr
        return _cp_val

    registry.add_field((ftype, "cutting_plane_%s_x" % basename),
                       function=_cp_vectors('x'),
                       units=field_units)
    registry.add_field((ftype, "cutting_plane_%s_y" % basename),
                       function=_cp_vectors('y'),
                       units=field_units)
    registry.add_field((ftype, "cutting_plane_%s_z" % basename),
                       function=_cp_vectors('z'),
                       units=field_units)

    def _divergence(field, data):
        ds = div_fac * just_one(data["index", "dx"])
        f  = data[xn][sl_right,1:-1,1:-1]/ds
        f -= data[xn][sl_left ,1:-1,1:-1]/ds
        ds = div_fac * just_one(data["index", "dy"])
        f += data[yn][1:-1,sl_right,1:-1]/ds
        f -= data[yn][1:-1,sl_left ,1:-1]/ds
        ds = div_fac * just_one(data["index", "dz"])
        f += data[zn][1:-1,1:-1,sl_right]/ds
        f -= data[zn][1:-1,1:-1,sl_left ]/ds
        new_field = data.ds.arr(np.zeros(data[xn].shape, dtype=np.float64),
                                f.units)        
        new_field[1:-1,1:-1,1:-1] = f
        return new_field

    def _divergence_abs(field, data):
        return np.abs(data[ftype, "%s_divergence" % basename])

    field_units = Unit(field_units, registry=registry.ds.unit_registry)
    div_units = field_units / registry.ds.unit_system["length"]

    registry.add_field((ftype, "%s_divergence" % basename),
                       function=_divergence,
                       units=div_units,
                       validators=[ValidateSpatial(1)])
    
    registry.add_field((ftype, "%s_divergence_absolute" % basename),
                       function=_divergence_abs,
                       units=div_units)

    def _tangential_over_magnitude(field, data):
        tr = data[ftype, "tangential_%s" % basename] / \
             data[ftype, "%s_magnitude" % basename]
        return np.abs(tr)
    registry.add_field((ftype, "tangential_over_%s_magnitude" % basename),
             function=_tangential_over_magnitude,
             take_log=False)

    def _cylindrical_radius_component(field, data):
        """The cylindrical radius component of the vector field

        Relative to the coordinate system defined by the *normal* vector,
        *center*, and *bulk_* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = resize_vector(data["index", 'cylindrical_theta'], vectors)
        return get_cyl_r_component(vectors, theta, normal)

    registry.add_field((ftype, "%s_cylindrical_radius" % basename),
                       function=_cylindrical_radius_component,
                       units=field_units,
                       validators=[ValidateParameter("normal")])

    def _cylindrical_radial(field, data):
        """This field is deprecated and will be removed in a future version"""
        return data[ftype, '%s_cylindrical_radius' % basename]

    registry.add_field((ftype, "cylindrical_radial_%s" % basename),
                       function=_cylindrical_radial,
                       units=field_units)

    def _cylindrical_radial_absolute(field, data):
        """This field is deprecated and will be removed in a future version"""
        return np.abs(data[ftype, '%s_cylindrical_radius' % basename])

    registry.add_field((ftype, "cylindrical_radial_%s_absolute" % basename),
                       function=_cylindrical_radial_absolute,
                       units=field_units,
                       validators=[ValidateParameter("normal")])

    def _cylindrical_theta_component(field, data):
        """The cylindrical theta component of the vector field

        Relative to the coordinate system defined by the *normal* vector,
        *center*, and *bulk_* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data["index", 'cylindrical_theta'].copy()
        theta = np.tile(theta, (3,) + (1,)*len(theta.shape))
        return get_cyl_theta_component(vectors, theta, normal)

    registry.add_field((ftype, "%s_cylindrical_theta" % basename),
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
                       function=_cylindrical_tangential,
                       units=field_units)

    registry.add_field((ftype, "cylindrical_tangential_%s_absolute" % basename),
                       function=_cylindrical_tangential_absolute,
                       units=field_units)

    def _cylindrical_z_component(field, data):
        """The cylindrical z component of the vector field

        Relative to the coordinate system defined by the *normal* vector,
        *center*, and *bulk_* field parameters.
        """
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn), "bulk_%s" % basename)
        return get_cyl_z_component(vectors, normal)

    registry.add_field((ftype, "%s_cylindrical_z" % basename),
                       function=_cylindrical_z_component,
                       units=field_units,
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center"),
                                   ValidateParameter("bulk_%s" % basename)])


def create_averaged_field(registry, basename, field_units, ftype="gas",
                          slice_info=None, validators=None, weight="cell_mass"):

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
            sl = [slice(i, nx-(2-i)), slice(j, ny-(2-j)), slice(k, nz-(2-k))]
            new_field += data[(ftype, basename)][sl] * data[(ftype, weight)][sl]
            weight_field += data[(ftype, weight)][sl]

        # Now some fancy footwork
        new_field2 = data.ds.arr(np.zeros((nx, ny, nz)), 
                                 data[(ftype, basename)].units)
        new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
        return new_field2

    registry.add_field((ftype, "averaged_%s" % basename),
                       function=_averaged_field,
                       units=field_units,
                       validators=validators)
