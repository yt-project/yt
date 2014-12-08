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

from yt.units.yt_array import YTArray

from .derived_field import \
    ValidateParameter, \
    ValidateSpatial

from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component, \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, get_sph_r, \
    get_sph_theta, get_sph_phi, \
    periodic_dist, euclidean_dist, \
    resize_vector

from yt.funcs import just_one

from yt.utilities.lib.misc_utilities import obtain_rvec, obtain_rv_vec

def create_magnitude_field(registry, basename, field_units,
                           ftype = "gas", slice_info = None,
                           validators = None, particle_type=False):

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    def _magnitude(field, data):
        fn = field_components[0]
        mag = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            mag += data[fn] * data[fn]
        return np.sqrt(mag)

    registry.add_field((ftype, "%s_magnitude" % basename),
                       function = _magnitude, units = field_units,
                       validators = validators, particle_type = particle_type)

def create_squared_field(registry, basename, field_units,
                         ftype = "gas", slice_info = None,
                         validators = None, particle_type=False):

    field_components = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    def _squared(field, data):
        fn = field_components[0]
        squared  = data[fn] * data[fn]
        for idim in range(1, registry.ds.dimensionality):
            fn = field_components[idim]
            squared += data[fn] * data[fn]
        return squared

    registry.add_field((ftype, "%s_squared" % basename),
                       function = _squared, units = field_units,
                       validators = validators, particle_type = particle_type)

def create_vector_fields(registry, basename, field_units,
                         ftype = "gas", slice_info = None):
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
    sl_center = slice(1, -1, None)

    xn, yn, zn = [(ftype, "%s_%s" % (basename, ax)) for ax in 'xyz']

    # Is this safe?
    if registry.ds.dimensionality < 3:
        zn = ("index", "zeros")
    if registry.ds.dimensionality < 2:
        yn = ("index", "zeros")

    create_magnitude_field(registry, basename, field_units,
                           ftype=ftype, slice_info=slice_info)
        
    def _radial(field, data):
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data['index', 'spherical_theta']
        phi   = data['index', 'spherical_phi']
        rv = get_sph_r_component(vectors, theta, phi, normal)
        # Now, anywhere that radius is in fact zero, we want to zero out our
        # return values.
        rv[np.isnan(theta)] = 0.0
        return rv
    def _radial_absolute(field, data):
        return np.abs(data[ftype, "radial_%s" % basename])

    def _tangential(field, data):
        return np.sqrt(data[ftype, "%s_magnitude" % basename]**2.0
                     - data[ftype, "radial_%s" % basename]**2.0)

    registry.add_field((ftype, "radial_%s" % basename),
                       function = _radial, units = field_units,
                       validators=[ValidateParameter("normal"),
                                   ValidateParameter("center")])
    registry.add_field((ftype, "radial_%s_absolute" % basename),
                       function = _radial_absolute, units = field_units)
    registry.add_field((ftype, "tangential_%s" % basename),
                       function=_tangential, units = field_units)

    def _cp_vectors(ax):
        def _cp_val(field, data):
            vec = data.get_field_parameter("cp_%s_vec" % (ax))
            bv = data.get_field_parameter("bulk_%s" % basename)
            if bv == None: bv = np.zeros(3)
            tr  = (data[xn] - bv[0]) * vec[0]
            tr += (data[yn] - bv[1]) * vec[1]
            tr += (data[zn] - bv[2]) * vec[2]
            return tr
        return _cp_val
    registry.add_field((ftype, "cutting_plane_%s_x" % basename),
                       function = _cp_vectors('x'), units = field_units)
    registry.add_field((ftype, "cutting_plane_%s_y" % basename),
                       function = _cp_vectors('y'), units = field_units)
    registry.add_field((ftype, "cutting_plane_%s_z" % basename),
                       function = _cp_vectors('z'), units = field_units)

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

    registry.add_field((ftype, "%s_divergence" % basename),
                       function = _divergence, units = "1/s",
                       validators = [ValidateSpatial(1)])
    registry.add_field((ftype, "%s_divergence_absolute" % basename),
                       function = _divergence_abs, units = "1/s")

    def _tangential_over_magnitude(field, data):
        tr = data[ftype, "tangential_%s" % basename] / \
             data[ftype, "%s_magnitude" % basename]
        return np.abs(tr)
    registry.add_field((ftype, "tangential_over_%s_magnitude" % basename),
             function=_tangential_over_magnitude,
             take_log=False)

    def _cylindrical_radial(field, data):
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = resize_vector(data["index", 'cylindrical_theta'], vectors)
        return get_cyl_r_component(vectors, theta, normal)

    def _cylindrical_radial_absolute(field, data):
        return np.abs(_cylindrical_radial(field, data))

    registry.add_field((ftype, "cylindrical_radial_%s" % basename),
             function=_cylindrical_radial,
             units=field_units,
             validators=[ValidateParameter("normal")])
    registry.add_field((ftype, "cylindrical_radial_%s_absolute" % basename),
             function=_cylindrical_radial_absolute,
             units=field_units,
             validators=[ValidateParameter("normal")])

    def _cylindrical_tangential(field, data):
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data["index", 'cylindrical_theta'].copy()
        theta = np.tile(theta, (3,) + (1,)*len(theta.shape))
        return get_cyl_theta_component(vectors, theta, normal)

    def _cylindrical_tangential_absolute(field, data):
        return np.abs(_cylindrical_tangential(field, data))

    registry.add_field(
             (ftype, "cylindrical_tangential_%s" % basename),
             function=_cylindrical_tangential,
             units=field_units,
             validators=[ValidateParameter("normal")])
    registry.add_field(
             (ftype, "cylindrical_tangential_%s_absolute" % basename),
             function=_cylindrical_tangential_absolute,
             units=field_units,
             validators=[ValidateParameter("normal")])

def create_averaged_field(registry, basename, field_units,
                          ftype = "gas", slice_info = None,
                          validators = None,
                          weight = "cell_mass"):

    if validators is None:
        validators = []
    validators += [ValidateSpatial(1, [(ftype, basename)])]

    def _averaged_field(field, data):
        nx, ny, nz = data[(ftype, basename)].shape
        new_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2), dtype=np.float64),
                                (just_one(data[(ftype, basename)]) *
                                 just_one(data[(ftype, weight)])).units)
        weight_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2), dtype=np.float64),
                                   data[(ftype, weight)].units)
        i_i, j_i, k_i = np.mgrid[0:3, 0:3, 0:3]

        for i, j, k in zip(i_i.ravel(), j_i.ravel(), k_i.ravel()):
            sl = [slice(i, nx-(2-i)), slice(j, ny-(2-j)), slice(k, nz-(2-k))]
            new_field += data[(ftype, basename)][sl] * \
              data[(ftype, weight)][sl]
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
