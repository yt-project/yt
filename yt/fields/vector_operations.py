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

from yt.data_objects.yt_array import YTArray

from .derived_field import \
    ValidateParameter, \
    ValidateSpatial

from .field_plugin_registry import \
    register_field_plugin

def create_vector_fields(registry, basename, field_units,
                         ftype = "gas", slice_info = None):

    xn, yn, zn = ["%s_%s" % (basename, ax) for ax in 'xyz']

    # Is this safe?
    if registry.pf.dimensionality < 3:
        zn = "zeros"
    if registry.pf.dimensionality < 2:
        yn = "zeros"

    def _magnitude(field, data):
        mag  = data[xn] * data[xn]
        mag += data[yn] * data[yn]
        mag += data[zn] * data[zn]
        return np.sqrt(mag)

    registry.add_field((ftype, "%s_magnitude" % basename),
                       function = _magnitude, units = field_units)

    def _radial(field, data):
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data['spherical_theta']
        phi   = data['spherical_phi']
        return get_sph_r_component(vectors, theta, phi, normal)
    def _radial_absolute(field, data):
        return np.abs(data["radial_%s" % basename])

    def _tangential(field, data):
        return np.sqrt(data["%s_magnitude" % basename]**2.0
                     - data["radial_%s" % basename]**2.0)

    registry.add_field((ftype, "radial_%s" % basename),
                       function = _radial, units = field_units)
    registry.add_field((ftype, "radial_%s_absolute" % basename),
                       function = _radial, units = field_units)
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
    def _divergence(field, data):
        ds = div_fac * just_one(data["dx"])
        f  = data[xn][sl_right,1:-1,1:-1]/ds
        f -= data[xn][sl_left ,1:-1,1:-1]/ds
        ds = div_fac * just_one(data["dy"])
        f += data[yn][1:-1,sl_right,1:-1]/ds
        f -= data[yn][1:-1,sl_left ,1:-1]/ds
        ds = div_fac * just_one(data["dz"])
        f += data[zn][1:-1,1:-1,sl_right]/ds
        f -= data[zn][1:-1,1:-1,sl_left ]/ds
        new_field = YTArray(np.zeros(data[xn].shape, dtype=np.float64),
                            '1/s')
        new_field[1:-1,1:-1,1:-1] = f
        return new_field
    def _divergence_abs(field, data):
        return np.abs(data[ftype, "%s_divergence" % basename])

    registry.add_field((ftype, "%s_divergence" % basename),
                       function = _divergence, units = "1/s")
    registry.add_field((ftype, "%s_divergence_absolute" % basename),
                       function = _divergence, units = "1/s")

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
        theta = resize_vector(data['cylindrical_theta'], vectors)
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
        theta = data['cylindrical_theta'].copy()
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

