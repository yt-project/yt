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

from yt.data_objects.yt_array import YTArray

def create_vector_fields(pf, registry, basename, field_units,
                         ftype = "gas", slice_info = None)

    xn, yn, zn = ["%s_%s" % (basename, ax) for ax in 'xyz']

    # Is this safe?
    if pf.dimensionality < 3:
        zn = "zeros"
    if pf.dimensionality < 2:
        yn = "zeros"

    def _magnitude(field, data):
        mag  = xn*xn
        mag += yn*yn
        mag += zn*zn
        return np.sqrt(mag)

    registry.add_field((ftype, "%s_magnitude" % basename),
                       function = _magnitude, units = units)

    def _radial(field, data):
        normal = data.get_field_parameter("normal")
        vectors = obtain_rv_vec(data, (xn, yn, zn),
                                "bulk_%s" % basename)
        theta = data['spherical_theta']
        phi   = data['spherical_phi']
        return get_sph_r_component(vectors, theta, phi, normal)
    def _radial_absolute(field, data):
        return np.abs(data["radial_%s" % basename)

    def _tangential(field, data):
        return np.sqrt(data["%s_magnitude" % basename]**2.0
                     - data["radial_%s" % basename]**2.0)

    registry.add_field((ftype, "radial_%s" % basename),
                       function = _radial, units = units)
    registry.add_field((ftype, "radial_%s_absolute" % basename),
                       function = _radial, units = units)
    registry.add_field((ftype, "tangential_%s" % basename),
                       function=_tangential, units=units)

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
                       function = _cp_vectors('x'), units=units)
    registry.add_field((ftype, "cutting_plane_%s_y" % basename),
                       function = _cp_vectors('y'), units=units)
    registry.add_field((ftype, "cutting_plane_%s_z" % basename),
                       function = _cp_vectors('z'), units=units)

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

def create_gradient_fields(pf, registry, basename, field_units,
                           ftype = "gas", slice_info = None)

