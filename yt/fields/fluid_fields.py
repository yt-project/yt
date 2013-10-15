"""
Here are some fields that are specific to fluids.



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
    periodic_dist, euclidean_dist

@register_field_plugin
def setup_fluid_fields(registry, ftype = "gas", slice_info = None):
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
        sl_left, sl_right, div_face = slice_info

    def _sound_speed(field, data):
        tr = data.pf.gamma * data[ftype, "pressure"] / data[ftype, "density"]
        return np.sqrt(tr)
    registry.add_field((ftype, "sound_speed"),
             function=_sound_speed,
             units="cm/s")

    def _radial_mach_number(field, data):
        """ M{|v|/t_sound} """
        tr = data[ftype, "radial_velocity"] / data[ftype, "sound_speed"]
        return np.abs(tr)
    registry.add_field((ftype, "radial_mach_number"),
             function=_radial_mach_number,
             units = "")

    def _mach_number(field, data):
        """ M{|v|/t_sound} """
        return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]
    registry.add_field((ftype, "mach_number"),
            function=_mach_number,
            units = "")

    def _courant_time_step(field, data):
        t1 = data["dx"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_x"]))
        t2 = data["dy"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_y"]))
        t3 = data["dz"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_z"]))
        tr = np.minimum(np.minimum(t1, t2), t3)
        return field.apply_units(tr)

    registry.add_field((ftype, "courant_time_step"),
             function=_courant_time_step,
             units="s")

    def _pressure(field, data):
        """ M{(Gamma-1.0)*rho*E} """
        tr = (data.pf.gamma - 1.0) \
           * (data[ftype, "density"] * data[ftype, "thermal_energy"])
        return tr

    registry.add_field((ftype, "pressure"),
             function=_pressure,
             units="dyne/cm**2")

    def _entropy(field, data):
        mw = data.get_field_parameter("mu")
        if mw is None:
            mw = 1.0
        mw = mh
        gammam1 = data.pf.gamma - 1.0
        tr = kboltz * data[ftype, "temperature"] / \
               ((data[ftype, "density"]/mw)**gammam1)
        # Gamma here needs some fancy units.
        # TODO: Add fancy units for Gamma!
        return field.apply_units(tr)
    registry.add_field((ftype, "entropy"),
             units="erg/K",
             function=_entropy)

    # This may not be appropriate to have an 'ftype' for.
    def _mean_molecular_weight(field,data):
        return (data[ftype, "density"] / (mh * data[ftype, "number_density"]))
    registry.add_field((ftype, "mean_molecular_weight"),
              function=_mean_molecular_weight,
              units=r"")

    def _vorticity_squared(field, data):
        # We need to set up stencils
        new_field = YTArray(np.zeros(data[ftype, "velocity_x"].shape), 'cm/s')
        dvzdy = (data[ftype, "velocity_z"][1:-1,sl_right,1:-1] -
                 data[ftype, "velocity_z"][1:-1,sl_left,1:-1]) \
                 / (div_fac*just_one(data["dy"]))
        dvydz = (data[ftype, "velocity_y"][1:-1,1:-1,sl_right] -
                 data[ftype, "velocity_y"][1:-1,1:-1,sl_left]) \
                 / (div_fac*just_one(data["dz"]))
        new_field[1:-1,1:-1,1:-1] += (dvzdy - dvydz)**2.0
        del dvzdy, dvydz
        dvxdz = (data[ftype, "velocity_x"][1:-1,1:-1,sl_right] -
                 data[ftype, "velocity_x"][1:-1,1:-1,sl_left]) \
                 / (div_fac*just_one(data["dz"]))
        dvzdx = (data[ftype, "velocity_z"][sl_right,1:-1,1:-1] -
                 data[ftype, "velocity_z"][sl_left,1:-1,1:-1]) \
                 / (div_fac*just_one(data["dx"]))
        new_field[1:-1,1:-1,1:-1] += (dvxdz - dvzdx)**2.0
        del dvxdz, dvzdx
        dvydx = (data[ftype, "velocity_y"][sl_right,1:-1,1:-1] -
                 data[ftype, "velocity_y"][sl_left,1:-1,1:-1]) \
                 / (div_fac*just_one(data["dx"]))
        dvxdy = (data[ftype, "velocity_x"][1:-1,sl_right,1:-1] -
                 data[ftype, "velocity_x"][1:-1,sl_left,1:-1]) \
                 / (div_fac*just_one(data["dy"]))
        new_field[1:-1,1:-1,1:-1] += (dvydx - dvxdy)**2.0
        del dvydx, dvxdy
        new_field = np.abs(new_field)
        return new_field
    registry.add_field((ftype, "vorticity_squared"),
             function=_vorticity_squared,
             validators=[ValidateSpatial(1,
                         ["velocity_x","velocity_y","velocity_z"])],
             units="s**-2")

    setup_gradient_fields(registry, (ftype, "pressure"), "dyne/cm**2",
                          slice_info)

    setup_gradient_fields(registry, (ftype, "density"), "g / cm**3",
                          slice_info)

def setup_gradient_fields(registry, field, field_units, slice_info = None):
    assert(isinstance(field, tuple))
    ftype, fname = field
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_face = slice_info

    slice_3d = [slice(1, -1), slice(1, -1), slice(1, -1)]

    def grad_func(axi, ax):
        slice_3dl = slice_3d[:]
        slice_3dr = slice_3d[:]
        slice_3dl[axi] = sl_left
        slice_3dr[axi] = sl_right
        def func(field, data):
            # We need to set up stencils
            # This is based on enzo parameters and should probably be changed.    
            new_field = YTArray(np.zeros(data[field].shape, dtype=np.float64),
                                field_units)
            ds = div_fac * just_one(data['dx'])
            new_field[slice_3d]  = data[field][slice_3dr]/ds
            new_field[slice_3d] -= data[field][slice_3dl]/ds
            return new_field

    for axi, ax in enumerate('xyz'):
        f = grad_func(axi, ax)
        registry.add_field((ftype, "%s_gradient_%s" % (fname, ax)),
                 function = f,
                 validators = [ValidateSpatial(1, [field])],
                 units = "%s / cm" % field_units)
    
    def _gradient_magnitude(field, data):
        return np.sqrt(data[ftype, "%s_gradient_x" % fname]**2 +
                       data[ftype, "%s_gradient_y" % fname]**2 +
                       data[ftype, "%s_gradient_z" % fname]**2)
    registry.add_field((ftype, "%s_gradient_magnitude" % fname),
             function = _gradient_magnitude,
             validators = [ValidateSpatial(1, [field])],
             units = "%s / cm" % field_units)
