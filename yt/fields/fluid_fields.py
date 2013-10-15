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

from yt.data_objects.derived_fields import \
    ValidateParameter

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

def setup_fluid_fields(registry, ftype = "gas"):
    
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

