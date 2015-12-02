"""
Magnetic field ... er, fields.



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
    ValidateParameter

from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.math_utils import \
    get_sph_theta_component, \
    get_sph_phi_component

@register_field_plugin
def setup_magnetic_field_fields(registry, ftype = "gas", slice_info = None):
    from yt.utilities.physical_constants import mu_0
    unit_system = registry.ds.unit_system
    if str(unit_system) == "mks":
        mag_units = unit_system["magnetic_field_mks"]
        mag_fac = 1.0/mu_0
    else:
        mag_fac = 1.0/(4.0*np.pi)
        mag_units = unit_system["magnetic_field_cgs"]

    def _magnetic_energy(field,data):
        return 0.5*mag_fac*(data[ftype,"magnetic_field_x"]**2 +
                            data[ftype,"magnetic_field_y"]**2 +
                            data[ftype,"magnetic_field_z"]**2)
    registry.add_field((ftype, "magnetic_energy"),
             function=_magnetic_energy,
             units=unit_system["pressure"])

    def _plasma_beta(field,data):
        return data[ftype,'pressure']/data[ftype,'magnetic_energy']
    registry.add_field((ftype, "plasma_beta"),
             function=_plasma_beta,
             units="")

    def _magnetic_pressure(field,data):
        return data[ftype,'magnetic_energy']
    registry.add_field((ftype, "magnetic_pressure"),
             function=_magnetic_pressure,
             units=unit_system["pressure"])

    def _magnetic_field_strength(field,data):
        return np.sqrt(2.*mag_fac*data[ftype,"magnetic_energy"])
    registry.add_field((ftype,"magnetic_field_strength"),
                       function=_magnetic_field_strength,
                       units = mag_units)

    def _magnetic_field_poloidal(field,data):
        normal = data.get_field_parameter("normal")
        d = data[ftype,'magnetic_field_x']
        Bfields = data.ds.arr(
                    [data[ftype,'magnetic_field_x'],
                     data[ftype,'magnetic_field_y'],
                     data[ftype,'magnetic_field_z']],
                     d.units)
        
        theta = data["index", 'spherical_theta']
        phi   = data["index", 'spherical_phi']
        
        return get_sph_theta_component(Bfields, theta, phi, normal)

    registry.add_field((ftype, "magnetic_field_poloidal"),
             function=_magnetic_field_poloidal,
             units=mag_units,
             validators=[ValidateParameter("normal")])

    def _magnetic_field_toroidal(field,data):
        normal = data.get_field_parameter("normal")
        d = data[ftype,'magnetic_field_x']
        Bfields = data.ds.arr(
                    [data[ftype,'magnetic_field_x'],
                     data[ftype,'magnetic_field_y'],
                     data[ftype,'magnetic_field_z']],
                     d.units)
        
        phi = data["index", 'spherical_phi']
        return get_sph_phi_component(Bfields, phi, normal)

    registry.add_field((ftype, "magnetic_field_toroidal"),
             function=_magnetic_field_toroidal,
             units=mag_units,
             validators=[ValidateParameter("normal")])

    def _alfven_speed(field,data):
        return data[ftype,'magnetic_field_strength']/np.sqrt(mag_fac*data[ftype,'density'])
    registry.add_field((ftype, "alfven_speed"), function=_alfven_speed,
                       units=mag_units)

    def _mach_alfven(field,data):
        return data[ftype,'velocity_magnitude']/data[ftype,'alfven_speed']
    registry.add_field((ftype, "mach_alfven"), function=_mach_alfven,
                       units="dimensionless")

