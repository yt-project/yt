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

from yt.units.yt_array import YTArray
from yt.utilities.lib.misc_utilities import \
    obtain_rvec, obtain_rv_vec
from yt.utilities.math_utils import resize_vector
from yt.utilities.cosmology import Cosmology
from yt.fields.derived_field import \
    ValidateParameter

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
def setup_magnetic_field_fields(registry, ftype = "gas", slice_info = None):
    def _magnetic_energy(field,data):
        """This assumes that your front end has provided Bx, By, Bz in
        units of Gauss. If you use MKS, make sure to write your own
        magnetic_energy field to deal with non-unitary \mu_0.
        """
        return (data[ftype,"magnetic_field_x"]**2 +
                data[ftype,"magnetic_field_y"]**2 +
                data[ftype,"magnetic_field_z"]**2)/(8*np.pi)
    registry.add_field((ftype, "magnetic_energy"),
             function=_magnetic_energy,
             units="erg / cm**3")

    def _plasma_beta(field,data):
        """This assumes that your front end has provided Bx, By, Bz in
        units of Gauss. If you use MKS, make sure to write your own
        plasma_beta field to deal with non-unitary \mu_0.
        """
        return data[ftype,'pressure']/data[ftype,'magnetic_energy']
    registry.add_field((ftype, "plasma_beta"),
             function=_plasma_beta,
             units="")

    def _magnetic_pressure(field,data):
        return data[ftype,'magnetic_energy']
    registry.add_field((ftype, "magnetic_pressure"),
             function=_magnetic_pressure,
             units="erg / cm**3")

    def _magnetic_field_strength(field,data):
        """This assumes that your front end has provided Bx, By, Bz in
        units of Gauss. If you use MKS, make sure to write your own
        PlasmaBeta field to deal with non-unitary \mu_0.
        """
        return np.sqrt(8.*np.pi*data[ftype,"magnetic_energy"])
    registry.add_field((ftype,"magnetic_field_strength"),
                       function=_magnetic_field_strength,
                       units = "gauss")

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
             units="gauss",
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
             units="gauss",
             validators=[ValidateParameter("normal")])

    def _alfven_speed(field,data):
        """This assumes that your front end has provided Bx, By, Bz in
        units of Gauss. If you use MKS, make sure to write your own
        alfven_speed field to deal with non-unitary \mu_0.
        """
        return data[ftype,'magnetic_field_strength']/np.sqrt(4.*np.pi*data[ftype,'density'])
    registry.add_field((ftype, "alfven_speed"), function=_alfven_speed,
                       units="cm/s")

    def _mach_alfven(field,data):
        return data[ftype,'velocity_magnitude']/data[ftype,'alfven_speed']
    registry.add_field((ftype, "mach_alfven"), function=_mach_alfven,
                       units="dimensionless")

