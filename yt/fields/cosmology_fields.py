"""
Cosmology related fields.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.cosmology import \
     Cosmology
from yt.utilities.physical_constants import \
     speed_of_light_cgs
    
@register_field_plugin
def setup_cosmology_fields(registry, ftype = "gas", slice_info = None):
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

    def _matter_density(field, data):
        return data[ftype, "density"] + \
          data[ftype, "dark_matter_density"]

    registry.add_field((ftype, "matter_density"),
                       function=_matter_density, 
                       units="g/cm**3")
        
    # rho_total / rho_cr(z).
    def _overdensity(field, data):
        # consider moving cosmology object to pf attribute
        co = Cosmology(hubble_constant=data.pf.hubble_constant,
                       omega_matter=data.pf.omega_matter,
                       omega_lambda=data.pf.omega_lambda)
        return data["matter_density"] / co.critical_density(data.pf.current_redshift)
    
    registry.add_field((ftype, "overdensity"),
                       function=_overdensity,
                       units="")

    # rho_baryon / <rho_baryon>
    def _baryon_overdensity(field, data):
        omega_baryon = data.get_field_parameter("omega_baryon")
        # consider moving cosmology object to pf attribute
        co = Cosmology(hubble_constant=data.pf.hubble_constant,
                       omega_matter=data.pf.omega_matter,
                       omega_lambda=data.pf.omega_lambda,
                       unit_registry=data.pf.unit_registry)
        # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
        # mean density(z) ~ omega_matter * (1 + z)^3
        return data["density"] / omega_baryon_now / co.critical_density(0.0) / \
          (1.0 + data.pf.hubble_constant)**3

    registry.add_field("baryon_overdensity",
                       function=_baryon_overdensity,
                       units="")

    # rho_matter / <rho_matter>
    def _matter_overdensity(field, data):
        # consider moving cosmology object to pf attribute
        co = Cosmology(hubble_constant=data.pf.hubble_constant,
                       omega_matter=data.pf.omega_matter,
                       omega_lambda=data.pf.omega_lambda,
                       unit_registry=data.pf.unit_registry)
        # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
        # mean density(z) ~ omega_matter * (1 + z)^3
        return data["density"] / data.pf.omega_matter / co.critical_density(0.0) / \
          (1.0 + data.pf.hubble_constant)**3

    registry.add_field("matter_overdensity",
                       function=_matter_overdensity,
                       units="")
    
    # Weak lensing convergence.
    # Eqn 4 of Metzler, White, & Loken (2001, ApJ, 547, 560).
    # This needs to be checked for accuracy.
    def _weak_lensing_convergence(field, data):
        # consider moving cosmology object to pf attribute
        co = Cosmology(hubble_constant=data.pf.hubble_constant,
                       omega_matter=data.pf.omega_matter,
                       omega_lambda=data.pf.omega_lambda,
                       unit_registry=data.pf.unit_registry)
        observer_redshift = data.get_field_parameter('observer_redshift')
        source_redshift = data.get_field_parameter('source_redshift')
        
        # observer to lens
        dl = co.angular_diameter_distance(observer_redshift, data.pf.current_redshift)
        # observer to source
        ds = co.angular_diameter_distance(observer_redshift, source_redshift)
        # lens to source
        dls = co.angular_diameter_distance(data.pf.current_redshift, source_redshift)

        # removed the factor of 1 / a to account for the fact that we are projecting 
        # with a proper distance.
        return 1.5 * (co.hubble_constant / speed_of_light_cgs)**2 * (dl * dls / ds) * \
          data["matter_overdensity"]
       
    registry.add_field("weak_lensing_convergence",
                       function=_weak_lensing_convergence,
                       units="1/cm")
