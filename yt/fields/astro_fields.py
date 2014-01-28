"""
Astronomy and astrophysics fields.



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

from yt.utilities.physical_constants import \
    mh, \
    me, \
    sigma_thompson, \
    clight, \
    kboltz, \
    G, \
    speed_of_light_cgs
    
@register_field_plugin
def setup_astro_fields(registry, ftype = "gas", slice_info = None):
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
    
    def _dynamical_time(field, data):
        """
        sqrt(3 pi / (16 G rho))
        """
        return np.sqrt(3.0 * np.pi / (16.0 * G * data["density"]))

    registry.add_field("dynamical_time",
                       function=_dynamical_time,
                       units="s")

    # def _jeans_mass(field, data):
    #     MJ_constant = (((5.0 * kboltz) / (G * mh)) ** (1.5)) * \
    #       (3.0 / (4.0 * np.pi)) ** (0.5)
    #     u = (MJ_constant * \
    #          ((data["temperature"] / data["mean_molecular_weight"])**(1.5)) * \
    #          (data["density"]**(-0.5)))
    #     return u

    # registry.add_field("jeans_mass",
    #                    function=_jeans_mass,
    #                    units="g")

    def _chandra_emissivity(field, data):
        logT0 = np.log10(data["temperature"]) - 7
        return ( data["number_density"].astype(np.float64)**2
                 * ( 10**(-0.0103 * logT0**8
                          +0.0417 * logT0**7
                          -0.0636 * logT0**6
                          +0.1149 * logT0**5
                          -0.3151 * logT0**4
                          +0.6655 * logT0**3
                          -1.1256 * logT0**2
                          +1.0026 * logT0**1
                          -0.6984 * logT0)
                     + data["metallicity"] * 10**(0.0305 * logT0**11
                                                  -0.0045 * logT0**10
                                                  -0.3620 * logT0**9
                                                  +0.0513 * logT0**8
                                                  +1.6669 * logT0**7
                                                  -0.3854 * logT0**6
                                                  -3.3604 * logT0**5
                                                  +0.4728 * logT0**4
                                                  +4.5774 * logT0**3
                                                  -2.3661 * logT0**2
                                                  -1.6667 * logT0**1
                                                  -0.2193 * logT0) ) )

    def _convert_chandra_emissivity(data):
        return 1.0  # 1.0e-23*0.76**2

    #add_field("chandra_emissivity", function=_chandra_emissivity,
    #          convert_function=_convert_chandra_emissivity,
    #          projection_conversion="1")

    def _xray_emissivity(field, data):
        return ( data["density"].astype(np.float64)**2
                 * data["temperature"]**0.5 )

    def _convert_xray_emissivity(data):
        return 2.168e60  #TODO: cnvert me to constants

    #add_field("xray_emissivity", function=_xray_emissivity,
    #          convert_function=_convert_xray_emissivity,
    #          projection_conversion="1")

    def _sz_kinetic(field, data):
        vel_axis = data.get_field_parameter("axis")
        if vel_axis > 2:
            raise NeedsParameter(["axis"])
        vel = data["velocity_%s" % ({0: "x", 1: "y", 2: "z"}[vel_axis])]
        return (vel * data["density"])

    def _convert_sz_kinetic(data):
        return 0.88 * sigma_thompson / mh / clight

    #add_field("sz_kinetic", function=_sz_kinetic,
    #          convert_function=_convert_sz_kinetic,
    #          validators=[ValidateParameter("axis")])

    def _szy(field, data):
        return data["density"] * data["temperature"]

    def _convert_szy(data):
        conv = 0.88 / mh * kboltz / (me * clight*clight) * sigma_thompson
        return conv

    #add_field("szy", function=_szy, convert_function=_convert_szy)
