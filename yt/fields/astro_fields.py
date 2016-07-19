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

from .derived_field import \
    ValidateParameter
from .field_exceptions import \
    NeedsParameter
from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.physical_constants import \
    mh, \
    me, \
    sigma_thompson, \
    clight, \
    kboltz, \
    G

@register_field_plugin
def setup_astro_fields(registry, ftype = "gas", slice_info = None):
    unit_system = registry.ds.unit_system
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
        return np.sqrt(3.0 * np.pi / (16.0 * G * data[ftype, "density"]))

    registry.add_field((ftype, "dynamical_time"),
                       function=_dynamical_time,
                       units=unit_system["time"])

    def _jeans_mass(field, data):
        MJ_constant = (((5.0 * kboltz) / (G * mh)) ** (1.5)) * \
          (3.0 / (4.0 * np.pi)) ** (0.5)
        u = (MJ_constant * \
             ((data[ftype, "temperature"] /
               data[ftype, "mean_molecular_weight"])**(1.5)) * \
             (data[ftype, "density"]**(-0.5)))
        return u

    registry.add_field((ftype, "jeans_mass"),
                       function=_jeans_mass,
                       units=unit_system["mass"])

    def _chandra_emissivity(field, data):
        logT0 = np.log10(data[ftype, "temperature"].to_ndarray().astype(np.float64)) - 7
        # we get rid of the units here since this is a fit and not an
        # analytical expression
        return data.ds.arr(data[ftype, "number_density"].to_ndarray().astype(np.float64)**2
                           * (10**(- 0.0103 * logT0**8 + 0.0417 * logT0**7
                                   - 0.0636 * logT0**6 + 0.1149 * logT0**5
                                   - 0.3151 * logT0**4 + 0.6655 * logT0**3
                                   - 1.1256 * logT0**2 + 1.0026 * logT0**1
                                   - 0.6984 * logT0)
                             + data[ftype, "metallicity"].to_ndarray() *
                             10**(  0.0305 * logT0**11 - 0.0045 * logT0**10
                                    - 0.3620 * logT0**9  + 0.0513 * logT0**8
                                    + 1.6669 * logT0**7  - 0.3854 * logT0**6
                                    - 3.3604 * logT0**5  + 0.4728 * logT0**4
                                    + 4.5774 * logT0**3  - 2.3661 * logT0**2
                                    - 1.6667 * logT0**1  - 0.2193 * logT0)),
                           "") # add correct units here

    registry.add_field((ftype, "chandra_emissivity"),
                       function=_chandra_emissivity,
                       units="") # add correct units here

    def _emission_measure(field, data):
        if data.has_field_parameter("X_H"):
            X_H = data.get_field_parameter("X_H")
        else:
            X_H = 0.76
        nenh = data["density"]*data["density"]
        nenh /= mh*mh
        nenh *= 0.5*(1.+X_H)*X_H*data["cell_volume"]
        return nenh
    
    registry.add_field((ftype, "emission_measure"),
                       function=_emission_measure,
                       units=unit_system["number_density"])

    def _xray_emissivity(field, data):
        # old scaling coefficient was 2.168e60
        return data.ds.arr(data[ftype, "density"].to_ndarray().astype(np.float64)**2
                           * data[ftype, "temperature"].to_ndarray()**0.5,
                           "") # add correct units here

    registry.add_field((ftype, "xray_emissivity"),
                       function=_xray_emissivity,
                       units="") # add correct units here

    def _mazzotta_weighting(field, data):
        # Spectroscopic-like weighting field for galaxy clusters
        # Only useful as a weight_field for temperature, metallicity, velocity
        return data["density"]*data["density"]*data["kT"]**-0.25/mh/mh

    registry.add_field((ftype,"mazzotta_weighting"),
                       function=_mazzotta_weighting,
                       units="keV**-0.25*cm**-6")

    def _sz_kinetic(field, data):
        scale = 0.88 * sigma_thompson / mh / clight
        vel_axis = data.get_field_parameter("axis")
        if vel_axis > 2:
            raise NeedsParameter(["axis"])
        vel = data[ftype, "velocity_%s" % ({0: "x", 1: "y", 2: "z"}[vel_axis])]
        # minus sign is because radial velocity is WRT viewer
        # See issue #1225
        return -scale * vel * data[ftype, "density"]

    registry.add_field((ftype, "sz_kinetic"),
                       function=_sz_kinetic,
                       units=unit_system["length"]**-1,
                       validators=[
                           ValidateParameter("axis", {'axis': [0, 1, 2]})])

    def _szy(field, data):
        scale = 0.88 / mh * kboltz / (me * clight*clight) * sigma_thompson
        return scale * data[ftype, "density"] * data[ftype, "temperature"]

    registry.add_field((ftype, "szy"),
                       function=_szy,
                       units=unit_system["length"]**-1)
