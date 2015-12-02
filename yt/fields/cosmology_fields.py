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

from .derived_field import \
    ValidateParameter
from .field_exceptions import \
    NeedsConfiguration, \
    NeedsParameter
from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.physical_constants import \
    speed_of_light_cgs

@register_field_plugin
def setup_cosmology_fields(registry, ftype = "gas", slice_info = None):
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

    def _matter_density(field, data):
        return data[ftype, "density"] + \
          data[ftype, "dark_matter_density"]

    registry.add_field((ftype, "matter_density"),
                       function=_matter_density,
                       units=unit_system["density"])

    def _matter_mass(field, data):
        return data[ftype, "matter_density"] * data["index", "cell_volume"]

    registry.add_field((ftype, "matter_mass"),
                       function=_matter_mass,
                       units=unit_system["mass"])

    # rho_total / rho_cr(z).
    def _overdensity(field, data):
        if not hasattr(data.ds, "cosmological_simulation") or \
          not data.ds.cosmological_simulation:
            raise NeedsConfiguration("cosmological_simulation", 1)
        co = data.ds.cosmology
        return data[ftype, "matter_density"] / \
          co.critical_density(data.ds.current_redshift)

    registry.add_field((ftype, "overdensity"),
                       function=_overdensity,
                       units="")

    # rho_baryon / <rho_baryon>
    def _baryon_overdensity(field, data):
        if not hasattr(data.ds, "cosmological_simulation") or \
          not data.ds.cosmological_simulation:
            raise NeedsConfiguration("cosmological_simulation", 1)
        omega_baryon = data.get_field_parameter("omega_baryon")
        if omega_baryon is None:
            raise NeedsParameter("omega_baryon")
        co = data.ds.cosmology
        # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
        # mean matter density(z) ~ omega_matter * (1 + z)^3
        return data[ftype, "density"] / omega_baryon / co.critical_density(0.0) / \
          (1.0 + data.ds.current_redshift)**3

    registry.add_field((ftype, "baryon_overdensity"),
                       function=_baryon_overdensity,
                       units="",
                       validators=[ValidateParameter("omega_baryon")])

    # rho_matter / <rho_matter>
    def _matter_overdensity(field, data):
        if not hasattr(data.ds, "cosmological_simulation") or \
          not data.ds.cosmological_simulation:
            raise NeedsConfiguration("cosmological_simulation", 1)
        co = data.ds.cosmology
        # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
        # mean density(z) ~ omega_matter * (1 + z)^3
        return data[ftype, "matter_density"] / data.ds.omega_matter / \
          co.critical_density(0.0) / \
          (1.0 + data.ds.current_redshift)**3

    registry.add_field((ftype, "matter_overdensity"),
                       function=_matter_overdensity,
                       units="")

    # r / r_vir
    def _virial_radius_fraction(field, data):
        virial_radius = data.get_field_parameter("virial_radius")
        return data["radius"] / virial_radius

    registry.add_field(("index", "virial_radius_fraction"),
                       function=_virial_radius_fraction,
                       validators=[ValidateParameter("virial_radius")],
                       units="")

    # Weak lensing convergence.
    # Eqn 4 of Metzler, White, & Loken (2001, ApJ, 547, 560).
    # This needs to be checked for accuracy.
    def _weak_lensing_convergence(field, data):
        if not hasattr(data.ds, "cosmological_simulation") or \
          not data.ds.cosmological_simulation:
            raise NeedsConfiguration("cosmological_simulation", 1)
        co = data.ds.cosmology
        observer_redshift = data.get_field_parameter('observer_redshift')
        source_redshift = data.get_field_parameter('source_redshift')

        # observer to lens
        dl = co.angular_diameter_distance(observer_redshift, data.ds.current_redshift)
        # observer to source
        ds = co.angular_diameter_distance(observer_redshift, source_redshift)
        # lens to source
        dls = co.angular_diameter_distance(data.ds.current_redshift, source_redshift)

        # removed the factor of 1 / a to account for the fact that we are projecting
        # with a proper distance.
        return (1.5 * (co.hubble_constant / speed_of_light_cgs)**2 * (dl * dls / ds) * \
          data[ftype, "matter_overdensity"]).in_units("1/cm")

    registry.add_field((ftype, "weak_lensing_convergence"),
                       function=_weak_lensing_convergence,
                       units=unit_system["length"]**-1,
        validators=[ValidateParameter("observer_redshift"),
                    ValidateParameter("source_redshift")])
