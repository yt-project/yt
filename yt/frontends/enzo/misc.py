import numpy as np

from yt.utilities.physical_ratios import (
    boltzmann_constant_erg_per_K,
    cm_per_mpc,
    mass_hydrogen_grams,
    newton_cgs,
    rho_crit_g_cm3_h2,
)


def cosmology_get_units(
    hubble_constant, omega_matter, box_size, initial_redshift, current_redshift
):
    """
    Return a dict of Enzo cosmological unit conversions.
    """

    zp1 = 1.0 + current_redshift
    zip1 = 1.0 + initial_redshift

    k = {}
    # For better agreement with values calculated by Enzo,
    # adopt the exact constants that are used there.

    time_scaling = np.sqrt(1 / (4 * np.pi * newton_cgs * rho_crit_g_cm3_h2))
    vel_scaling = cm_per_mpc / time_scaling
    temp_scaling = mass_hydrogen_grams / boltzmann_constant_erg_per_K * vel_scaling**2

    k["utim"] = time_scaling / np.sqrt(omega_matter) / hubble_constant / zip1**1.5
    k["urho"] = rho_crit_g_cm3_h2 * omega_matter * hubble_constant**2 * zp1**3
    k["uxyz"] = cm_per_mpc * box_size / hubble_constant / zp1
    k["uaye"] = 1.0 / zip1
    k["uvel"] = vel_scaling * box_size * np.sqrt(omega_matter) * np.sqrt(zip1)
    k["utem"] = temp_scaling * (box_size**2) * omega_matter * zip1
    k["aye"] = zip1 / zp1
    return k
