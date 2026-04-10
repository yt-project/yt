import numpy as np
from numpy.testing import assert_allclose

from yt.utilities import physical_constants as pc
from yt.utilities import physical_ratios as pr


def test_physical_ratios_are_self_consistent():
    assert_allclose(pr.mass_hydrogen_grams, 1.007947 * pr.amu_grams)
    assert_allclose(pr.cm_per_mpc * pr.mpc_per_cm, 1.0)
    assert_allclose(pr.sec_per_day / pr.sec_per_hr, 24.0)
    assert_allclose(pr.sec_per_hr / pr.sec_per_min, 60.0)
    assert_allclose(pr.mass_earth_grams * 328900.56, pr.mass_sun_grams)
    assert_allclose(pr.planck_time_s * pr.speed_of_light_cm_per_s, pr.planck_length_cm)
    assert_allclose(pr._primordial_mass_fraction["H"], 0.76)
    assert_allclose(pr._primordial_mass_fraction["He"], 0.24)


def test_physical_constants_aliases_and_em_relations():
    assert pc.mass_electron is pc.me
    assert pc.mass_hydrogen is pc.mp
    assert pc.mass_hydrogen is pc.mh
    assert pc.speed_of_light is pc.c
    assert pc.speed_of_light is pc.clight
    assert pc.planck_constant is pc.hcgs

    assert_allclose(pc.mu_0.to_value("N/A**2"), 4.0e-7 * np.pi)
    assert_allclose(
        (pc.eps_0 * pc.clight**2 * pc.mu_0).to_value("dimensionless"),
        1.0,
    )
    assert_allclose(pc.hbar.to_value("erg*s"), 0.5 * pc.hcgs.to_value("erg*s") / np.pi)
    assert_allclose(pc.Na.to_value("1/g"), 1.0 / pr.amu_grams)
