"""
Test cosmology calculator.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.testing import \
     assert_almost_equal, \
     assert_rel_equal, \
     assert_equal
from yt.units.yt_array import \
     YTArray, \
     YTQuantity
from yt.utilities.cosmology import \
     Cosmology

def z_from_t_analytic(my_time, hubble_constant=0.7,
                      omega_matter=0.3, omega_lambda=0.7):
    """
    Compute the redshift from time after the big bang.  This is based on
    Enzo's CosmologyComputeExpansionFactor.C, but altered to use physical
    units.
    """

    hubble_constant = YTQuantity(hubble_constant, "100*km/s/Mpc")
    omega_curvature = 1.0 - omega_matter - omega_lambda

    OMEGA_TOLERANCE = 1e-5
    ETA_TOLERANCE = 1.0e-10

    # Convert the time to Time * H0.

    if not isinstance(my_time, YTArray):
        my_time = YTArray(my_time, "s")

    t0 = (my_time.in_units("s") *
          hubble_constant.in_units("1/s")).to_ndarray()

    # For a flat universe with omega_matter = 1, it's easy.

    if np.fabs(omega_matter-1) < OMEGA_TOLERANCE and \
      omega_lambda < OMEGA_TOLERANCE:
        a = np.power(1.5 * t0, 2.0/3.0)

    # For omega_matter < 1 and omega_lambda == 0 see
    # Peebles 1993, eq. 13-3, 13-10.
    # Actually, this is a little tricky since we must solve an equation
    # of the form eta - np.sinh(eta) + x = 0..

    elif omega_matter < 1 and omega_lambda < OMEGA_TOLERANCE:
        x = 2*t0*np.power(1.0 - omega_matter, 1.5) / omega_matter

        # Compute eta in a three step process, first from a third-order
        # Taylor expansion of the formula above, then use that in a fifth-order
        # approximation.  Then finally, iterate on the formula itself, solving for
        # eta.  This works well because parts 1 & 2 are an excellent approximation
        # when x is small and part 3 converges quickly when x is large.

        eta = np.power(6*x, 1.0/3.0)                # part 1
        eta = np.power(120*x/(20+eta*eta), 1.0/3.0) # part 2
        mask = np.ones(eta.size, dtype=bool)
        max_iter = 1000
        for i in range(max_iter):                   # part 3
            eta_old = eta[mask]
            eta[mask] = np.arcsinh(eta[mask] + x[mask])
            mask[mask] = np.fabs(eta[mask]-eta_old) >= ETA_TOLERANCE
            if not mask.any():
                break
            if (i == max_iter-1):
                raise RuntimeError(
                    "No convergence after %d iterations." % i)

        # Now use eta to compute the expansion factor (eq. 13-10, part 2).

        a = omega_matter/(2.0*(1.0 - omega_matter))*\
            (np.cosh(eta) - 1.0)

    # For flat universe, with non-zero omega_lambda, see eq. 13-20.

    elif np.fabs(omega_curvature) < OMEGA_TOLERANCE and \
      omega_lambda > OMEGA_TOLERANCE:
        a = np.power(omega_matter / (1 - omega_matter), 1.0/3.0) * \
          np.power(np.sinh(1.5 * np.sqrt(1.0 - omega_matter) * \
                           t0), 2.0/3.0)

    else:
        raise NotImplementedError

    redshift = (1.0/a) - 1.0

    return redshift

def t_from_z_analytic(z, hubble_constant=0.7,
                      omega_matter=0.3, omega_lambda=0.7):
    """
    Compute the age of the Universe from redshift.  This is based on Enzo's
    CosmologyComputeTimeFromRedshift.C, but altered to use physical units.
    """

    hubble_constant = YTQuantity(hubble_constant, "100*km/s/Mpc")
    omega_curvature = 1.0 - omega_matter - omega_lambda

    # For a flat universe with omega_matter = 1, things are easy.

    if omega_matter == 1.0 and omega_lambda == 0.0:
        t0 = 2.0/3.0/np.power(1+z, 1.5)

    # For omega_matter < 1 and omega_lambda == 0 see
    # Peebles 1993, eq. 13-3, 13-10.

    elif omega_matter < 1 and omega_lambda == 0:
        eta = np.arccosh(1 + 2*(1-omega_matter)/omega_matter/(1+z))
        t0 = omega_matter/ \
          (2*np.power(1.0-omega_matter, 1.5))*\
          (np.sinh(eta) - eta)

    # For flat universe, with non-zero omega_lambda, see eq. 13-20.

    elif np.fabs(omega_curvature) < 1.0e-3 and omega_lambda != 0:
        t0 = 2.0/3.0/np.sqrt(1-omega_matter)*\
          np.arcsinh(np.sqrt((1-omega_matter)/omega_matter)/ \
                     np.power(1+z, 1.5))

    else:
        raise NotImplementedError("%s, %s, %s" %
                                  (hubble_constant, omega_matter, omega_lambda))

    # Now convert from Time * H0 to time.

    my_time = t0 / hubble_constant

    return my_time

def test_z_t_roundtrip():
    """
    Make sure t_from_z and z_from_t are consistent.

    """

    co = Cosmology()
    # random sample in log(a) from -6 to 6
    my_random = np.random.RandomState(6132305)
    la = 12 * my_random.random_sample(10000) - 6
    z1 = 1 / np.power(10, la) - 1
    t = co.t_from_z(z1)
    z2 = co.z_from_t(t)
    assert_rel_equal(z1, z2, 4)

def test_z_t_analytic():
    """
    Test z/t conversions against analytic solutions.
    """

    cosmos = (
        {'hubble_constant': 0.7, 'omega_matter': 0.3, 'omega_lambda': 0.7},
        {'hubble_constant': 0.7, 'omega_matter': 1.0, 'omega_lambda': 0.0},
        {'hubble_constant': 0.7, 'omega_matter': 0.3, 'omega_lambda': 0.0},
    )

    for cosmo in cosmos:
        omega_curvature = 1 - cosmo['omega_matter'] - cosmo['omega_lambda']
        co = Cosmology(omega_curvature=omega_curvature, **cosmo)
        # random sample in log(a) from -6 to 6
        my_random = np.random.RandomState(10132324)
        la = 12 * my_random.random_sample(1000) - 6
        z = 1 / np.power(10, la) - 1

        t_an = t_from_z_analytic(z, **cosmo).to('Gyr')
        t_co = co.t_from_z(z).to('Gyr')

        assert_rel_equal(
            t_an, t_co, 4,
            err_msg='t_from_z does not match analytic version for cosmology %s.' % cosmo)

        # random sample in log(t/t0) from -3 to 1
        t0 = np.power(10, 4 * my_random.random_sample(1000) - 3)
        t = (t0 / co.hubble_constant).to('Gyr')

        z_an = z_from_t_analytic(t, **cosmo)
        z_co = co.z_from_t(t)

        # compare scale factors since z approaches 0
        assert_rel_equal(
            1 / (1 + z_an), 1 / (1 + z_co), 5,
            err_msg='z_from_t does not match analytic version for cosmology %s.' % cosmo)

def test_cosmology_calculator_answers():
    """
    Test cosmology calculator functions against previously calculated values.
    """

    cosmos = (
        {'omega_matter': 0.3, 'omega_lambda': 0.7, 'omega_radiation': 0.0},
        {'omega_matter': 1.0, 'omega_lambda': 0.0, 'omega_radiation': 0.0},
        {'omega_matter': 0.3, 'omega_lambda': 0.0, 'omega_radiation': 0.0},
        {'omega_matter': 0.2999, 'omega_lambda': 0.7, 'omega_radiation': 1e-4},
    )

    all_answers = (
        [4282.7494, 1174.71744826, 1174.71744826, 6.79029852, -137.44380354,
         -137.44380354, 9876.12942342, 1113.19802768, 421.60592442, 123.24771803],
        [4282.7494, 661.92695892, 661.92695892, 1.21483926, -135.88421154,
         -135.88421154, 6266.65604101, 627.61305115, 1088.01528882, 197.98989873],
        [4282.7494, 933.24209295, 938.4206189, 48.25005669, -55.51112604,
         -55.51112604, 10177.58171302, 883.64108565, 707.20993773, 159.62455951],
        [4282.7494, 1174.26084215, 1174.26084215, 6.78238355, -137.51190671,
         -137.51190671, 9873.75444152, 1112.76996783, 421.71472595, 123.26361994]
    )

    z   = 1
    z_i = 2
    z_f = 3

    for cosmo, answers in zip(cosmos, all_answers):
        omega_curvature = 1 - cosmo['omega_matter'] - \
          cosmo['omega_lambda'] - cosmo['omega_radiation']
        co = Cosmology(hubble_constant=0.7,
                       omega_curvature=omega_curvature, **cosmo)

        values = [
            co.hubble_distance().to('Mpc'),
            co.comoving_radial_distance(z_i, z_f).to('Mpc'),
            co.comoving_transverse_distance(z_i, z_f).to('Mpc'),
            co.comoving_volume(z_i, z_f).to('Gpc**3'),
            co.angular_diameter_distance(z_i, z_f).to('Mpc'),
            co.angular_scale(z_i, z_f).to('Mpc/radian'),
            co.luminosity_distance(z_i, z_f).to('Mpc'),
            co.lookback_time(z_i, z_f).to('Myr'),
            co.critical_density(z).to('Msun/kpc**3'),
            co.hubble_parameter(z).to('km/s/Mpc')
        ]

        for value, answer in zip(values, answers):
            assert_almost_equal(value.d, answer, 8)

def test_dark_factor():
    """
    Test that dark factor returns same value for when not
    being used and when w_0 = -1 and w_z = 0.
    """

    co = Cosmology(w_0=-1, w_a=0,
                   use_dark_factor=False)

    assert_equal(co.get_dark_factor(0), 1.0)
    co.use_dark_factor = True
    assert_equal(co.get_dark_factor(0), 1.0)
