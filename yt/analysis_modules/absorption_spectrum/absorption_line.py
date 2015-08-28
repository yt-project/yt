"""
Absorption line generating functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.utilities.physical_constants import \
    charge_proton_cgs, \
    mass_electron_cgs, \
    speed_of_light_cgs
from yt.utilities.on_demand_imports import _scipy, NotAModule

special = _scipy.special
tau_factor = None
_cs = None


def voigt_scipy(a, u):
    x = np.asarray(u).astype(np.float64)
    y = np.asarray(a).astype(np.float64)
    return special.wofz(x + 1j * y).real


def voigt_old(a, u):
    """
    NAME:
        VOIGT
    PURPOSE:
        Implementation of Voigt function
    CATEGORY:
            Math
    CALLING SEQUENCE:
            voigt=Voigt(a,u)
    INPUTS:
            A = Voigt "A" parameter.
            U = Frequency in units of the Doppler frequency.

            The line profile "Phi(v)", the doppler width
            "Delv", the voigt parameter "a", and the frequency "u"
            are given by:

            Phi(v) =  Voigt(a,u)/[ Delv * sqrt(pi) ]
            Delv   =  Vo/c * sqrt[ 2kT/m ]
            u      =  V - Vo / Delv
            a      =  GAMMA / [ Delv * 4pi ]
            Gamma  =  Gu + Gl + 2*Vcol
            "Gu" and "Gl" are the widths of the upper and lower states
            "Vcol" is the collisions per unit time
            "Vo" is the line center frequency

    OUTPUTS:
            An array of the same type as u
    RESTRICTIONS:
            U must be an array, a should not be. Also this procedure is only
            valid for the region a<1.0, u<4.0 or a<1.8(u+1), u>4, which should
            be most astrophysical conditions (see the article below for further
            comments
    PROCEDURE:
            Follows procedure in Armstrong JQSRT 7, 85 (1967)
            also the same as the intrinsic in the previous version of IDL
    MODIFICATION HISTORY:
            J. Murthy, Mar 1990 (adapted from the FORTRAN program of Armstrong)
                      Sep 1990 (better overflow checking)
    """
    x = np.asarray(u).astype(np.float64)
    y = np.asarray(a).astype(np.float64)

    # Hummer's Chebyshev Coefficients
    c = (0.1999999999972224, -0.1840000000029998,   0.1558399999965025,
         -0.1216640000043988,  0.0877081599940391,  -0.0585141248086907,
         0.0362157301623914, -0.0208497654398036,   0.0111960116346270,
         -0.56231896167109e-2, 0.26487634172265e-2, -0.11732670757704e-2,
         0.4899519978088e-3, -0.1933630801528e-3,   0.722877446788e-4,
         -0.256555124979e-4,   0.86620736841e-5,    -0.27876379719e-5,
         0.8566873627e-6,    -0.2518433784e-6,      0.709360221e-7,
         -0.191732257e-7,      0.49801256e-8,       -0.12447734e-8,
         0.2997777e-9,       -0.696450e-10,         0.156262e-10,
         -0.33897e-11,         0.7116e-12,          -0.1447e-12,
         0.285e-13,          -0.55e-14,             0.10e-14,
         -0.2e-15)

    y2 = y * y

    # limits are y<1.,  x<4 or y<1.8(x+1),  x>4 (no checking performed)
    u1 = np.exp(-x * x + y2) * np.cos(2. * x * y)

    # Clenshaw's Algorithm
    bno1 = np.zeros(x.shape)
    bno2 = np.zeros(x.shape)
    x1 = np.clip((x / 5.), -np.inf, 1.)
    coef = 4. * x1 * x1 - 2.
    for i in range(33, -1, -1):
        bn = coef * bno1 - bno2 + c[i]
        bno2 = np.copy(bno1)
        bno1 = np.copy(bn)

    f = x1 * (bn - bno2)
    dno1 = 1. - 2. * x * f
    dno2 = f

    q = np.abs(x) > 5
    if q.any():
        x14 = np.power(np.clip(x[q], -np.inf, 500.),  14)
        x12 = np.power(np.clip(x[q], -np.inf, 1000.), 12)
        x10 = np.power(np.clip(x[q], -np.inf, 5000.), 10)
        x8 = np.power(np.clip(x[q], -np.inf, 50000.), 8)
        x6 = np.power(np.clip(x[q], -np.inf, 1.e6),   6)
        x4 = np.power(np.clip(x[q], -np.inf, 1.e9),   4)
        x2 = np.power(np.clip(x[q], -np.inf, 1.e18),  2)
        dno1[q] = -(0.5 / x2 + 0.75 / x4 + 1.875 / x6 +
                    6.5625 / x8 + 29.53125 / x10 +
                    162.4218 / x12 + 1055.7421 / x14)
        dno2[q] = (1. - dno1[q]) / (2. * x[q])

    funct = y * dno1
    if (y > 1.e-8).any():
        q = 1.0
        yn = y
        for i in range(2, 51):
            dn = (x * dno1 + dno2) * (-2. / i)
            dno2 = dno1
            dno1 = dn
            if (i % 2) == 1:
                q = -q
                yn = yn * y2
                g = dn.astype(np.float64) * yn
                funct = funct + q * g
                if np.max(np.abs(g / funct)) <= 1.e-8:
                    break

    k1 = u1 - 1.12837917 * funct
    k1 = k1.astype(np.float64).clip(0)
    return k1


def tau_profile(lambda_0, f_value, gamma, v_doppler, column_density,
                delta_v=None, delta_lambda=None,
                lambda_bins=None, n_lambda=12000, dlambda=0.01):
    r"""
    Create an optical depth vs. wavelength profile for an
    absorption line using a voigt profile.

    Parameters
    ----------

    lambda_0 : float in angstroms
       central wavelength.
    f_value : float
       absorption line f-value.
    gamma : float
       absorption line gamma value.
    v_doppler : float in cm/s
       doppler b-parameter.
    column_density : float in cm^-2
       column density.
    delta_v : float in cm/s
       velocity offset from lambda_0.
       Default: None (no shift).
    delta_lambda : float in angstroms
        wavelength offset.
        Default: None (no shift).
    lambda_bins : array in angstroms
        wavelength array for line deposition.  If None, one will be
        created using n_lambda and dlambda.
        Default: None.
    n_lambda : int
        size of lambda bins to create if lambda_bins is None.
        Default: 12000.
    dlambda : float in angstroms
        lambda bin width in angstroms if lambda_bins is None.
        Default: 0.01.

    """
    global tau_factor
    if tau_factor is None:
        tau_factor = (
            np.sqrt(np.pi) * charge_proton_cgs ** 2 /
            (mass_electron_cgs * speed_of_light_cgs)
        ).in_cgs().d

    global _cs
    if _cs is None:
        _cs = speed_of_light_cgs.d[()]

    # shift lambda_0 by delta_v
    if delta_v is not None:
        lam1 = lambda_0 * (1 + delta_v / _cs)
    elif delta_lambda is not None:
        lam1 = lambda_0 + delta_lambda
    else:
        lam1 = lambda_0

    # conversions
    nudop = 1e8 * v_doppler / lam1   # doppler width in Hz

    # create wavelength
    if lambda_bins is None:
        lambda_bins = lam1 + \
            np.arange(n_lambda, dtype=np.float) * dlambda - \
            n_lambda * dlambda / 2  # wavelength vector (angstroms)

    # tau_0
    tau_X = tau_factor * column_density * f_value / v_doppler
    tau0 = tau_X * lambda_0 * 1e-8

    # dimensionless frequency offset in units of doppler freq
    x = _cs / v_doppler * (lam1 / lambda_bins - 1.0)
    a = gamma / (4.0 * np.pi * nudop)               # damping parameter
    phi = voigt(a, x)                               # line profile
    tauphi = tau0 * phi              # profile scaled with tau0

    return (lambda_bins, tauphi)

if isinstance(special, NotAModule):
    voigt = voigt_old
else:
    voigt = voigt_scipy
