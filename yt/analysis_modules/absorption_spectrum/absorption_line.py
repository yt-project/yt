"""
Absorption line generating functions.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as na

def voigt(a,u):
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
            U must be an array, a should not be. Also this procedure is only valid
            for the region a<1.0, u<4.0 or a<1.8(u+1), u>4, which should be most 
            astrophysical conditions (see the article below for further comments
    PROCEDURE:
            Follows procedure in Armstrong JQSRT 7, 85 (1967)
            also the same as the intrinsic in the previous version of IDL
    MODIFICATION HISTORY:
            J. Murthy, Mar 1990 (adapted from the FORTRAN program of Armstrong)
                      Sep 1990 (better overflow checking)
    """
    x = na.asarray(u).astype(na.float64)
    y = na.asarray(a).astype(na.float64)

    w = na.array([0.462243670,   0.286675505,   0.109017206, 
                  0.0248105209,  0.00324377334, 0.000228338636, 
                  7.80255648e-6, 1.08606937e-7, 4.39934099e-10, 
                  2.22939365e-13])

    t = na.array([0.245340708, 0.737473729, 1.23407622, 1.73853771, 
                  2.25497400,  2.78880606,  3.34785457, 3.94476404, 
                  4.60368245,  5.38748089])

    # Hummer's Chebyshev Coefficients
    c = ( 0.1999999999972224, -0.1840000000029998,   0.1558399999965025, 
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
    u1 = na.exp(-x * x + y2) * na.cos(2. * x * y)

    # Clenshaw's Algorithm
    bno1 = na.zeros(x.shape)
    bno2 = na.zeros(x.shape)
    x1 = na.clip((x / 5.), -na.inf, 1.)
    coef = 4. * x1 * x1 - 2.
    for i in range(33, -1, -1):
        bn = coef * bno1 - bno2 + c[i]
        bno2 = na.copy(bno1)
        bno1 = na.copy(bn)

    f = x1 * (bn - bno2)
    dno1 = 1. - 2. * x * f
    dno2 = f

    q = na.abs(x) > 5
    if q.any():
        x14 = na.power(na.clip(x[q], -na.inf, 500.),  14)
        x12 = na.power(na.clip(x[q], -na.inf, 1000.), 12)
        x10 = na.power(na.clip(x[q], -na.inf, 5000.), 10)
        x8  = na.power(na.clip(x[q], -na.inf, 50000.), 8)
        x6  = na.power(na.clip(x[q], -na.inf, 1.e6),   6)
        x4  = na.power(na.clip(x[q], -na.inf, 1.e9),   4)
        x2  = na.power(na.clip(x[q], -na.inf, 1.e18),  2)
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
                g = dn.astype(na.float64) * yn
                funct = funct + q * g
                if na.max(na.abs(g / funct)) <= 1.e-8: break

    k1 = u1 - 1.12837917 * funct
    k1 = k1.astype(na.float64).clip(0)
    return k1

def tau_profile(lam0, fval, gamma, vkms, column_density, 
                deltav=None, delta_lambda=None,
                lambda_bins=None, n_lambda=12000, dlambda=0.01):
    """
    Create an optical depth vs. wavelength profile for an 
    absorption line using a voigt profile.
    :param lam0 (float): central wavelength (angstroms).
    :param fval (float): f-value.
    :param gamma (float): gamma value.
    :param vkms (float): doppler b-parameter.
    :param column_density (float): column density (cm^-2).
    :param deltav (float): velocity offset from lam0 (km/s).
    Default: None (no shift).
    :param delta_lambda (float): wavelength offset in angstroms.
    Default: None (no shift).
    :param lambda_bins (array): array of wavelengths in angstroms.
    Default: None
    :param n_lambda (float): size of lambda bins to create 
    array if lambda_bins is None.  Default: 12000
    :param dlambda (float): lambda bin width if lambda_bins is 
    None. Default: 0.01
    """

    ## constants
    me = 1.6726231e-24 / 1836.        # grams mass electron 
    e = 4.8032e-10                    # esu 
    c = 2.99792456e5                  # km/s
    ccgs = c * 1.e5                   # cm/s 

    ## shift lam0 by deltav
    if deltav is not None:
        lam1 = lam0 * (1 + deltav / c)
    elif delta_lambda is not None:
        lam1 = lam0 + delta_lambda
    else:
        lam1 = lam0

    ## conversions
    vdop = vkms * 1.e5                # in cm/s
    lam0cgs = lam0 / 1.e8             # rest wavelength in cm
    lam1cgs = lam1 / 1.e8             # line wavelength in cm
    nu1 = ccgs / lam1cgs              # line freq in Hz
    nudop = vdop / ccgs * nu1         # doppler width in Hz
    lamdop = vdop / ccgs * lam1       # doppler width in Ang

    ## create wavelength
    if lambda_bins is None:
        lambda_bins = lam1 + \
            na.arange(n_lambda, dtype=na.float) * dlambda - \
            n_lambda * dlambda / 2    # wavelength vector (angstroms)
    nua = ccgs / (lambda_bins / 1.e8) # frequency vector (Hz)

    ## tau_0
    tau_X = na.sqrt(na.pi) * e**2 / (me * ccgs) * \
        column_density * fval / vdop
    tau1 = tau_X * lam1cgs
    tau0 = tau_X * lam0cgs

    # dimensionless frequency offset in units of doppler freq
    x = (nua - nu1) / nudop
    a = gamma / (4 * na.pi * nudop)   # damping parameter 
    phi = voigt(a, x)                 # profile
    tauphi = tau0 * phi               # profile scaled with tau0

    return (lambda_bins, tauphi)
