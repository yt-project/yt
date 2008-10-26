"""
Cosmology calculator based on http://www.kempner.net/cosmic.php.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Britton Smith.  All Rights Reserved.

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

c_kms = 2.99792458e5 # c in km/s
G = 6.67259e-8 # cgs
kmPerMpc = 3.08567758e19

class Cosmology(object):
    def __init__(self, HubbleConstantNow = 71.0,
                 OmegaMatterNow = 0.27,
                 OmegaLambdaNow = 0.73,
                 OmegaCurvatureNow = 0.0):
        self.HubbleConstantNow = HubbleConstantNow
        self.OmegaMatterNow = OmegaMatterNow
        self.OmegaLambdaNow = OmegaLambdaNow
        self.OmegaCurvatureNow = OmegaCurvatureNow

    def HubbleDistance(self):
        return (c_kms / self.HubbleConstantNow)

    def ComovingRadialDistance(self,z_i,z_f):
        return self.HubbleDistance() * \
            romberg(self.InverseExpansionFactor,z_i,z_f)

    def ComovingTransverseDistance(self,z_i,z_f):
         if (self.OmegaCurvatureNow > 0):
             return (self.HubbleDistance() / na.sqrt(self.OmegaCurvatureNow) * 
                     na.sinh(na.sqrt(self.OmegaCurvatureNow) * 
                          self.ComovingRadialDistance(z_i,z_f) / 
                          self.HubbleDistance()))
         elif (self.OmegaCurvatureNow < 0):
             return (self.HubbleDistance() / na.sqrt(na.fabs(self.OmegaCurvatureNow)) * 
                     sin(na.sqrt(na.fabs(self.OmegaCurvatureNow)) * 
                         self.ComovingRadialDistance(z_i,z_f) / self.HubbleDistance()))
         else:
             return self.ComovingRadialDistance(z_i,z_f)

    def ComovingVolume(self,z_i,z_f):
        if (self.OmegaCurvatureNow > 0):
             return (2 * na.pi * na.power(self.HubbleDistance(), 3) / self.OmegaCurvatureNow * 
                     (self.ComovingTransverseDistance(z_i,z_f) / self.HubbleDistance() * 
                      na.sqrt(1 + self.OmegaCurvatureNow * 
                           sqr(self.ComovingTransverseDistance(z_i,z_f) / 
                               self.HubbleDistance())) - 
                      ana.sinh(na.fabs(self.OmegaCurvatureNow) * 
                            self.ComovingTransverseDistance(z_i,z_f) / 
                            self.HubbleDistance()) / na.sqrt(self.OmegaCurvatureNow)) / 1e9)
        elif (self.OmegaCurvatureNow < 0):
             return (2 * na.pi * na.power(self.HubbleDistance(), 3) / 
                     na.fabs(self.OmegaCurvatureNow) * 
                     (self.ComovingTransverseDistance(z_i,z_f) / self.HubbleDistance() * 
                      na.sqrt(1 + self.OmegaCurvatureNow * 
                           sqr(self.ComovingTransverseDistance(z_i,z_f) / 
                               self.HubbleDistance())) - 
                      asin(na.fabs(self.OmegaCurvatureNow) * 
                           self.ComovingTransverseDistance(z_i,z_f) / 
                           self.HubbleDistance()) / 
                      na.sqrt(na.fabs(self.OmegaCurvatureNow))) / 1e9)
        else:
             return (4 * na.pi * na.power(self.ComovingTransverseDistance(z_i,z_f), 3) / 
                     3 / 1e9)

    def AngularDiameterDistance(self,z_i,z_f):
        return (self.ComovingTransverseDistance(0,z_f) / (1 + z_f) - 
                self.ComovingTransverseDistance(0,z_i) / (1 + z_i))

    def LuminosityDistance(self,z_i,z_f):
        return (self.ComovingTransverseDistance(0,z_f) * (1 + z_f) - 
                self.ComovingTransverseDistance(0,z_i) * (1 + z_i))

    def LookbackTime(self,z_i,z_f):
        return (romberg(self.AgeIntegrand,z_i,z_f) / self.HubbleConstantNow * kmPerMpc)

    def UniverseAge(self,z):
        return (romberg(self.AgeIntegrand,z,1000) / self.HubbleConstantNow * kmPerMpc)

    def AngularScale_1arcsec_kpc(self,z_i,z_f):
        return (self.AngularDiameterDistance(z_i,z_f) / 648. * na.pi)

    def CriticalDensity(self,z):
        return (3.0 / 8.0 / na.pi * sqr(self.HubbleConstantNow / kmPerMpc) / G *
                (self.OmegaLambdaNow + ((1 + z)**3.0) * self.OmegaMatterNow))

    def AgeIntegrand(self,z):
        return (1 / (z + 1) / self.ExpansionFactor(z))

    def ExpansionFactor(self,z):
        return na.sqrt(self.OmegaMatterNow * ((1 + z)**3.0) + 
                    self.OmegaCurvatureNow * na.sqrt(1 + z) + 
                    self.OmegaLambdaNow)

    def InverseExpansionFactor(self,z):
        return 1 / self.ExpansionFactor(z)

def sqr(x):
    return (x**2.0)

def romberg(f, a, b, eps=1e-8):
    """Approximate the definite integral of f from a to b by Romberg's method.
    eps is the desired accuracy."""
    R = [[0.5 * (b - a) * (f(a) + f(b))]]  # R[0][0]
    n = 1
    while True:
        h = float(b - a) / 2 ** n
        R.append([None] * (n + 1))  # Add an empty row.
        # for proper limits
        R[n][0] = 0.5*R[n-1][0] + h*sum(f(a+(2*k-1)*h) for k in xrange(1, 2**(n-1)+1))
        for m in xrange(1, n+1):
            R[n][m] = R[n][m-1] + (R[n][m-1] - R[n-1][m-1]) / (4 ** m - 1)
        if abs(R[n][n-1] - R[n][n]) < eps:
            return R[n][n]
        n += 1
