"""
Cosmology calculator based on http://www.kempner.net/cosmic.php.  Also,
conversion functions between time and redshift pulled from Enzo.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
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

import numpy as np

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
             return (self.HubbleDistance() / np.sqrt(self.OmegaCurvatureNow) * 
                     np.sinh(np.sqrt(self.OmegaCurvatureNow) * 
                          self.ComovingRadialDistance(z_i,z_f) / 
                          self.HubbleDistance()))
         elif (self.OmegaCurvatureNow < 0):
             return (self.HubbleDistance() / np.sqrt(np.fabs(self.OmegaCurvatureNow)) * 
                     sin(np.sqrt(np.fabs(self.OmegaCurvatureNow)) * 
                         self.ComovingRadialDistance(z_i,z_f) / self.HubbleDistance()))
         else:
             return self.ComovingRadialDistance(z_i,z_f)

    def ComovingVolume(self,z_i,z_f):
        if (self.OmegaCurvatureNow > 0):
             return (2 * np.pi * np.power(self.HubbleDistance(), 3) / self.OmegaCurvatureNow * 
                     (self.ComovingTransverseDistance(z_i,z_f) / self.HubbleDistance() * 
                      np.sqrt(1 + self.OmegaCurvatureNow * 
                           sqr(self.ComovingTransverseDistance(z_i,z_f) / 
                               self.HubbleDistance())) - 
                      anp.sinh(np.fabs(self.OmegaCurvatureNow) * 
                            self.ComovingTransverseDistance(z_i,z_f) / 
                            self.HubbleDistance()) / np.sqrt(self.OmegaCurvatureNow)) / 1e9)
        elif (self.OmegaCurvatureNow < 0):
             return (2 * np.pi * np.power(self.HubbleDistance(), 3) / 
                     np.fabs(self.OmegaCurvatureNow) * 
                     (self.ComovingTransverseDistance(z_i,z_f) / self.HubbleDistance() * 
                      np.sqrt(1 + self.OmegaCurvatureNow * 
                           sqr(self.ComovingTransverseDistance(z_i,z_f) / 
                               self.HubbleDistance())) - 
                      asin(np.fabs(self.OmegaCurvatureNow) * 
                           self.ComovingTransverseDistance(z_i,z_f) / 
                           self.HubbleDistance()) / 
                      np.sqrt(np.fabs(self.OmegaCurvatureNow))) / 1e9)
        else:
             return (4 * np.pi * np.power(self.ComovingTransverseDistance(z_i,z_f), 3) / 
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
        return (self.AngularDiameterDistance(z_i,z_f) / 648. * np.pi)

    def CriticalDensity(self,z):
        return (3.0 / 8.0 / np.pi * sqr(self.HubbleConstantNow / kmPerMpc) / G *
                (self.OmegaLambdaNow + ((1 + z)**3.0) * self.OmegaMatterNow))

    def AgeIntegrand(self,z):
        return (1 / (z + 1) / self.ExpansionFactor(z))

    def ExpansionFactor(self,z):
        return np.sqrt(self.OmegaMatterNow * ((1 + z)**3.0) + 
                    self.OmegaCurvatureNow * np.sqrt(1 + z) + 
                    self.OmegaLambdaNow)

    def InverseExpansionFactor(self,z):
        return 1 / self.ExpansionFactor(z)

    def PathLengthFunction(self,z):
        return ((1 + z)**2) * self.InverseExpansionFactor(z)

    def PathLength(self, z_i, z_f):
        return romberg(self.PathLengthFunction, z_i, z_f)

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

class EnzoCosmology(object):
    def __init__(self, HubbleConstantNow = 71.0,
                 OmegaMatterNow = 0.27,
                 OmegaLambdaNow = 0.73,
                 OmegaCurvatureNow = 0.0,
                 InitialRedshift = 99.0):
        self.HubbleConstantNow = HubbleConstantNow
        self.OmegaMatterNow = OmegaMatterNow
        self.OmegaLambdaNow = OmegaLambdaNow
        self.OmegaCurvatureNow = OmegaCurvatureNow
        self.InitialRedshift = InitialRedshift
        self.InitialTime = self.ComputeTimeFromRedshift(self.InitialRedshift)
        self.TimeUnits = self.ComputeTimeUnits()

    def ComputeTimeUnits(self):
        """
        Taken from CosmologyGetUnits.C in Enzo.
        """
        # Changed 2.52e17 to 2.52e19 because H_0 is in km/s/Mpc, 
        # instead of 100 km/s/Mpc.
        return 2.52e19 / np.sqrt(self.OmegaMatterNow) / \
            self.HubbleConstantNow / np.power(1 + self.InitialRedshift,1.5)

    def ComputeRedshiftFromTime(self,time):
        """
        Compute the redshift from time after the big bang.  This is based on
        Enzo's CosmologyComputeExpansionFactor.C, but altered to use physical
        units.
        """

        OmegaCurvatureNow = 1.0 - self.OmegaMatterNow - self.OmegaLambdaNow

        OMEGA_TOLERANCE = 1e-5
        ETA_TOLERANCE = 1.0e-10

        # Convert the time to Time * H0.
 
        TimeHubble0 = time * self.HubbleConstantNow / kmPerMpc
 
        # 1) For a flat universe with OmegaMatterNow = 1, it's easy.
 
        if ((np.fabs(self.OmegaMatterNow-1) < OMEGA_TOLERANCE) and
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            a = np.power(time/self.InitialTime,2.0/3.0)
 
        # 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
        #    Actually, this is a little tricky since we must solve an equation
        #    of the form eta - np.sinh(eta) + x = 0..
 
        if ((self.OmegaMatterNow < 1) and 
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            x = 2*TimeHubble0*np.power(1.0 - self.OmegaMatterNow, 1.5) / \
                self.OmegaMatterNow;
 
            # Compute eta in a three step process, first from a third-order
            # Taylor expansion of the formula above, then use that in a fifth-order
            # approximation.  Then finally, iterate on the formula itself, solving for
            # eta.  This works well because parts 1 & 2 are an excellent approximation
            # when x is small and part 3 converges quickly when x is large. 
 
            eta = np.power(6*x,1.0/3.0)                # part 1
            eta = np.power(120*x/(20+eta*eta),1.0/3.0) # part 2
            for i in range(40):                      # part 3
                eta_old = eta
                eta = np.arcsinh(eta + x)
                if (np.fabs(eta-eta_old) < ETA_TOLERANCE): 
                    break
                if (i == 39):
                    print "No convergence after %d iterations." % i
 
            # Now use eta to compute the expansion factor (eq. 13-10, part 2).
 
            a = self.OmegaMatterNow/(2.0*(1.0 - self.OmegaMatterNow))*\
                (np.cosh(eta) - 1.0)

        # 3) For OmegaMatterNow > 1 and OmegaLambdaNow == 0, use sin/cos.
        #    Easy, but skip it for now.
 
        if ((self.OmegaMatterNow > 1) and 
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            print "Never implemented in Enzo, not implemented here."
            return 0
 
        # 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20.
 
        if ((np.fabs(OmegaCurvatureNow) < OMEGA_TOLERANCE) and
            (self.OmegaLambdaNow > OMEGA_TOLERANCE)):
            a = np.power(self.OmegaMatterNow / (1 - self.OmegaMatterNow),1.0/3.0) * \
                np.power(np.sinh(1.5 * np.sqrt(1.0 - self.OmegaMatterNow)*\
                                     TimeHubble0),2.0/3.0)


        redshift = (1.0/a) - 1.0

        return redshift

    def ComputeTimeFromRedshift(self,z):
        """
        Compute the time from redshift.  This is based on Enzo's
        CosmologyComputeTimeFromRedshift.C, but altered to use physical units.
        """
        OmegaCurvatureNow = 1.0 - self.OmegaMatterNow - self.OmegaLambdaNow
 
        # 1) For a flat universe with OmegaMatterNow = 1, things are easy.
 
        if ((self.OmegaMatterNow == 1.0) and (self.OmegaLambdaNow == 0.0)):
            TimeHubble0 = 2.0/3.0/np.power(1+z,1.5)
 
        # 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
 
        if ((self.OmegaMatterNow < 1) and (self.OmegaLambdaNow == 0)):
            eta = np.arccosh(1 + 2*(1-self.OmegaMatterNow)/self.OmegaMatterNow/(1+z))
            TimeHubble0 = self.OmegaMatterNow/(2*np.power(1.0-self.OmegaMatterNow, 1.5))*\
                (np.sinh(eta) - eta)
 
        # 3) For OmegaMatterNow > 1 and OmegaLambdaNow == 0, use sin/cos.
 
        if ((self.OmegaMatterNow > 1) and (self.OmegaLambdaNow == 0)):
            eta = np.acos(1 - 2*(1-self.OmegaMatterNow)/self.OmegaMatterNow/(1+z))
            TimeHubble0 = self.OmegaMatterNow/(2*np.power(1.0-self.OmegaMatterNow, 1.5))*\
                (eta - np.sin(eta))
 
        # 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20.
 
        if ((np.fabs(OmegaCurvatureNow) < 1.0e-3) and (self.OmegaLambdaNow != 0)):
            TimeHubble0 = 2.0/3.0/np.sqrt(1-self.OmegaMatterNow)*\
                np.arcsinh(np.sqrt((1-self.OmegaMatterNow)/self.OmegaMatterNow)/ \
                               np.power(1+z,1.5))
  
        # Now convert from Time * H0 to time.
  
        time = TimeHubble0 / (self.HubbleConstantNow/kmPerMpc)
    
        return time
