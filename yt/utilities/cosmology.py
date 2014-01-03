"""
Cosmology calculator based on http://www.kempner.net/cosmic.php.  Also,
conversion functions between time and redshift pulled from Enzo.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.data_objects.yt_array import \
     YTArray, \
     YTQuantity

from yt.utilities.physical_constants import \
    gravitational_constant_cgs, \
    km_per_cm, \
    pc_per_mpc, \
    km_per_pc, \
    speed_of_light_cgs

c_kms = speed_of_light_cgs * km_per_cm # c in km/s
G = gravitational_constant_cgs
kmPerMpc = km_per_pc * pc_per_mpc

class Cosmology(object):
    def __init__(self, hubble_constant = 0.71,
                 omega_matter = 0.27,
                 omega_lambda = 0.73,
                 omega_curvature = 0.0):
        self.hubble_constant = YTQuantity(hubble_constant, "100*km/s/Mpc")
        self.omega_matter = omega_matter
        self.omega_lambda = omega_lambda
        self.omega_curvature = omega_curvature

    def hubble_distance(self):
        return (speed_of_light_cgs / self.hubble_constant).in_cgs()

    def comoving_radial_distance(self, z_i, z_f):
        return (self.hubble_distance() *
                trapzint(self.inverse_expansion_factor, z_i, z_f)).in_cgs()

    def comoving_transverse_distance(self, z_i, z_f):
         if (self.omega_curvature > 0):
             return (self.hubble_distance() / np.sqrt(self.omega_curvature) * 
                     np.sinh(np.sqrt(self.omega_curvature) * 
                          self.comoving_radial_distance(z_i, z_f) /
                          self.hubble_distance())).in_cgs()
         elif (self.omega_curvature < 0):
             return (self.hubble_distance() /
                     np.sqrt(np.fabs(self.omega_curvature)) * 
                     np.sin(np.sqrt(np.fabs(self.omega_curvature)) * 
                            self.comoving_radial_distance(z_i, z_f) /
                            self.hubble_distance())).in_cgs()
         else:
             return self.comoving_radial_distance(z_i, z_f)

    def comoving_volume(self, z_i, z_f):
        if (self.omega_curvature > 0):
             return (2 * np.pi * np.power(self.hubble_distance(), 3) /
                     self.omega_curvature * 
                     (self.comoving_transverse_distance(z_i, z_f) /
                      self.hubble_distance() * 
                      np.sqrt(1 + self.omega_curvature * 
                           sqr(self.comoving_transverse_distance(z_i, z_f) /
                               self.hubble_distance())) - 
                      np.sinh(np.fabs(self.omega_curvature) * 
                            self.comoving_transverse_distance(z_i, z_f) /
                            self.hubble_distance()) /
                            np.sqrt(self.omega_curvature)) / 1e9).in_cgs()
        elif (self.omega_curvature < 0):
             return (2 * np.pi * np.power(self.hubble_distance(), 3) /
                     np.fabs(self.omega_curvature) * 
                     (self.comoving_transverse_distance(z_i, z_f) /
                      self.hubble_distance() * 
                      np.sqrt(1 + self.omega_curvature * 
                           sqr(self.comoving_transverse_distance(z_i, z_f) /
                               self.hubble_distance())) - 
                      np.arcsin(np.fabs(self.omega_curvature) * 
                           self.comoving_transverse_distance(z_i, z_f) /
                           self.hubble_distance()) /
                      np.sqrt(np.fabs(self.omega_curvature))) / 1e9).in_cgs()
        else:
             return (4 * np.pi *
                     np.power(self.comoving_transverse_distance(z_i, z_f), 3) /\
                     3 / 1e9).in_cgs()

    def angular_diameter_distance(self, z_i, z_f):
        return (self.comoving_transverse_distance(0, z_f) / (1 + z_f) - 
                self.comoving_transverse_distance(0, z_i) / (1 + z_i)).in_cgs()

    def luminosity_distance(self, z_i, z_f):
        return (self.comoving_transverse_distance(0, z_f) * (1 + z_f) - 
                self.comoving_transverse_distance(0, z_i) * (1 + z_i)).in_cgs()

    def lookback_time(self, z_i, z_f):
        return (trapzint(self.age_integrand, z_i, z_f) / \
                self.hubble_constant).in_cgs()
    
    def hubble_time(self, z, z_inf=1e6):
        return (trapzint(self.age_integrand, z, z_inf) /
                self.hubble_constant).in_cgs()

    def angular_scale_1arcsec_kpc(self, z_i, z_f):
        return (self.angular_diameter_distance(z_i, z_f) /
                648. * np.pi).in_cgs()

    def critical_density(self, z):
        return (3.0 / 8.0 / np.pi * 
                self.hubble_constant**2 / G *
                ((1 + z)**3.0 * self.omega_matter + 
                 self.omega_lambda)).in_cgs()

    def age_integrand(self, z):
        return (1 / (z + 1) / self.expansion_factor(z))

    def expansion_factor(self, z):
        return np.sqrt(self.omega_matter * ((1 + z)**3.0) + 
                    self.omega_curvature * np.sqrt(1 + z) + 
                    self.omega_lambda)

    def inverse_expansion_factor(self, z):
        return 1 / self.expansion_factor(z)

    def path_length_function(self, z):
        return ((1 + z)**2) * self.inverse_expansion_factor(z)

    def path_length(self, z_i, z_f):
        return trapzint(self.path_length_function, z_i, z_f)

    def z_from_t(self, my_time):
        """
        Compute the redshift from time after the big bang.  This is based on
        Enzo's CosmologyComputeExpansionFactor.C, but altered to use physical
        units.
        """

        omega_curvature = 1.0 - self.omega_matter - self.omega_lambda

        OMEGA_TOLERANCE = 1e-5
        ETA_TOLERANCE = 1.0e-10

        # Convert the time to Time * H0.

        if not isinstance(my_time, YTArray):
            my_time = YTQuantity(my_time, "s")

        t0 = (my_time.in_units("s") *
              self.hubble_constant.in_units("1/s")).to_ndarray()
 
        # 1) For a flat universe with omega_matter = 1, it's easy.
 
        if ((np.fabs(self.omega_matter-1) < OMEGA_TOLERANCE) and
            (self.omega_lambda < OMEGA_TOLERANCE)):
            a = np.power(my_time/self.InitialTime, 2.0/3.0)
 
        # 2) For omega_matter < 1 and omega_lambda == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
        #    Actually, this is a little tricky since we must solve an equation
        #    of the form eta - np.sinh(eta) + x = 0..
 
        if ((self.omega_matter < 1) and 
            (self.omega_lambda < OMEGA_TOLERANCE)):
            x = 2*t0*np.power(1.0 - self.omega_matter, 1.5) / \
                self.omega_matter;
 
            # Compute eta in a three step process, first from a third-order
            # Taylor expansion of the formula above, then use that in a fifth-order
            # approximation.  Then finally, iterate on the formula itself, solving for
            # eta.  This works well because parts 1 & 2 are an excellent approximation
            # when x is small and part 3 converges quickly when x is large. 
 
            eta = np.power(6*x, 1.0/3.0)                # part 1
            eta = np.power(120*x/(20+eta*eta), 1.0/3.0) # part 2
            for i in range(40):                      # part 3
                eta_old = eta
                eta = np.arcsinh(eta + x)
                if (np.fabs(eta-eta_old) < ETA_TOLERANCE): 
                    break
                if (i == 39):
                    print "No convergence after %d iterations." % i
 
            # Now use eta to compute the expansion factor (eq. 13-10, part 2).
 
            a = self.omega_matter/(2.0*(1.0 - self.omega_matter))*\
                (np.cosh(eta) - 1.0)

        # 3) For omega_matter > 1 and omega_lambda == 0, use sin/cos.
        #    Easy, but skip it for now.
 
        if ((self.omega_matter > 1) and 
            (self.omega_lambda < OMEGA_TOLERANCE)):
            print "Never implemented in Enzo, not implemented here."
            return 0
 
        # 4) For flat universe, with non-zero omega_lambda, see eq. 13-20.
 
        if ((np.fabs(omega_curvature) < OMEGA_TOLERANCE) and
            (self.omega_lambda > OMEGA_TOLERANCE)):
            a = np.power(self.omega_matter / 
                         (1 - self.omega_matter), 1.0/3.0) * \
                np.power(np.sinh(1.5 * np.sqrt(1.0 - self.omega_matter)*\
                                     t0), 2.0/3.0)


        redshift = (1.0/a) - 1.0

        return redshift

    def t_from_z(self, z):
        """
        Compute the time from redshift.  This is based on Enzo's
        CosmologyComputeTimeFromRedshift.C, but altered to use physical units.
        """
        omega_curvature = 1.0 - self.omega_matter - self.omega_lambda
 
        # 1) For a flat universe with omega_matter = 1, things are easy.
 
        if ((self.omega_matter == 1.0) and (self.omega_lambda == 0.0)):
            t0 = 2.0/3.0/np.power(1+z, 1.5)
 
        # 2) For omega_matter < 1 and omega_lambda == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
 
        if ((self.omega_matter < 1) and (self.omega_lambda == 0)):
            eta = np.arccosh(1 + 
                             2*(1-self.omega_matter)/self.omega_matter/(1+z))
            t0 = self.omega_matter/ \
              (2*np.power(1.0-self.omega_matter, 1.5))*\
              (np.sinh(eta) - eta)
 
        # 3) For omega_matter > 1 and omega_lambda == 0, use sin/cos.
 
        if ((self.omega_matter > 1) and (self.omega_lambda == 0)):
            eta = np.arccos(1 - 2*(1-self.omega_matter)/self.omega_matter/(1+z))
            t0 = self.omega_matter/(2*np.power(1.0-self.omega_matter, 1.5))*\
                (eta - np.sin(eta))
 
        # 4) For flat universe, with non-zero omega_lambda, see eq. 13-20.
 
        if ((np.fabs(omega_curvature) < 1.0e-3) and (self.omega_lambda != 0)):
            t0 = 2.0/3.0/np.sqrt(1-self.omega_matter)*\
                np.arcsinh(np.sqrt((1-self.omega_matter)/self.omega_matter)/ \
                               np.power(1+z, 1.5))
  
        # Now convert from Time * H0 to time.
  
        my_time = t0 / self.hubble_constant
    
        return my_time.in_cgs()

def trapzint(f, a, b, bins=10000):
    zbins = np.logspace(np.log10(a + 1), np.log10(b + 1), bins) - 1
    return np.trapz(f(zbins[:-1]), x=zbins[:-1], dx=np.diff(zbins))

class EnzoCosmology(object):
    def __init__(self, HubbleConstantNow = 71.0,
                 omega_matter = 0.27,
                 omega_lambda = 0.73,
                 omega_curvature = 0.0,
                 InitialRedshift = 99.0):
        self.HubbleConstantNow = HubbleConstantNow
        self.omega_matter = omega_matter
        self.omega_lambda = omega_lambda
        self.omega_curvature = omega_curvature
        self.InitialRedshift = InitialRedshift
        self.InitialTime = self.ComputeTimeFromRedshift(self.InitialRedshift)
        self.TimeUnits = self.ComputeTimeUnits()

    def ComputeTimeUnits(self):
        """
        Taken from CosmologyGetUnits.C in Enzo.
        """
        # Changed 2.52e17 to 2.52e19 because H_0 is in km/s/Mpc, 
        # instead of 100 km/s/Mpc.
        # TODO: Move me to physical_units
        return 2.52e19 / np.sqrt(self.omega_matter) / \
            self.HubbleConstantNow / np.power(1 + self.InitialRedshift,1.5)

    def ComputeRedshiftFromTime(self,time):
        """
        Compute the redshift from time after the big bang.  This is based on
        Enzo's CosmologyComputeExpansionFactor.C, but altered to use physical
        units.
        """

        omega_curvature = 1.0 - self.omega_matter - self.omega_lambda

        OMEGA_TOLERANCE = 1e-5
        ETA_TOLERANCE = 1.0e-10

        # Convert the time to Time * H0.
 
        TimeHubble0 = time * self.HubbleConstantNow / kmPerMpc
 
        # 1) For a flat universe with omega_matter = 1, it's easy.
 
        if ((np.fabs(self.omega_matter-1) < OMEGA_TOLERANCE) and
            (self.omega_lambda < OMEGA_TOLERANCE)):
            a = np.power(time/self.InitialTime,2.0/3.0)
 
        # 2) For omega_matter < 1 and omega_lambda == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
        #    Actually, this is a little tricky since we must solve an equation
        #    of the form eta - np.sinh(eta) + x = 0..
 
        if ((self.omega_matter < 1) and 
            (self.omega_lambda < OMEGA_TOLERANCE)):
            x = 2*TimeHubble0*np.power(1.0 - self.omega_matter, 1.5) / \
                self.omega_matter;
 
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
 
            a = self.omega_matter/(2.0*(1.0 - self.omega_matter))*\
                (np.cosh(eta) - 1.0)

        # 3) For omega_matter > 1 and omega_lambda == 0, use sin/cos.
        #    Easy, but skip it for now.
 
        if ((self.omega_matter > 1) and 
            (self.omega_lambda < OMEGA_TOLERANCE)):
            print "Never implemented in Enzo, not implemented here."
            return 0
 
        # 4) For flat universe, with non-zero omega_lambda, see eq. 13-20.
 
        if ((np.fabs(omega_curvature) < OMEGA_TOLERANCE) and
            (self.omega_lambda > OMEGA_TOLERANCE)):
            a = np.power(self.omega_matter / (1 - self.omega_matter),1.0/3.0) * \
                np.power(np.sinh(1.5 * np.sqrt(1.0 - self.omega_matter)*\
                                     TimeHubble0),2.0/3.0)


        redshift = (1.0/a) - 1.0

        return redshift

    def ComputeTimeFromRedshift(self,z):
        """
        Compute the time from redshift.  This is based on Enzo's
        CosmologyComputeTimeFromRedshift.C, but altered to use physical units.
        """
        omega_curvature = 1.0 - self.omega_matter - self.omega_lambda
 
        # 1) For a flat universe with omega_matter = 1, things are easy.
 
        if ((self.omega_matter == 1.0) and (self.omega_lambda == 0.0)):
            TimeHubble0 = 2.0/3.0/np.power(1+z,1.5)
 
        # 2) For omega_matter < 1 and omega_lambda == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
 
        if ((self.omega_matter < 1) and (self.omega_lambda == 0)):
            eta = np.arccosh(1 + 2*(1-self.omega_matter)/self.omega_matter/(1+z))
            TimeHubble0 = self.omega_matter/(2*np.power(1.0-self.omega_matter, 1.5))*\
                (np.sinh(eta) - eta)
 
        # 3) For omega_matter > 1 and omega_lambda == 0, use sin/cos.
 
        if ((self.omega_matter > 1) and (self.omega_lambda == 0)):
            eta = np.arccos(1 - 2*(1-self.omega_matter)/self.omega_matter/(1+z))
            TimeHubble0 = self.omega_matter/(2*np.power(1.0-self.omega_matter, 1.5))*\
                (eta - np.sin(eta))
 
        # 4) For flat universe, with non-zero omega_lambda, see eq. 13-20.
 
        if ((np.fabs(omega_curvature) < 1.0e-3) and (self.omega_lambda != 0)):
            TimeHubble0 = 2.0/3.0/np.sqrt(1-self.omega_matter)*\
                np.arcsinh(np.sqrt((1-self.omega_matter)/self.omega_matter)/ \
                               np.power(1+z,1.5))
  
        # Now convert from Time * H0 to time.
  
        time = TimeHubble0 / (self.HubbleConstantNow/kmPerMpc)
    
        return time
