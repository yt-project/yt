"""
Cosmology calculator.
Cosmology calculator based originally on http://www.kempner.net/cosmic.php 
and featuring time and redshift conversion functions from Enzo.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013-2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import functools
import numpy as np

from yt.units import dimensions
from yt.units.unit_registry import \
     UnitRegistry
from yt.units.yt_array import \
     YTArray, \
     YTQuantity

from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G, \
    speed_of_light_cgs

class Cosmology(object):
    r"""
    Create a cosmology calculator to compute cosmological distances and times.

    For an explanation of the various cosmological measures, see, for example 
    Hogg (1999, http://xxx.lanl.gov/abs/astro-ph/9905116).

    WARNING: Cosmological distance calculations return values that are either
    in the comoving or proper frame, depending on the specific quantity.  For
    simplicity, the proper and comoving frames are set equal to each other
    within the cosmology calculator.  This means that for some distance value,
    x, x.to("Mpc") and x.to("Mpccm") will be the same.  The user should take
    care to understand which reference frame is correct for the given calculation.

    Parameters
    ----------
    hubble_constant : float
        The Hubble parameter at redshift zero in units of 100 km/s/Mpc.
        Default: 0.71.
    omega_matter : the fraction of the energy density of the Universe in 
        matter at redshift zero.
        Default: 0.27.
    omega_lambda : the fraction of the energy density of the Universe in 
        a cosmological constant.
        Default: 0.73.
    omega_curvature : the fraction of the energy density of the Universe in 
        curvature.
        Default: 0.0.
    unit_system : :class:`yt.units.unit_systems.UnitSystem`, optional
        The units system to use when making calculations. If not specified,
        cgs units are assumed.

    Examples
    --------

    >>> from yt.utilities.cosmology import Cosmology
    >>> co = Cosmology()
    >>> print(co.hubble_time(0.0).in_units("Gyr"))

    """
    def __init__(self, hubble_constant = 0.71,
                 omega_matter = 0.27,
                 omega_lambda = 0.73,
                 omega_curvature = 0.0,
                 unit_registry = None,
                 unit_system = "cgs"):
        self.omega_matter = float(omega_matter)
        self.omega_lambda = float(omega_lambda)
        self.omega_curvature = float(omega_curvature)
        if unit_registry is None:
            unit_registry = UnitRegistry()
            unit_registry.modify("h", hubble_constant)
            for my_unit in ["m", "pc", "AU", "au"]:
                new_unit = "%scm" % my_unit
                # technically not true, but distances here are actually comoving
                unit_registry.add(new_unit, unit_registry.lut[my_unit][0],
                                  dimensions.length, "\\rm{%s}/(1+z)" % my_unit)
        self.unit_registry = unit_registry
        self.hubble_constant = self.quan(hubble_constant, "100*km/s/Mpc")
        self.unit_system = unit_system

    def hubble_distance(self):
        r"""
        The distance corresponding to c / h, where c is the speed of light 
        and h is the Hubble parameter in units of 1 / time.
        """
        return self.quan((speed_of_light_cgs / self.hubble_constant)).in_base(self.unit_system)

    def comoving_radial_distance(self, z_i, z_f):
        r"""
        The comoving distance along the line of sight to on object at redshift, 
        z_f, viewed at a redshift, z_i.

        Parameters
        ----------
        z_i : float
            The redshift of the observer.
        z_f : float
            The redshift of the observed object.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.comoving_radial_distance(0., 1.).in_units("Mpccm"))
        
        """
        return (self.hubble_distance() *
                trapzint(self.inverse_expansion_factor, z_i, z_f)).in_base(self.unit_system)

    def comoving_transverse_distance(self, z_i, z_f):
        r"""
        When multiplied by some angle, the distance between two objects 
        observed at redshift, z_f, with an angular separation given by that 
        angle, viewed by an observer at redshift, z_i (Hogg 1999).

        Parameters
        ----------
        z_i : float
            The redshift of the observer.
        z_f : float
            The redshift of the observed object.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.comoving_transverse_distance(0., 1.).in_units("Mpccm"))
        
        """
        if (self.omega_curvature > 0):
            return (self.hubble_distance() / np.sqrt(self.omega_curvature) * 
                    np.sinh(np.sqrt(self.omega_curvature) * 
                            self.comoving_radial_distance(z_i, z_f) /
                            self.hubble_distance())).in_base(self.unit_system)
        elif (self.omega_curvature < 0):
            return (self.hubble_distance() /
                    np.sqrt(np.fabs(self.omega_curvature)) * 
                    np.sin(np.sqrt(np.fabs(self.omega_curvature)) * 
                           self.comoving_radial_distance(z_i, z_f) /
                           self.hubble_distance())).in_base(self.unit_system)
        else:
            return self.comoving_radial_distance(z_i, z_f)

    def comoving_volume(self, z_i, z_f):
        r"""
        "The comoving volume is the volume measure in which number densities 
        of non-evolving objects locked into Hubble flow are constant with 
        redshift." -- Hogg (1999)
        
        Parameters
        ----------
        z_i : float
            The lower redshift of the interval.
        z_f : float
            The higher redshift of the interval.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.comoving_volume(0., 1.).in_units("Gpccm**3"))

        """
        if (self.omega_curvature > 0):
             return (2 * np.pi * np.power(self.hubble_distance(), 3) /
                     self.omega_curvature * 
                     (self.comoving_transverse_distance(z_i, z_f) /
                      self.hubble_distance() * 
                      np.sqrt(1 + self.omega_curvature * 
                           np.sqrt(self.comoving_transverse_distance(z_i, z_f) /
                               self.hubble_distance())) - 
                      np.sinh(np.fabs(self.omega_curvature) * 
                            self.comoving_transverse_distance(z_i, z_f) /
                            self.hubble_distance()) /
                            np.sqrt(self.omega_curvature))).in_base(self.unit_system)
        elif (self.omega_curvature < 0):
             return (2 * np.pi * np.power(self.hubble_distance(), 3) /
                     np.fabs(self.omega_curvature) * 
                     (self.comoving_transverse_distance(z_i, z_f) /
                      self.hubble_distance() * 
                      np.sqrt(1 + self.omega_curvature * 
                           np.sqrt(self.comoving_transverse_distance(z_i, z_f) /
                               self.hubble_distance())) - 
                      np.arcsin(np.fabs(self.omega_curvature) * 
                           self.comoving_transverse_distance(z_i, z_f) /
                           self.hubble_distance()) /
                      np.sqrt(np.fabs(self.omega_curvature)))).in_base(self.unit_system)
        else:
             return (4 * np.pi *
                     np.power(self.comoving_transverse_distance(z_i, z_f), 3) /\
                     3).in_base(self.unit_system)

    def angular_diameter_distance(self, z_i, z_f):
        r"""
        Following Hogg (1999), the angular diameter distance is 'the ratio of 
        an object's physical transverse size to its angular size in radians.'

        Parameters
        ----------
        z_i : float
            The redshift of the observer.
        z_f : float
            The redshift of the observed object.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.angular_diameter_distance(0., 1.).in_units("Mpc"))
        
        """
        
        return (self.comoving_transverse_distance(0, z_f) / (1 + z_f) - 
                self.comoving_transverse_distance(0, z_i) / (1 + z_i)).in_base(self.unit_system)

    def angular_scale(self, z_i, z_f):
        r"""
        The proper transverse distance between two points at redshift z_f 
        observed at redshift z_i per unit of angular separation.

        Parameters
        ----------
        z_i : float
            The redshift of the observer.
        z_f : float
            The redshift of the observed object.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.angular_scale(0., 1.).in_units("kpc / arcsec"))
        
        """

        scale = self.angular_diameter_distance(z_i, z_f) / \
          self.quan(1, "radian")
        return scale.in_base(self.unit_system)
    
    def luminosity_distance(self, z_i, z_f):
        r"""
        The distance that would be inferred from the inverse-square law of 
        light and the measured flux and luminosity of the observed object.

        Parameters
        ----------
        z_i : float
            The redshift of the observer.
        z_f : float
            The redshift of the observed object.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.luminosity_distance(0., 1.).in_units("Mpc"))
        
        """

        return (self.comoving_transverse_distance(0, z_f) * (1 + z_f) - 
                self.comoving_transverse_distance(0, z_i) * (1 + z_i)).in_base(self.unit_system)

    def lookback_time(self, z_i, z_f):
        r"""
        The difference in the age of the Universe between the redshift interval 
        z_i to z_f.

        Parameters
        ----------
        z_i : float
            The lower redshift of the interval.
        z_f : float
            The higher redshift of the interval.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.lookback_time(0., 1.).in_units("Gyr"))

        """
        return (trapzint(self.age_integrand, z_i, z_f) / \
                self.hubble_constant).in_base(self.unit_system)
    
    def hubble_time(self, z, z_inf=1e6):
        r"""
        The age of the Universe at a given redshift.

        Parameters
        ----------
        z : float
            Redshift.
        z_inf : float
            The upper bound of the integral of the age integrand.
            Default: 1e6.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.hubble_time(0.).in_units("Gyr"))

        See Also
        --------

        t_from_z

        """
        return (trapzint(self.age_integrand, z, z_inf) /
                self.hubble_constant).in_base(self.unit_system)

    def critical_density(self, z):
        r"""
        The density required for closure of the Universe at a given 
        redshift in the proper frame.

        Parameters
        ----------
        z : float
            Redshift.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.critical_density(0.).in_units("g/cm**3"))
        >>> print(co.critical_density(0).in_units("Msun/Mpc**3"))
        
        """
        return (3.0 / 8.0 / np.pi * 
                self.hubble_constant**2 / G *
                ((1 + z)**3.0 * self.omega_matter + 
                 self.omega_lambda)).in_base(self.unit_system)

    def hubble_parameter(self, z):
        r"""
        The value of the Hubble parameter at a given redshift.

        Parameters
        ----------
        z: float
            Redshift.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.hubble_parameter(1.0).in_units("km/s/Mpc"))

        """
        return self.hubble_constant.in_base(self.unit_system) * self.expansion_factor(z)

    def age_integrand(self, z):
        return (1 / (z + 1) / self.expansion_factor(z))

    def expansion_factor(self, z):
        r"""
        The ratio between the Hubble parameter at a given redshift and 
        redshift zero.

        This is also the primary function integrated to calculate the 
        cosmological distances.
        
        """
        return np.sqrt(self.omega_matter * ((1 + z)**3.0) + 
                       self.omega_curvature * ((1 + z)**2.0) + 
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

        Parameters
        ----------
        my_time : float
            Age of the Universe in seconds.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.z_from_t(4.e17))

        """

        omega_curvature = 1.0 - self.omega_matter - self.omega_lambda

        OMEGA_TOLERANCE = 1e-5
        ETA_TOLERANCE = 1.0e-10

        # Convert the time to Time * H0.

        if not isinstance(my_time, YTArray):
            my_time = self.quan(my_time, "s")

        t0 = (my_time.in_units("s") *
              self.hubble_constant.in_units("1/s")).to_ndarray()
 
        # 1) For a flat universe with omega_matter = 1, it's easy.
 
        if ((np.fabs(self.omega_matter-1) < OMEGA_TOLERANCE) and
            (self.omega_lambda < OMEGA_TOLERANCE)):
            a = np.power(my_time/self.initial_time, 2.0/3.0)
 
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
                    print("No convergence after %d iterations." % i)
 
            # Now use eta to compute the expansion factor (eq. 13-10, part 2).
 
            a = self.omega_matter/(2.0*(1.0 - self.omega_matter))*\
                (np.cosh(eta) - 1.0)

        # 3) For omega_matter > 1 and omega_lambda == 0, use sin/cos.
        #    Easy, but skip it for now.
 
        if ((self.omega_matter > 1) and 
            (self.omega_lambda < OMEGA_TOLERANCE)):
            print("Never implemented in Enzo, not implemented here.")
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
        Compute the age of the Universe from redshift.  This is based on Enzo's
        CosmologyComputeTimeFromRedshift.C, but altered to use physical units.  
        Similar to hubble_time, but using an analytical function.

        Parameters
        ----------
        z : float
            Redshift.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.t_from_z(0.).in_units("Gyr"))

        See Also
        --------

        hubble_time
        
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
    
        return my_time.in_base(self.unit_system)

    _arr = None
    @property
    def arr(self):
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry = self.unit_registry)
        return self._arr
    
    _quan = None
    @property
    def quan(self):
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity,
                registry = self.unit_registry)
        return self._quan

def trapzint(f, a, b, bins=10000):
    zbins = np.logspace(np.log10(a + 1), np.log10(b + 1), bins) - 1
    return np.trapz(f(zbins[:-1]), x=zbins[:-1], dx=np.diff(zbins))
