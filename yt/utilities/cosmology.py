import functools

import numpy as np

from yt.units import dimensions
from yt.units.unit_object import Unit  # type: ignore
from yt.units.unit_registry import UnitRegistry  # type: ignore
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.physical_constants import (
    gravitational_constant_cgs as G,
    speed_of_light_cgs,
)


class Cosmology:
    r"""
    Create a cosmology calculator to compute cosmological distances and times.

    For an explanation of the various cosmological measures, see, for example
    Hogg (1999, https://arxiv.org/abs/astro-ph/9905116).

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
    omega_radiation : the fraction of the energy density of the Universe in
        relativistic matter at redshift zero.
    omega_curvature : the fraction of the energy density of the Universe in
        curvature.
        Default: 0.0.
    unit_system : :class:`yt.units.unit_systems.UnitSystem`, optional
        The units system to use when making calculations. If not specified,
        cgs units are assumed.
    use_dark_factor: Bool, optional
        The flag to either use the cosmological constant (False, default)
        or to use the parameterization of w(a) as given in Linder 2002. This,
        along with w_0 and w_a, only matters in the function expansion_factor.
    w_0 : float, optional
        The Linder 2002 parameterization of w(a) is: w(a) = w_0 + w_a(1 - a).
        w_0 is w(a = 1). Only matters if use_dark_factor = True. Default is None.
        Cosmological constant case corresponds to w_0 = -1.
    w_a : float, optional
        See w_0. w_a is the derivative of w(a) evaluated at a = 1. Cosmological
        constant case corresponds to w_a = 0. Default is None.

    Examples
    --------

    >>> from yt.utilities.cosmology import Cosmology
    >>> co = Cosmology()
    >>> print(co.t_from_z(0.0).in_units("Gyr"))

    """

    def __init__(
        self,
        hubble_constant=0.71,
        omega_matter=0.27,
        omega_lambda=0.73,
        omega_radiation=0.0,
        omega_curvature=0.0,
        unit_registry=None,
        unit_system="cgs",
        use_dark_factor=False,
        w_0=-1.0,
        w_a=0.0,
    ):
        self.omega_matter = float(omega_matter)
        self.omega_radiation = float(omega_radiation)
        self.omega_lambda = float(omega_lambda)
        self.omega_curvature = float(omega_curvature)
        hubble_constant = float(hubble_constant)
        if unit_registry is None:
            unit_registry = UnitRegistry(unit_system=unit_system)
            unit_registry.add("h", hubble_constant, dimensions.dimensionless, r"h")
            for my_unit in ["m", "pc", "AU", "au"]:
                new_unit = f"{my_unit}cm"
                my_u = Unit(my_unit, registry=unit_registry)
                # technically not true, but distances here are actually comoving
                unit_registry.add(
                    new_unit,
                    my_u.base_value,
                    dimensions.length,
                    "\\rm{%s}/(1+z)" % my_unit,
                    prefixable=True,
                )
        self.unit_registry = unit_registry
        self.hubble_constant = self.quan(hubble_constant, "100*km/s/Mpc")
        self.unit_system = unit_system

        # For non-standard dark energy. If false, use default cosmological constant
        # This only affects the expansion_factor function.
        self.use_dark_factor = use_dark_factor
        self.w_0 = w_0
        self.w_a = w_a

    def hubble_distance(self):
        r"""
        The distance corresponding to c / h, where c is the speed of light
        and h is the Hubble parameter in units of 1 / time.
        """
        return self.quan(speed_of_light_cgs / self.hubble_constant).in_base(
            self.unit_system
        )

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
        >>> print(co.comoving_radial_distance(0.0, 1.0).in_units("Mpccm"))

        """
        return (
            self.hubble_distance() * trapzint(self.inverse_expansion_factor, z_i, z_f)
        ).in_base(self.unit_system)

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
        >>> print(co.comoving_transverse_distance(0.0, 1.0).in_units("Mpccm"))

        """
        if self.omega_curvature > 0:
            return (
                self.hubble_distance()
                / np.sqrt(self.omega_curvature)
                * np.sinh(
                    np.sqrt(self.omega_curvature)
                    * self.comoving_radial_distance(z_i, z_f)
                    / self.hubble_distance()
                )
            ).in_base(self.unit_system)
        elif self.omega_curvature < 0:
            return (
                self.hubble_distance()
                / np.sqrt(np.fabs(self.omega_curvature))
                * np.sin(
                    np.sqrt(np.fabs(self.omega_curvature))
                    * self.comoving_radial_distance(z_i, z_f)
                    / self.hubble_distance()
                )
            ).in_base(self.unit_system)
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
        >>> print(co.comoving_volume(0.0, 1.0).in_units("Gpccm**3"))

        """
        if self.omega_curvature > 1e-10:
            return (
                2
                * np.pi
                * np.power(self.hubble_distance(), 3)
                / self.omega_curvature
                * (
                    self.comoving_transverse_distance(z_i, z_f)
                    / self.hubble_distance()
                    * np.sqrt(
                        1
                        + self.omega_curvature
                        * np.sqrt(
                            self.comoving_transverse_distance(z_i, z_f)
                            / self.hubble_distance()
                        )
                    )
                    - np.sinh(
                        np.fabs(self.omega_curvature)
                        * self.comoving_transverse_distance(z_i, z_f)
                        / self.hubble_distance()
                    )
                    / np.sqrt(self.omega_curvature)
                )
            ).in_base(self.unit_system)
        elif self.omega_curvature < -1e-10:
            return (
                2
                * np.pi
                * np.power(self.hubble_distance(), 3)
                / np.fabs(self.omega_curvature)
                * (
                    self.comoving_transverse_distance(z_i, z_f)
                    / self.hubble_distance()
                    * np.sqrt(
                        1
                        + self.omega_curvature
                        * np.sqrt(
                            self.comoving_transverse_distance(z_i, z_f)
                            / self.hubble_distance()
                        )
                    )
                    - np.arcsin(
                        np.fabs(self.omega_curvature)
                        * self.comoving_transverse_distance(z_i, z_f)
                        / self.hubble_distance()
                    )
                    / np.sqrt(np.fabs(self.omega_curvature))
                )
            ).in_base(self.unit_system)
        else:
            return (
                4 * np.pi * np.power(self.comoving_transverse_distance(z_i, z_f), 3) / 3
            ).in_base(self.unit_system)

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
        >>> print(co.angular_diameter_distance(0.0, 1.0).in_units("Mpc"))

        """

        return (
            self.comoving_transverse_distance(0, z_f) / (1 + z_f)
            - self.comoving_transverse_distance(0, z_i) / (1 + z_i)
        ).in_base(self.unit_system)

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
        >>> print(co.angular_scale(0.0, 1.0).in_units("kpc / arcsec"))

        """

        scale = self.angular_diameter_distance(z_i, z_f) / self.quan(1, "radian")
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
        >>> print(co.luminosity_distance(0.0, 1.0).in_units("Mpc"))

        """

        return (
            self.comoving_transverse_distance(0, z_f) * (1 + z_f)
            - self.comoving_transverse_distance(0, z_i) * (1 + z_i)
        ).in_base(self.unit_system)

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
        >>> print(co.lookback_time(0.0, 1.0).in_units("Gyr"))

        """
        return (trapzint(self.age_integrand, z_i, z_f) / self.hubble_constant).in_base(
            self.unit_system
        )

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
        >>> print(co.critical_density(0.0).in_units("g/cm**3"))
        >>> print(co.critical_density(0).in_units("Msun/Mpc**3"))

        """
        return (3.0 * self.hubble_parameter(z) ** 2 / 8.0 / np.pi / G).in_base(
            self.unit_system
        )

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
        return 1.0 / (z + 1) / self.expansion_factor(z)

    def expansion_factor(self, z):
        r"""
        The ratio between the Hubble parameter at a given redshift and
        redshift zero.

        This is also the primary function integrated to calculate the
        cosmological distances.

        """

        # Use non-standard dark energy
        if self.use_dark_factor:
            dark_factor = self.get_dark_factor(z)

        # Use default cosmological constant
        else:
            dark_factor = 1.0

        zp1 = 1 + z
        return np.sqrt(
            self.omega_matter * zp1**3
            + self.omega_curvature * zp1**2
            + self.omega_radiation * zp1**4
            + self.omega_lambda * dark_factor
        )

    def inverse_expansion_factor(self, z):
        return 1.0 / self.expansion_factor(z)

    def path_length_function(self, z):
        return ((1 + z) ** 2) * self.inverse_expansion_factor(z)

    def path_length(self, z_i, z_f):
        return trapzint(self.path_length_function, z_i, z_f)

    def t_from_a(self, a):
        """
        Compute the age of the Universe for a given scale factor.

        Parameters
        ----------
        a : float
            Scale factor.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.t_from_a(1.0).in_units("Gyr"))

        """

        # Interpolate from a table of log(a) vs. log(t)
        la = np.log10(a)
        la_i = min(-6, np.asarray(la).min() - 3)
        la_f = np.asarray(la).max()
        bins_per_dex = 1000
        n_bins = int((la_f - la_i) * bins_per_dex + 1)
        la_bins = np.linspace(la_i, la_f, n_bins)
        z_bins = 1.0 / np.power(10, la_bins) - 1

        # Integrate in redshift.
        lt = trapezoid_cumulative_integral(self.age_integrand, z_bins)

        # Add a minus sign because we've switched the integration limits.
        table = InterpTable(la_bins[1:], np.log10(-lt))
        t = np.power(10, table(la))

        return (t / self.hubble_constant).in_base(self.unit_system)

    def t_from_z(self, z):
        """
        Compute the age of the Universe for a given redshift.

        Parameters
        ----------
        z : float
            Redshift.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.t_from_z(0.0).in_units("Gyr"))

        """

        return self.t_from_a(1.0 / (1.0 + z))

    def a_from_t(self, t):
        """
        Compute the scale factor for a given age of the Universe.

        Parameters
        ----------
        t : YTQuantity or float
            Time since the Big Bang.  If a float is given, units are
            assumed to be seconds.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.a_from_t(4.0e17))

        """

        if not isinstance(t, YTArray):
            t = self.arr(t, "s")
        lt = np.log10((t * self.hubble_constant).to(""))

        # Interpolate from a table of log(a) vs. log(t)
        # Make initial guess for bounds and widen if necessary.
        la_i = -6
        la_f = 6
        bins_per_dex = 1000
        iter = 0

        while True:
            good = True
            n_bins = int((la_f - la_i) * bins_per_dex + 1)
            la_bins = np.linspace(la_i, la_f, n_bins)
            z_bins = 1.0 / np.power(10, la_bins) - 1

            # Integrate in redshift.
            lt_bins = trapezoid_cumulative_integral(self.age_integrand, z_bins)

            # Add a minus sign because we've switched the integration limits.
            table = InterpTable(np.log10(-lt_bins), la_bins[1:])
            la = table(lt)
            # We want to have the la_bins lower bound be decently
            # below the minimum calculated la values.

            laa = np.asarray(la)
            if laa.min() < la_i + 2:
                la_i -= 3
                good = False
            if laa.max() > la_f:
                la_f = laa.max() + 1
                good = False
            if good:
                break
            iter += 1
            if iter > 10:
                raise RuntimeError("a_from_t calculation did not converge!")

        a = np.power(10, table(lt))
        return a

    def z_from_t(self, t):
        """
        Compute the redshift for a given age of the Universe.

        Parameters
        ----------
        t : YTQuantity or float
            Time since the Big Bang.  If a float is given, units are
            assumed to be seconds.

        Examples
        --------

        >>> from yt.utilities.cosmology import Cosmology
        >>> co = Cosmology()
        >>> print(co.z_from_t(4.0e17))

        """

        a = self.a_from_t(t)
        return 1.0 / a - 1.0

    def get_dark_factor(self, z):
        """
        This function computes the additional term that enters the expansion factor
        when using non-standard dark energy. See Dolag et al 2004 eq. 7 for ref (but
        note that there's a typo in his eq. There should be no negative sign).

        At the moment, this only works using the parameterization given in Linder 2002
        eq. 7: w(a) = w0 + wa(1 - a) = w0 + wa * z / (1+z). This gives rise to an
        analytic expression.
        It is also only functional for Gadget simulations, at the moment.

        Parameters
        ----------
        z:  float
            Redshift
        """

        # Get value of scale factor a corresponding to redshift z
        scale_factor = 1.0 / (1.0 + z)

        # Evaluate exponential using Linder02 parameterization
        dark_factor = np.power(
            scale_factor, -3.0 * (1.0 + self.w_0 + self.w_a)
        ) * np.exp(-3.0 * self.w_a * (1.0 - scale_factor))

        return dark_factor

    _arr = None

    @property
    def arr(self):
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(YTArray, registry=self.unit_registry)
        return self._arr

    _quan = None

    @property
    def quan(self):
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(YTQuantity, registry=self.unit_registry)
        return self._quan


def trapzint(f, a, b, bins=10000):
    zbins = np.logspace(np.log10(a + 1), np.log10(b + 1), bins) - 1
    return np.trapz(f(zbins[:-1]), x=zbins[:-1], dx=np.diff(zbins))


def trapezoid_cumulative_integral(f, x):
    """
    Perform cumulative integration using the trapezoid rule.
    """

    fy = f(x)
    return (0.5 * (fy[:-1] + fy[1:]) * np.diff(x)).cumsum()


class InterpTable:
    """
    Generate a function to linearly interpolate from provided arrays.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __call__(self, val):
        i = np.clip(np.digitize(val, self.x) - 1, 0, self.x.size - 2)
        slope = (self.y[i + 1] - self.y[i]) / (self.x[i + 1] - self.x[i])
        return slope * (val - self.x[i]) + self.y[i]
