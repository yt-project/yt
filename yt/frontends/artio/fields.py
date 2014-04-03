"""
ARTIO-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.units.yt_array import \
    YTArray

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs, \
    boltzmann_constant_cgs, \
    amu_cgs

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
# NOTE: ARTIO uses momentum density.
mom_units = "code_mass / (code_length**2 * code_time)"
en_units = "code_mass*code_velocity**2/code_length**3"

class ARTIOFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("HVAR_GAS_DENSITY", (rho_units, ["density"], None)),
        ("HVAR_GAS_ENERGY", (en_units, ["total_energy"], None)),
        ("HVAR_INTERNAL_ENERGY", (en_units, ["thermal_energy"], None)),
        ("HVAR_PRESSURE", ("", ["pressure"], None)), # Unused
        ("HVAR_MOMENTUM_X", (mom_units, ["momentum_x"], None)),
        ("HVAR_MOMENTUM_Y", (mom_units, ["momentum_y"], None)),
        ("HVAR_MOMENTUM_Z", (mom_units, ["momentum_z"], None)),
        ("HVAR_GAMMA", ("", ["gamma"], None)),
        ("HVAR_METAL_DENSITY_Ia", (rho_units, ["metal_ia_density"], None)),
        ("HVAR_METAL_DENSITY_II", (rho_units, ["metal_ii_density"], None)),
        ("VAR_POTENTIAL", ("", ["potential"], None)),
        ("VAR_POTENTIAL_HYDRO", ("", ["gas_potential"], None)),
    )

    known_particle_fields = (
        ("POSITION_X", ("code_length", ["particle_position_x"], None)),
        ("POSITION_Y", ("code_length", ["particle_position_y"], None)),
        ("POSITION_Z", ("code_length", ["particle_position_z"], None)),
        ("VELOCITY_X", (vel_units, ["particle_velocity_x"], None)),
        ("VELOCITY_Y", (vel_units, ["particle_velocity_y"], None)),
        ("VELOCITY_Z", (vel_units, ["particle_velocity_z"], None)),
        ("MASS", ("code_mass", ["particle_mass"], None)),
        ("PID", ("", ["particle_index"], None)),
        ("SPECIES", ("", ["particle_type"], None)),
        ("BIRTH_TIME", ("code_time", ["creation_time"], None)),
        ("INITIAL_MASS", ("code_mass", ["initial_mass"], None)),
        ("METALLICITY_SNIa", ("", ["metallicity_snia"], None)),
        ("METALLICITY_SNII", ("", ["metallicity_snii"], None)),
    )

    def setup_fluid_fields(self):
        def _get_vel(axis):
            def velocity(field, data):
                return data["momentum_%s" % axis]/data["density"]
            return velocity
        for ax in 'xyz':
            self.add_field(("gas", "velocity_%s" % ax),
                           function = _get_vel(ax),
                           units = "cm/s")

        def _temperature(field, data):
            tr = data["thermal_energy"]/data["density"]
            # We want this to match *exactly* what ARTIO would compute
            # internally.  We therefore use the exact values that are internal
            # to ARTIO, rather than yt's own internal constants.
            mH  = 1.007825*amu_cgs
            mHe = 4.002602*amu_cgs
            Yp    = 0.24
            XH    = 1.0 - Yp
            XHe   = 0.25*Yp
            mb = XH*mH + XHe*mHe
            wmu   = 4.0/(8.0-5.0*Yp)
            # Note that we have gamma = 5.0/3.0 here
            tr *= (data["gamma"] - 1.0)
            tr *= wmu
            tr *= mb/boltzmann_constant_cgs
            return tr
        # TODO: The conversion factor here needs to be addressed, as previously
        # it was set as:
        # unit_T = unit_v**2.0*mb / constants.k
        self.add_field(("gas", "temperature"), function = _temperature,
                       units = "K")

        def _metal_density(field, data):
            tr = data["metal_ia_density"]
            tr += data["metal_ii_density"]
            return tr
        self.add_field(("gas","metal_density"),
                       function=_metal_density,
                       units="g/cm**3",
                       take_log=True)

    def setup_particle_fields(self, ptype):

        def _particle_age(field, data):
            return b2t(data[ptype,"creation_time"])
        self.add_field((ptype, "particle_age"), function=_particle_age, units="s",
                  particle_type=True)

        super(ARTIOFieldInfo, self).setup_particle_fields(ptype)

#stolen from frontends/art/
#All of these functions are to convert from hydro time var to
#proper time
sqrt = np.sqrt
sign = np.sign


def find_root(f, a, b, tol=1e-6):
    c = (a+b)/2.0
    last = -np.inf
    assert(sign(f(a)) != sign(f(b)))
    while np.abs(f(c)-last) > tol:
        last = f(c)
        if sign(last) == sign(f(b)):
            b = c
        else:
            a = c
        c = (a+b)/2.0
    return c


def quad(fintegrand, xmin, xmax, n=1e4):
    spacings = np.logspace(np.log10(xmin), np.log10(xmax), n)
    integrand_arr = fintegrand(spacings)
    val = np.trapz(integrand_arr, dx=np.diff(spacings))
    return val


def a2b(at, Om0=0.27, Oml0=0.73, h=0.700):
    def f_a2b(x):
        val = 0.5*sqrt(Om0) / x**3.0
        val /= sqrt(Om0/x**3.0 + Oml0 + (1.0-Om0-Oml0)/x**2.0)
        return val
    #val, err = si.quad(f_a2b,1,at)
    val = quad(f_a2b, 1, at)
    return val


def b2a(bt, **kwargs):
    #converts code time into expansion factor
    #if Om0 ==1and OmL == 0 then b2a is (1 / (1-td))**2
    #if bt < -190.0 or bt > -.10:  raise 'bt outside of range'
    f_b2a = lambda at: a2b(at, **kwargs)-bt
    return find_root(f_b2a, 1e-4, 1.1)
    #return so.brenth(f_b2a,1e-4,1.1)
    #return brent.brent(f_b2a)


def a2t(at, Om0=0.27, Oml0=0.73, h=0.700):
    integrand = lambda x: 1./(x*sqrt(Oml0+Om0*x**-3.0))
    #current_time,err = si.quad(integrand,0.0,at,epsabs=1e-6,epsrel=1e-6)
    current_time = quad(integrand, 1e-4, at)
    #spacings = np.logspace(-5,np.log10(at),1e5)
    #integrand_arr = integrand(spacings)
    #current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
    current_time *= 9.779/h
    return current_time


def b2t(tb, n=1e2, logger=None, **kwargs):
    tb = np.array(tb)
    if len(np.atleast_1d(tb)) == 1: 
        return a2t(b2a(tb))
    if tb.shape == ():
        return None 
    if len(tb) < n:
        n = len(tb)
    age_min = a2t(b2a(tb.max(), **kwargs), **kwargs)
    age_max = a2t(b2a(tb.min(), **kwargs), **kwargs)
    tbs = -1.*np.logspace(np.log10(-tb.min()),
                          np.log10(-tb.max()), n)
    ages = []
    for i, tbi in enumerate(tbs):
        ages += a2t(b2a(tbi)),
        if logger:
            logger(i)
    ages = np.array(ages)
    fb2t = np.interp(tb, tbs, ages)
    #fb2t = interp1d(tbs,ages)
    return fb2t*1e9*31556926


def spread_ages(ages, logger=None, spread=.0e7*365*24*3600):
    #stars are formed in lumps; spread out the ages linearly
    da = np.diff(ages)
    assert np.all(da <= 0)
    #ages should always be decreasing, and ordered so
    agesd = np.zeros(ages.shape)
    idx, = np.where(da < 0)
    idx += 1
    #mark the right edges
    #spread this age evenly out to the next age
    lidx = 0
    lage = 0
    for i in idx:
        n = i-lidx
        #n stars affected
        rage = ages[i]
        lage = max(rage-spread, 0.0)
        agesd[lidx:i] = np.linspace(lage, rage, n)
        lidx = i
        #lage=rage
        if logger:
            logger(i)
    #we didn't get the last iter
    i = ages.shape[0]-1
    n = i-lidx
    #n stars affected
    rage = ages[i]
    lage = max(rage-spread, 0.0)
    agesd[lidx:i] = np.linspace(lage, rage, n)
    return agesd
