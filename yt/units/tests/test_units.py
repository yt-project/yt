"""
Test symbolic unit handling.




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from numpy.testing import \
    assert_array_almost_equal_nulp, \
    assert_raises, assert_equal
from nose.tools import assert_true
import operator
from sympy import Symbol
from yt.testing import \
    fake_random_ds, assert_allclose_units, \
    assert_almost_equal
from yt.units.unit_registry import UnitRegistry
from yt.units import electrostatic_unit, elementary_charge
from yt.units.unit_object import default_unit_registry

# dimensions
from yt.units.dimensions import \
    mass, length, time, temperature, energy, magnetic_field, power, rate
# functions
from yt.units.unit_object import get_conversion_factor
# classes
from yt.units.unit_object import Unit, UnitParseError, InvalidUnitOperation
# objects
from yt.units.unit_lookup_table import \
    default_unit_symbol_lut, unit_prefixes, prefixable_units
import yt.units.unit_symbols as unit_symbols
# unit definitions
from yt.utilities.physical_ratios import \
    cm_per_pc, sec_per_year, cm_per_km, cm_per_mpc, \
    mass_sun_grams

def test_no_conflicting_symbols():
    """
    Check unit symbol definitions for conflicts.

    """
    full_set = set(default_unit_symbol_lut.keys())

    # go through all possible prefix combos
    for symbol in default_unit_symbol_lut.keys():
        if symbol in prefixable_units:
            keys = unit_prefixes.keys()
        else:
            keys = [symbol]
        for prefix in keys:
            new_symbol = "%s%s" % (prefix, symbol)

            # test if we have seen this symbol
            if new_symbol in full_set:
                print("Duplicate symbol: %s" % new_symbol)
                raise RuntimeError

            full_set.add(new_symbol)

def test_dimensionless():
    """
    Create dimensionless unit and check attributes.

    """
    u1 = Unit()

    assert_true(u1.is_dimensionless)
    assert_true(u1.expr == 1)
    assert_true(u1.base_value == 1)
    assert_true(u1.dimensions == 1)

    u2 = Unit("")

    assert_true(u2.is_dimensionless)
    assert_true(u2.expr == 1)
    assert_true(u2.base_value == 1)
    assert_true(u2.dimensions == 1)

    assert_equal(u1.latex_repr, '')
    assert_equal(u2.latex_repr, '')

#
# Start init tests
#

def test_create_from_string():
    """
    Create units with strings and check attributes.

    """

    u1 = Unit("g * cm**2 * s**-2")
    assert_true(u1.dimensions == energy)
    assert_true(u1.base_value == 1.0)

    # make sure order doesn't matter
    u2 = Unit("cm**2 * s**-2 * g")
    assert_true(u2.dimensions == energy)
    assert_true(u2.base_value == 1.0)

    # Test rationals
    u3 = Unit("g**0.5 * cm**-0.5 * s**-1")
    assert_true(u3.dimensions == magnetic_field)
    assert_true(u3.base_value == 1.0)

    # sqrt functions
    u4 = Unit("sqrt(g)/sqrt(cm)/s")
    assert_true(u4.dimensions == magnetic_field)
    assert_true(u4.base_value == 1.0)

    # commutative sqrt function
    u5 = Unit("sqrt(g/cm)/s")
    assert_true(u5.dimensions == magnetic_field)
    assert_true(u5.base_value == 1.0)

    # nonzero CGS conversion factor
    u6 = Unit("Msun/pc**3")
    assert_true(u6.dimensions == mass/length**3)
    assert_array_almost_equal_nulp(np.array([u6.base_value]), np.array([mass_sun_grams/cm_per_pc**3]))

    assert_raises(UnitParseError, Unit, 'm**m')
    assert_raises(UnitParseError, Unit, 'm**g')
    assert_raises(UnitParseError, Unit, 'm+g')
    assert_raises(UnitParseError, Unit, 'm-g')


def test_create_from_expr():
    """
    Create units from sympy Exprs and check attributes.

    """
    pc_cgs = cm_per_pc
    yr_cgs = sec_per_year

    # Symbol expr
    s1 = Symbol("pc", positive=True)
    s2 = Symbol("yr", positive=True)
    # Mul expr
    s3 = s1 * s2
    # Pow expr
    s4  = s1**2 * s2**(-1)

    u1 = Unit(s1)
    u2 = Unit(s2)
    u3 = Unit(s3)
    u4 = Unit(s4)

    assert_true(u1.expr == s1)
    assert_true(u2.expr == s2)
    assert_true(u3.expr == s3)
    assert_true(u4.expr == s4)

    assert_allclose_units(u1.base_value, pc_cgs, 1e-12)
    assert_allclose_units(u2.base_value, yr_cgs, 1e-12)
    assert_allclose_units(u3.base_value, pc_cgs * yr_cgs, 1e-12)
    assert_allclose_units(u4.base_value, pc_cgs**2 / yr_cgs, 1e-12)

    assert_true(u1.dimensions == length)
    assert_true(u2.dimensions == time)
    assert_true(u3.dimensions == length * time)
    assert_true(u4.dimensions == length**2 / time)


def test_create_with_duplicate_dimensions():
    """
    Create units with overlapping dimensions. Ex: km/Mpc.

    """

    u1 = Unit("erg * s**-1")
    u2 = Unit("km/s/Mpc")
    km_cgs = cm_per_km
    Mpc_cgs = cm_per_mpc

    assert_true(u1.base_value == 1)
    assert_true(u1.dimensions == power)

    assert_allclose_units(u2.base_value, km_cgs / Mpc_cgs, 1e-12)
    assert_true(u2.dimensions == rate)

def test_create_new_symbol():
    """
    Create unit with unknown symbol.

    """
    u1 = Unit("abc", base_value=42, dimensions=(mass/time))

    assert_true(u1.expr == Symbol("abc", positive=True))
    assert_true(u1.base_value == 42)
    assert_true(u1.dimensions == mass / time)

    u1 = Unit("abc", base_value=42, dimensions=length**3)

    assert_true(u1.expr == Symbol("abc", positive=True))
    assert_true(u1.base_value == 42)
    assert_true(u1.dimensions == length**3)

    u1 = Unit("abc", base_value=42, dimensions=length*(mass*length))

    assert_true(u1.expr == Symbol("abc", positive=True))
    assert_true(u1.base_value == 42)
    assert_true( u1.dimensions == length**2*mass)

    assert_raises(UnitParseError, Unit, 'abc', base_value=42,
                  dimensions=length**length)
    assert_raises(UnitParseError, Unit, 'abc', base_value=42,
                  dimensions=length**(length*length))
    assert_raises(UnitParseError, Unit, 'abc', base_value=42,
                  dimensions=length-mass)
    assert_raises(UnitParseError, Unit, 'abc', base_value=42,
                  dimensions=length+mass)

def test_create_fail_on_unknown_symbol():
    """
    Fail to create unit with unknown symbol, without base_value and dimensions.

    """
    try:
        Unit(Symbol("jigawatts"))
    except UnitParseError:
        assert_true(True)
    else:
        assert_true(False)

def test_create_fail_on_bad_symbol_type():
    """
    Fail to create unit with bad symbol type.

    """
    try:
        Unit([1])  # something other than Expr and str
    except UnitParseError:
        assert_true(True)
    else:
        assert_true(False)

def test_create_fail_on_bad_dimensions_type():
    """
    Fail to create unit with bad dimensions type.

    """
    try:
        Unit("a", base_value=1, dimensions="(mass)")
    except UnitParseError:
        assert_true(True)
    else:
        assert_true(False)


def test_create_fail_on_dimensions_content():
    """
    Fail to create unit with bad dimensions expr.

    """
    a = Symbol("a")

    try:
        Unit("a", base_value=1, dimensions=a)
    except UnitParseError:
        pass
    else:
        assert_true(False)


def test_create_fail_on_base_value_type():
    """
    Fail to create unit with bad base_value type.

    """
    try:
        Unit("a", base_value="a", dimensions=(mass/time))
    except UnitParseError:
        assert_true(True)
    else:
        assert_true(False)

#
# End init tests
#

def test_string_representation():
    """
    Check unit string representation.

    """
    pc = Unit("pc")
    Myr = Unit("Myr")
    speed = pc / Myr
    dimensionless = Unit()

    assert_true(str(pc) == "pc")
    assert_true(str(Myr) == "Myr")
    assert_true(str(speed) == "pc/Myr")
    assert_true(repr(speed) == "pc/Myr")
    assert_true(str(dimensionless) == "dimensionless")

#
# Start operation tests
#

def test_multiplication():
    """
    Multiply two units.

    """
    msun_cgs = mass_sun_grams
    pc_cgs = cm_per_pc

    # Create symbols
    msun_sym = Symbol("Msun", positive=True)
    pc_sym = Symbol("pc", positive=True)
    s_sym = Symbol("s", positive=True)

    # Create units
    u1 = Unit("Msun")
    u2 = Unit("pc")

    # Mul operation
    u3 = u1 * u2

    assert_true(u3.expr == msun_sym * pc_sym)
    assert_allclose_units(u3.base_value, msun_cgs * pc_cgs, 1e-12)
    assert_true(u3.dimensions == mass * length)

    # Pow and Mul operations
    u4 = Unit("pc**2")
    u5 = Unit("Msun * s")

    u6 = u4 * u5

    assert_true(u6.expr == pc_sym**2 * msun_sym * s_sym)
    assert_allclose_units(u6.base_value, pc_cgs**2 * msun_cgs, 1e-12)
    assert_true(u6.dimensions == length**2 * mass * time)


def test_division():
    """
    Divide two units.

    """
    pc_cgs = cm_per_pc
    km_cgs = cm_per_km

    # Create symbols
    pc_sym = Symbol("pc", positive=True)
    km_sym = Symbol("km", positive=True)
    s_sym = Symbol("s", positive=True)

    # Create units
    u1 = Unit("pc")
    u2 = Unit("km * s")

    u3 = u1 / u2

    assert_true(u3.expr == pc_sym / (km_sym * s_sym))
    assert_allclose_units(u3.base_value, pc_cgs / km_cgs, 1e-12)
    assert_true(u3.dimensions == 1 / time)


def test_power():
    """
    Take units to some power.

    """
    from sympy import nsimplify

    pc_cgs = cm_per_pc
    mK_cgs = 1e-3
    u1_dims = mass * length**2 * time**-3 * temperature**4
    u1 = Unit("g * pc**2 * s**-3 * mK**4")

    u2 = u1**2

    assert_true(u2.dimensions == u1_dims**2)
    assert_allclose_units(u2.base_value, (pc_cgs**2 * mK_cgs**4)**2, 1e-12)

    u3 = u1**(-1.0/3)

    assert_true(u3.dimensions == nsimplify(u1_dims**(-1.0/3)))
    assert_allclose_units(u3.base_value, (pc_cgs**2 * mK_cgs**4)**(-1.0/3), 1e-12)


def test_equality():
    """
    Check unit equality with different symbols, but same dimensions and base_value.

    """
    u1 = Unit("km * s**-1")
    u2 = Unit("m * ms**-1")

    assert_true(u1 == u2)

#
# End operation tests.
#

def test_base_equivalent():
    """
    Check base equivalent of a unit.

    """
    Msun_cgs = mass_sun_grams
    Mpc_cgs = cm_per_mpc

    u1 = Unit("Msun * Mpc**-3")
    u2 = Unit("g * cm**-3")
    u3 = u1.get_base_equivalent()

    assert_true(u2.expr == u3.expr)
    assert_true(u2 == u3)

    assert_allclose_units(u1.base_value, Msun_cgs / Mpc_cgs**3, 1e-12)
    assert_true(u2.base_value == 1)
    assert_true(u3.base_value == 1)

    mass_density = mass / length**3

    assert_true(u1.dimensions == mass_density)
    assert_true(u2.dimensions == mass_density)
    assert_true(u3.dimensions == mass_density)

    assert_allclose_units(get_conversion_factor(u1, u3)[0], Msun_cgs / Mpc_cgs**3, 1e-12)

def test_is_code_unit():
    ds = fake_random_ds(64, nprocs=1)
    u1 = Unit('code_mass', registry=ds.unit_registry)
    u2 = Unit('code_mass/code_length', registry=ds.unit_registry)
    u3 = Unit('code_velocity*code_mass**2', registry=ds.unit_registry)
    u4 = Unit('code_time*code_mass**0.5', registry=ds.unit_registry)
    u5 = Unit('code_mass*g', registry=ds.unit_registry)
    u6 = Unit('g/cm**3')

    assert_true(u1.is_code_unit)
    assert_true(u2.is_code_unit)
    assert_true(u3.is_code_unit)
    assert_true(u4.is_code_unit)
    assert_true(not u5.is_code_unit)
    assert_true(not u6.is_code_unit)

def test_temperature_offsets():
    u1 = Unit('degC')
    u2 = Unit('degF')

    assert_raises(InvalidUnitOperation, operator.mul, u1, u2)
    assert_raises(InvalidUnitOperation, operator.truediv, u1, u2)

def test_latex_repr():
    ds = fake_random_ds(64, nprocs=1)

    # create a fake comoving unit
    ds.unit_registry.add('pccm', ds.unit_registry.lut['pc'][0]/(1+2), length,
                         "\\rm{pc}/(1+z)")

    test_unit = Unit('Mpccm', registry=ds.unit_registry)
    assert_almost_equal(test_unit.base_value, cm_per_mpc/3)
    assert_equal(test_unit.latex_repr, r'\rm{Mpc}/(1+z)')

    test_unit = Unit('code_mass', registry=ds.unit_registry)
    assert_equal(test_unit.latex_repr, '\\rm{code\\ mass}')

    test_unit = Unit('code_mass/code_length**3', registry=ds.unit_registry)
    assert_equal(test_unit.latex_repr,
                 '\\frac{\\rm{code\\ mass}}{\\rm{code\\ length}^{3}}')

    test_unit = Unit('cm**-3', base_value=1.0, registry=ds.unit_registry)
    assert_equal(test_unit.latex_repr, '\\frac{1}{\\rm{cm}^{3}}')

    test_unit = Unit('m_geom/l_geom**3')
    assert_equal(test_unit.latex_repr, '\\frac{1}{M_\\odot^{2}}')

    test_unit = Unit('1e9*cm')
    assert_equal(test_unit.latex_repr, '1.0 \\times 10^{9}\\ \\rm{cm}')

def test_latitude_longitude():
    lat = unit_symbols.lat
    lon = unit_symbols.lon
    deg = unit_symbols.deg
    assert_equal(lat.units.base_offset, 90.0)
    assert_equal((deg*90.0).in_units("lat").value, 0.0)
    assert_equal((deg*180).in_units("lat").value, -90.0)
    assert_equal((lat*0.0).in_units("deg"), deg*90.0)
    assert_equal((lat*-90).in_units("deg"), deg*180)

    assert_equal(lon.units.base_offset, -180.0)
    assert_equal((deg*0.0).in_units("lon").value, -180.0)
    assert_equal((deg*90.0).in_units("lon").value, -90.0)
    assert_equal((deg*180).in_units("lon").value, 0.0)
    assert_equal((deg*360).in_units("lon").value, 180.0)

    assert_equal((lon*-180.0).in_units("deg"), deg*0.0)
    assert_equal((lon*-90.0).in_units("deg"), deg*90.0)
    assert_equal((lon*0.0).in_units("deg"), deg*180.0)
    assert_equal((lon*180.0).in_units("deg"), deg*360)

def test_registry_json():
    reg = UnitRegistry()
    json_reg = reg.to_json()
    unserialized_reg = UnitRegistry.from_json(json_reg)

    assert_equal(reg.lut, unserialized_reg.lut)

def test_creation_from_ytarray():
    u1 = Unit(electrostatic_unit)
    assert_equal(str(u1), 'esu')
    assert_equal(u1, Unit('esu'))
    assert_equal(u1, electrostatic_unit.units)

    u2 = Unit(elementary_charge)
    assert_equal(str(u2), '4.8032056e-10*esu')
    assert_equal(u2, Unit('4.8032056e-10*esu'))
    assert_equal(u1, elementary_charge.units)

    assert_equal((u1/u2).base_value, electrostatic_unit/elementary_charge)

    assert_raises(UnitParseError, Unit, [1, 2, 3]*elementary_charge)

def test_list_same_dimensions():
    reg = default_unit_registry
    for name1, u1 in reg.unit_objs.items():
        for name2 in reg.list_same_dimensions(u1):
            if name2 == name1: continue
            if name2 in reg.unit_objs:
                dim2 = reg.unit_objs[name2].dimensions
            else:
                _, dim2, _, _ = reg.lut[name2]
            assert_true(u1.dimensions is dim2)
