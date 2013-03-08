"""
Test symbolic unit handling.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley

License:
  Copyright (C) 2012 Casey W. Stark.  All Rights Reserved.

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

import nose
from numpy.testing import assert_approx_equal
from sympy import Symbol

# base dimensions
from yt.utilities.units import mass, length, time, temperature
# functions
from yt.utilities.units import get_conversion_factor, make_symbols_positive
# classes
from yt.utilities.units import Unit, UnitParseError
# objects
from yt.utilities.units import unit_symbols_dict, unit_prefixes


def test_no_conflicting_symbols():
    """
    Check unit symbol definitions for conflicts.

    """
    full_set = set(unit_symbols_dict.keys())

    # go through all possible prefix combos
    for symbol in unit_symbols_dict.keys():
        for prefix in unit_prefixes.keys():
            new_symbol = "%s%s" % (prefix, symbol)

            # test if we have seen this symbol
            if new_symbol in full_set:
                print "Duplicate symbol: %s" % new_symbol
                assert False

            full_set.add(new_symbol)


def test_dimensionless():
    """
    Create dimensionless unit and check attributes.

    """
    u1 = Unit()

    assert u1.is_dimensionless
    assert u1.expr == 1
    assert u1.cgs_value == 1
    assert u1.dimensions == 1

#
# Start init tests
#

def test_create_from_string():
    """
    Create units with strings and check attributes.

    """
    from yt.utilities.units import energy

    u1 = Unit("g * cm**2 * s**-2")
    assert u1.dimensions == energy
    assert u1.cgs_value == 1.0

    # make sure order doesn't matter
    u2 = Unit("cm**2 * s**-2 * g")
    assert u2.dimensions == energy
    assert u2.cgs_value == 1.0


def test_create_from_expr():
    """
    Create units from sympy Exprs and check attributes.

    """
    pc_cgs = 3.08568e18
    yr_cgs = 31536000

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

    assert u1.expr == s1
    assert u2.expr == s2
    assert u3.expr == s3
    assert u4.expr == s4

    assert u1.cgs_value == pc_cgs
    assert u2.cgs_value == yr_cgs
    assert u3.cgs_value == pc_cgs * yr_cgs
    assert u4.cgs_value == pc_cgs**2 / yr_cgs

    assert u1.dimensions == length
    assert u2.dimensions == time
    assert u3.dimensions == length * time
    assert u4.dimensions == length**2 / time


def test_create_with_duplicate_dimensions():
    """
    Create units with overlapping dimensions. Ex: km/Mpc.

    """
    from yt.utilities.units import power, rate

    u1 = Unit("erg * s**-1")
    u2 = Unit("km/s/Mpc")
    km_cgs = 1e5
    Mpc_cgs = 3.08568e24

    print u2.expr
    print u2.dimensions
    print u2.cgs_value

    assert u1.cgs_value == 1
    assert u1.dimensions == power

    assert u2.cgs_value == km_cgs / Mpc_cgs
    assert u2.dimensions == rate


def test_create_new_symbol():
    """
    Create unit with unknown symbol.

    """
    u1 = Unit("abc", cgs_value=42, dimensions=(mass/time))

    assert u1.expr == Symbol("abc", positive=True)
    assert u1.cgs_value == 42
    assert u1.dimensions == mass / time


def test_create_fail_on_unknown_symbol():
    """
    Fail to create unit with unknown symbol, without cgs_value and dimensions.

    """
    try:
        u1 = Unit(Symbol("jigawatts"))
    except UnitParseError:
        pass
    else:
        assert False


def test_create_fail_on_bad_symbol_type():
    """
    Fail to create unit with bad symbol type.

    """
    try:
        u1 = Unit([1])  # something other than Expr and str
    except UnitParseError:
        pass
    else:
        assert False


def test_create_fail_on_bad_dimensions_type():
    """
    Fail to create unit with bad dimensions type.

    """
    try:
        u1 = Unit("a", cgs_value=1, dimensions="(mass)")
    except UnitParseError:
        pass
    else:
        assert False


def test_create_fail_on_dimensions_content():
    """
    Fail to create unit with bad dimensions expr.

    """
    a = Symbol("a")

    try:
        u1 = Unit("a", cgs_value=1, dimensions=a)
    except UnitParseError:
        pass
    else:
        assert False


def test_create_fail_on_cgs_value_type():
    """
    Fail to create unit with bad cgs_value type.

    """
    try:
        u1 = Unit("a", cgs_value="a", dimensions=(mass/time))
    except UnitParseError:
        pass
    else:
        assert False

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

    assert str(pc) == "pc"
    assert str(Myr) == "Myr"
    assert str(speed) == "pc/Myr"
    assert repr(speed) == "pc/Myr"
    assert str(dimensionless) == "(dimensionless)"

#
# Start operation tests
#

def test_multiplication():
    """
    Multiply two units.

    """
    msun_cgs = 1.98892e33
    pc_cgs = 3.08568e18

    # Create symbols
    msun_sym = Symbol("Msun", positive=True)
    pc_sym = Symbol("pc", positive=True)
    s_sym = Symbol("s", positive=True)

    # Create units
    u1 = Unit("Msun")
    u2 = Unit("pc")

    # Mul operation
    u3 = u1 * u2

    assert u3.expr == msun_sym * pc_sym
    assert u3.cgs_value == msun_cgs * pc_cgs
    assert u3.dimensions == mass * length

    # Pow and Mul operations
    u4 = Unit("pc**2")
    u5 = Unit("Msun * s")

    u6 = u4 * u5

    assert u6.expr == pc_sym**2 * msun_sym * s_sym
    assert u6.cgs_value == pc_cgs**2 * msun_cgs
    assert u6.dimensions == length**2 * mass * time


def test_division():
    """
    Divide two units.

    """
    pc_cgs = 3.08568e18
    km_cgs = 1e5

    # Create symbols
    pc_sym = Symbol("pc", positive=True)
    km_sym = Symbol("km", positive=True)
    s_sym = Symbol("s", positive=True)

    # Create units
    u1 = Unit("pc")
    u2 = Unit("km * s")

    u3 = u1 / u2

    assert u3.expr == pc_sym / (km_sym * s_sym)
    assert u3.cgs_value == pc_cgs / km_cgs
    assert u3.dimensions == 1 / time


def test_power():
    """
    Take units to some power.

    """
    from sympy import nsimplify

    pc_cgs = 3.08568e18
    mK_cgs = 1e-3
    u1_dims = mass * length**2 * time**-3 * temperature**4
    u1 = Unit("g * pc**2 * s**-3 * mK**4")

    u2 = u1**2

    assert u2.dimensions == u1_dims**2
    assert u2.cgs_value == (pc_cgs**2 * mK_cgs**4)**2

    u3 = u1**(-1.0/3)

    assert u3.dimensions == nsimplify(u1_dims**(-1.0/3))
    assert u3.cgs_value == (pc_cgs**2 * mK_cgs**4)**(-1.0/3)


def test_equality():
    """
    Check unit equality with different symbols, but same dimensions and cgs_value.

    """
    u1 = Unit("km * s**-1")
    u2 = Unit("m * ms**-1")

    assert u1 == u2

#
# End operation tests.
#

def test_cgs_equivalent():
    """
    Check cgs equivalent of a unit.

    """
    Msun_cgs = 1.98892e33
    Mpc_cgs = 3.08568e24

    u1 = Unit("Msun * Mpc**-3")
    u2 = Unit("g * cm**-3")
    u3 = u1.get_cgs_equivalent()

    assert u2.expr == u3.expr
    assert u2 == u3

    assert_approx_equal( u1.cgs_value, Msun_cgs / Mpc_cgs**3, significant=16)
    assert u2.cgs_value == 1
    assert u3.cgs_value == 1

    mass_density = mass / length**3

    assert u1.dimensions == mass_density
    assert u2.dimensions == mass_density
    assert u3.dimensions == mass_density

    assert_approx_equal( get_conversion_factor(u1, u3), Msun_cgs / Mpc_cgs**3,
                         significant=16 )
