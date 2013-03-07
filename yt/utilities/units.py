"""
Symbolic unit handling.

Author: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley

Homepage: http://yt-project.org/
License:
  Copyright (C) 2012, 2013 Casey W. Stark.  All Rights Reserved.

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

from sympy import Expr, Mul, Number, Pow, Rational, Symbol
from sympy import nsimplify, posify, sympify
from sympy.parsing.sympy_parser import parse_expr

# The base dimensions
mass = Symbol("(mass)", positive=True)
length = Symbol("(length)", positive=True)
time = Symbol("(time)", positive=True)
temperature = Symbol("(temperature)", positive=True)
base_dimensions = [mass, length, time, temperature]

### Misc. dimensions
rate = 1 / time

dimensionless = sympify(1)

velocity     = length / time
acceleration = length / time**2
jerk         = length / time**3
snap         = length / time**4
crackle      = length / time**5
pop          = length / time**6

momentum = mass * velocity
force    = mass * acceleration
energy   = force * length
power    = energy / time
charge   = (energy * length)**Rational(1, 2)  # proper 1/2 power

electric_field = charge / length**2
magnetic_field = electric_field

# The key is the symbol, the value is a tuple with the conversion factor to
# cgs, the dimensionality, and
default_unit_symbol_LUT = {
    # base
    "g":  (1.0, mass, r"\rm{g}"),
    #"cm": (1.0, length, r"\rm{cm}"),  # duplicate with meter below...
    "s":  (1.0, time, r"\rm{s}"),
    "K":  (1.0, temperature, r"\rm{K}"),

    # other cgs
    "dyne": (1.0, force, r"\rm{dyne}"),
    "erg":  (1.0, energy, r"\rm{erg}"),
    "esu":  (1.0, charge, r"\rm{esu}"),

    # some SI
    "m": (1.0e2, length, r"\rm{m}"),
    "J": (1.0e7, energy, r"\rm{J}"),
    "Hz": (1.0, rate, r"\rm{Hz}"),

    # times
    "min": (60.0, time, r"\rm{min}"),
    "hr":  (3600.0, time, r"\rm{hr}"),
    "day": (86400.0, time, r"\rm{day}"),
    "yr":  (31536000.0, time, r"\rm{yr}"),

    # Solar units
    "Msun": (1.98892e33, mass, r"M_{\odot}"),
    "Rsun": (6.96e10, length, r"R_{\odot}"),
    "Lsun": (3.9e33, power, r"L_{\odot}"),
    "Tsun": (5870.0, temperature, r"T_{\odot}"),
    "Zsun": (1.0, dimensionless, r"Z_{\odot}"),

    # astro distances
    "AU": (1.49598e13, length, r"\rm{AU}"),
    "ly": (9.46053e17, length, r"\rm{ly}"),
    "pc": (3.08568e18, length, r"\rm{pc}"),

    # other astro
    "H_0": (2.3e-18, rate, r"H_0"),  # check cf

    # other energy units
    "eV": (1.6021766e-12, energy, r"\rm{eV}"),

    # electric stuff
    "gauss": (1.0, magnetic_field, r"\rm{G}"),
}

# This dictionary formatting from magnitude package, credit to Juan Reyero.
unit_prefixes = {
    'Y': 1e24,   # yotta
    'Z': 1e21,   # zetta
    'E': 1e18,   # exa
    'P': 1e15,   # peta
    'T': 1e12,   # tera
    'G': 1e9,    # giga
    'M': 1e6,    # mega
    'k': 1e3,    # kilo
    'd': 1e1,    # deci
    'c': 1e-2,   # centi
    'm': 1e-3,   # mili
    'u': 1e-6,   # micro
    'n': 1e-9,   # nano
    'p': 1e-12,  # pico
    'f': 1e-15,  # femto
    'a': 1e-18,  # atto
    'z': 1e-21,  # zepto
    'y': 1e-24,  # yocto
}


class UnitSymbolRegistry:

    def __init__(self, add_default_symbols=True):
        self.lookup_dict = {}

        if add_default_symbols:
            self.lookup_dict.update(default_symbol_lookup_dict)




class UnitParseError(Exception):
    pass

class UnitOperationError(Exception):
    pass


class Unit(Expr):
    """
    A symbolic unit, using sympy functionality. We only add "dimensions" so that
    sympy understands relations between different units.

    """

    is_positive = True    # make sqrt(m**2) --> m
    is_commutative = True
    is_number = False

    __slots__ = ["expr", "cgs_value", "dimensions", "is_atomic"]

    def __new__(cls, unit_expr=sympify(1), cgs_value=None, dimensions=None,
                **assumptions):
        """
        Create a new unit. May be an atomic unit (like a gram) or a combination
        of other units (like g / cm**3).

        Parameters
        ----------
        unit_expr : string or sympy.core.expr.Expr
            The unit symbol expression.
        cgs_value : float
            This unit's value in cgs.
        dimensions : sympy.core.expr.Expr
            An expression representing the dimensionality of this unit. This
            should be a product (sympy.core.mul.Mul object) of mass, length,
            time, and temperature objects to various powers.

        """
        # if we have a string, parse into an expression
        if isinstance(unit_expr, str):
            if not unit_expr:
                # Bug catch...
                # if unit_expr is an empty string, parse_expr fails hard...
                unit_expr = "1"
            unit_expr = parse_expr(unit_expr)

        if not isinstance(unit_expr, Expr):
            raise UnitParseError("Unit representation must be a string or sympy Expr. %s is a %s." % (unit_expr, type(unit_expr)))
        # done with argument checking...

        # sympify, make positive symbols, and nsimplify the expr
        unit_expr = sympify(unit_expr)
        unit_expr = make_symbols_positive(unit_expr)
        unit_expr = nsimplify(unit_expr)

        # see if the unit is atomic.
        is_atomic = False
        if isinstance(unit_expr, Symbol):
            is_atomic = True

        # Did the user supply cgs_value and dimensions?
        if cgs_value and not dimensions or dimensions and not cgs_value:
            raise Exception("If you provide cgs_vale or dimensions, you must provide both. cgs_value is %s, dimensions is %s." % (cgs_value, dimensions))

        if cgs_value and dimensions:
            # check that cgs_value is a float or can be converted to one
            try:
                cgs_value = float(cgs_value)
            except ValueError:
                raise UnitParseError("Please provide a float for the cgs_value kwarg. I got '%s'." % cgs_value)
            # check that dimensions is valid
            dimensions = verify_dimensions(dimensions)
            # save the values
            this_cgs_value, this_dimensions = cgs_value, dimensions

        else:  # lookup the unit symbols
            this_cgs_value, this_dimensions = \
                get_unit_data_from_expr(unit_expr)

        # Trick to get dimensions powers as Rationals
        this_dimensions = nsimplify(this_dimensions)

        # create obj with superclass construct
        obj = Expr.__new__(cls, **assumptions)

        # attach attributes to obj
        obj.expr = unit_expr
        obj.is_atomic = is_atomic
        obj.cgs_value = this_cgs_value
        obj.dimensions = this_dimensions

        # return `obj` so __init__ can handle it.
        return obj

    ### some sympy conventions I guess
    def __getnewargs__(self):
        return (self.expr, self.is_atomic, self.cgs_value, self.dimensions)

    def __hash__(self):
        return super(Unit, self).__hash__()

    def _hashable_content(self):
        return (self.expr, self.is_atomic, self.cgs_value, self.dimensions)
    ### end sympy conventions

    def __repr__(self):
        if self.expr == 1:
            return "(dimensionless)"
        return str(self.expr)

    def __str__(self):
        if self.expr == 1:
            return "(dimensionless)"
        return str(self.expr)

    # for sympy.printing
    def _sympystr(self, *args):
        return str(self.expr)

    ### override sympy operations
    def __mul__(self, u):
        """ Multiply Unit with u (Unit object). """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to multiply Unit object with '%s'. This behavior is undefined." % u)

        return Unit(self.expr * u.expr, self.cgs_value * u.cgs_value,
                    self.dimensions * u.dimensions)

    def __div__(self, u):
        """ Divide Unit by u (Unit object). """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to divide Unit object by '%s'. This behavior is undefined." % u)

        return Unit(self.expr / u.expr, self.cgs_value / u.cgs_value,
                    self.dimensions / u.dimensions)

    def __pow__(self, p):
        """ Take Unit to power p (float). """
        try:
            p = sympify(p)
        except ValueError:
            raise UnitOperationError("Tried to take Unit object to the power '%s'. I could not cast this to a float." % p)

        return Unit(self.expr**p, self.cgs_value**p, self.dimensions**p)

    ### Comparison operators
    def same_dimensions_as(self, other_unit):
        """ Test if dimensions are the same. """
        return (self.dimensions / other_unit.dimensions) == 1

    def __eq__(self, u):
        """ Test unit equality. """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to test equality between Unit object and '%s'. This behavior is undefined." % u)

        return (self.cgs_value == u.cgs_value and self.dimensions == u.dimensions)

    @property
    def is_dimensionless(self):
        return self.dimensions == 1

    def get_cgs_equivalent(self):
        """ Create and return dimensionally-equivalent cgs units. """
        cgs_units_string = "g**(%s) * cm**(%s) * s**(%s) * K**(%s)" % \
            (self.dimensions.expand().as_coeff_exponent(mass)[1],
             self.dimensions.expand().as_coeff_exponent(length)[1],
             self.dimensions.expand().as_coeff_exponent(time)[1],
             self.dimensions.expand().as_coeff_exponent(temperature)[1])
        return Unit(cgs_units_string, 1, self.dimensions)

    def get_conversion_factor(self, other_units):
        return get_conversion_factor(self, other_units)


def make_symbols_positive(expr):
    """
    Grabs all symbols from expr, makes new positive symbols with the same names,
    and substitutes them back into the expression.

    """
    expr_symbols = expr.atoms(Symbol)  # grab all symbols

    # Replace one at a time
    for s in expr_symbols:
        # replace this symbol with a positive version
        expr = expr.subs(s, Symbol(s.name, positive=True))

    return expr


# @todo: simpler method that doesn't use recursion would be better...
# We could check if dimensions.atoms are all numbers or symbols, but we should
# check the functions also...
def verify_dimensions(d):
    """
    Make sure that `d` is a valid dimension expression. It must consist of only
    the base dimension symbols, to powers, multiplied together. If valid, return
    the simplified expression. If not, raise an Exception.

    """
    # in the case of a Number of Symbol, we can just return
    if isinstance(d, Number):
        return d
    elif isinstance(d, Symbol):
        if d in base_dimensions:
            return d
        else:
            raise UnitParseError("dimensionality expression contains an unknown symbol '%s'." % d)

    # validate args of a Pow or Mul separately
    elif isinstance(d, Pow):
        return verify_dimensions(d.args[0])**verify_dimensions(d.args[1])

    elif isinstance(d, Mul):
        total_mul = 1
        for arg in d.args:
            total_mul *= verify_dimensions(arg)
        return total_mul

    # should never get here
    raise UnitParseError("Bad dimensionality expression '%s'." % d)


def get_unit_data_from_expr(unit_expr):
    """
    Gets total cgs_value and dimensions from a valid unit expression.

    """
    # The simplest case first
    if isinstance(unit_expr, Unit):
        return (unit_expr.cgs_value, unit_expr.dimensions)

    # Now for the sympy possibilities
    if isinstance(unit_expr, Symbol):
        return lookup_unit_symbol(str(unit_expr))

    if isinstance(unit_expr, Number):
        return (1, 1)

    if isinstance(unit_expr, Pow):
        unit_data = get_unit_data_from_expr(unit_expr.args[0])
        power = unit_expr.args[1]
        return (unit_data[0]**power, unit_data[1]**power)

    if isinstance(unit_expr, Mul):
        cgs_value = 1
        dimensions = 1
        for expr in unit_expr.args:
            unit_data = get_unit_data_from_expr(expr)
            cgs_value *= unit_data[0]
            dimensions *= unit_data[1]

        return (cgs_value, dimensions)

    raise UnitParseError("Cannot parse for unit data from '%s'. Please supply an expression of only Unit, Symbol, Pow, and Mul objects." % str(unit_expr))


def lookup_unit_symbol(symbol_str):
    """ Searches for the unit data typle corresponding to the given symbol. """

    if symbol_str in default_unit_symbols_LUT:
        # lookup successful, return the tuple directly
        return default_unit_symbols_LUT[symbol_str]

    # could still be a known symbol with a prefix
    possible_prefix = symbol_str[0]
    if possible_prefix in unit_prefixes:
        # the first character could be a prefix, check the rest of the symbol
        symbol_wo_prefix = symbol_str[1:]

        if symbol_wo_prefix in default_unit_symbols_LUT:
            # lookup successful, it's a symbol with a prefix
            unit_data = default_unit_symbols_LUT[symbol_wo_prefix]
            prefix_value = unit_prefixes[possible_prefix]

            # don't forget to account for the prefix value!
            return (unit_data[0] * prefix_value, unit_data[1])

    # no dice
    raise UnitParseError("Could not find unit symbol '%s'. Please supply the dimensions and cgs value when creating this object." % symbol_str)


def get_conversion_factor(old_units, new_units):
    """
    Get the conversion factor between two units of equivalent dimensions. This
    is the number you multiply data by to convert from values in `old_units` to
    values in `new_units`.

    Parameters
    ----------
    old_units: str or Unit object
        The current units.
    new_units : str or Unit object
        The units we want.

    Returns
    -------
    conversion_factor : float
        `old_units / new_units`

    """
    # Make units out of strings if we need to.
    if isinstance(old_units, str):
        old_units = Unit(old_units)
    if isinstance(new_units, str):
        new_units = Unit(new_units)

    if not old_units.same_dimensions_as(new_units):
        raise UnitOperationError("Cannot convert from %s to %s because the dimensions do not match: %s and %s" % (old_units, new_units, old_units.dimensions, new_units.dimensions))

    return old_units.cgs_value / new_units.cgs_value


def convert_values(values, old_units, new_units):
    """
    Take data given in old units and convert to values in new units.

    Parameters
    ----------
    values : array_like
        The number or array we will convert.
    old_units : str or Unit object
        The units values are supplied in.
    new_units : str or Unit object
        The units values will be returned in.

    Returns values in new units.

    """
    return values * get_conversion_factor(old_units, new_units)
