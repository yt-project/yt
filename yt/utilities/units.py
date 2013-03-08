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

# Define a sympy one object.
sympy_one = sympify(1)

#
# Exceptions
#

class SymbolNotFoundError(Exception):
    pass

class UnitParseError(Exception):
    pass

class UnitOperationError(Exception):
    pass

#
# Base dimensions
#

mass = Symbol("(mass)", positive=True)
length = Symbol("(length)", positive=True)
time = Symbol("(time)", positive=True)
temperature = Symbol("(temperature)", positive=True)

base_dimensions = [mass, length, time, temperature]

#
# Derived dimensions
#

dimensionless = sympify(1)

rate = 1 / time

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

#
# The default/basic unit symbol lookup table.
#
# Lookup a unit symbol with the symbol string, and provide a tuple with the
# conversion factor to cgs and dimensionality.
#

default_unit_symbol_lut = {
    # base
    "g":  (1.0, mass),
    #"cm": (1.0, length, r"\rm{cm}"),  # duplicate with meter below...
    "s":  (1.0, time),
    "K":  (1.0, temperature),

    # other cgs
    "dyne": (1.0, force),
    "erg":  (1.0, energy),
    "esu":  (1.0, charge),

    # some SI
    "m": (1.0e2, length),
    "J": (1.0e7, energy),
    "Hz": (1.0, rate),

    # dimensionless stuff
    "h": (1.0, dimensionless),

    # times
    "min": (60.0, time),
    "hr":  (3600.0, time),
    "day": (86400.0, time),
    "yr":  (31536000.0, time),

    # Solar units
    "Msun": (1.98892e33, mass),
    "msun": (1.98892e33, mass),
    "Rsun": (6.96e10, length),
    "Lsun": (3.9e33, power),
    "Tsun": (5870.0, temperature),
    "Zsun": (1.0, dimensionless),

    # astro distances
    "AU": (1.49598e13, length),
    "ly": (9.46053e17, length),
    "pc": (3.08568e18, length),

    # electric stuff
    "gauss": (1.0, magnetic_field),

    # other astro
    "H_0": (2.3e-18, rate),  # check cf

    # misc
    "eV": (1.6021766e-12, energy),
    "amu": (1.6605402e-24, mass),

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


class UnitRegistry:

    def __init__(self, add_default_symbols=True):
        self.lut = {}

        if add_default_symbols:
            self.lut.update(default_unit_symbol_lut)

    def add(self, symbol, cgs_value, dimensions):
        """
        Add a symbol to this registry.

        """
        # Validate
        if not isinstance(cgs_value, float):
            raise UnitParseError("cgs_value must be a float, got a %s." % type(cgs_value))

        _validate_dimensions(dimensions)

        # Add to lut
        self.lut.update( {symbol: (cgs_value, dimensions)} )

    def remove(self, symbol):
        """
        Remove the entry for the unit matching `symbol`.

        """
        if symbol not in self.lut:
            raise SymbolNotFoundError("Tried to remove the symbol '%s', but it does not exist in this registry." % symbol)

        del self.lut[symbol]


default_unit_registry = UnitRegistry()

class Unit(Expr):
    """
    A symbolic unit, using sympy functionality. We only add "dimensions" so
    that sympy understands relations between different units.

    """

    # Set some assumptions for sympy.
    is_positive = True    # make sqrt(m**2) --> m
    is_commutative = True
    is_number = False

    # Extra attributes
    __slots__ = ["expr", "is_atomic", "cgs_value", "dimensions", "registry"]

    def __new__(cls, unit_expr=sympify(1), cgs_value=None, dimensions=None,
                registry=None, **assumptions):
        """
        Create a new unit. May be an atomic unit (like a gram) or combinations
        of atomic units (like g / cm**3).

        Parameters
        ----------
        unit_expr : Unit object, sympy.core.expr.Expr object, or str
            The symbolic unit expression.
        cgs_value : float
            The unit's value in cgs.
        dimensions : sympy.core.expr.Expr
            A sympy expression representing the dimensionality of this unit.
            It must contain only mass, length, time, and temperature symbols.
        registry : UnitRegistry object
            The unit registry we use to interpret unit symbols.

        """
        # Simplest case. If user passes a Unit object, just use the expr.
        if isinstance(unit_expr, Unit):
            # grab the unit object's sympy expression.
            unit_expr = unit_expr.expr
        # If we have a string, have sympy parse it into an Expr.
        elif isinstance(unit_expr, str):
            if not unit_expr:
                # Bug catch...
                # if unit_expr is an empty string, parse_expr fails hard...
                unit_expr = "1"
            unit_expr = parse_expr(unit_expr)
        # Make sure we have an Expr at this point.
        if not isinstance(unit_expr, Expr):
            raise UnitParseError("Unit representation must be a string or sympy Expr. %s has type %s." % (unit_expr, type(unit_expr)))

        if registry is None:
            # Caller did not set the registry, so use the default.
            registry = default_unit_registry

        # done with argument checking...

        # sympify, make positive symbols, and nsimplify the expr
        unit_expr = sympify(unit_expr)
        unit_expr = _make_symbols_positive(unit_expr)
        unit_expr = nsimplify(unit_expr)

        # see if the unit is atomic.
        is_atomic = False
        if isinstance(unit_expr, Symbol):
            is_atomic = True

        #
        # check cgs_value and dimensions
        #

        if cgs_value is not None and dimensions is not None:
            # check that cgs_value is a float or can be converted to one
            try:
                cgs_value = float(cgs_value)
            except ValueError:
                raise UnitParseError("Could not use cgs_value as a float. cgs_value is '%s' (type %s)." % (cgs_value, type(cgs_value)) )

            # check that dimensions is valid
            _validate_dimensions( sympify(dimensions) )
        else:
            # lookup the unit symbols
            cgs_value, dimensions = _get_unit_data_from_expr(unit_expr,
                                                             registry.lut)

        # Sympy trick to get dimensions powers as Rationals
        dimensions = nsimplify(dimensions)

        # Create obj with superclass construct.
        obj = Expr.__new__(cls, **assumptions)

        # Attach attributes to obj.
        obj.expr = unit_expr
        obj.is_atomic = is_atomic
        obj.cgs_value = cgs_value
        obj.dimensions = dimensions
        obj.registry = registry

        # Return `obj` so __init__ can handle it.
        return obj

    ### Some sympy conventions
    def __getnewargs__(self):
        return (self.expr, self.is_atomic, self.cgs_value, self.dimensions,
                self.registry)

    def __hash__(self):
        return super(Unit, self).__hash__()

    def _hashable_content(self):
        return (self.expr, self.is_atomic, self.cgs_value, self.dimensions,
                self.registry)
    ### end sympy conventions

    def __repr__(self):
        if self.expr == sympy_one:
            return "(dimensionless)"
        # @todo: don't use dunder method?
        return self.expr.__repr__()

    def __str__(self):
        if self.expr == sympy_one:
            return "(dimensionless)"
        # @todo: don't use dunder method?
        return self.expr.__str__()

    # for sympy.printing
    def _sympystr(self, *args):
        return str(self.expr)

    #
    # Start unit operations
    #

    def __mul__(self, u):
        """ Multiply Unit with u (Unit object). """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to multiply a Unit object with '%s' (type %s). This behavior is undefined." % (u, type(u)) )

        return Unit(self.expr * u.expr,
                    cgs_value=(self.cgs_value * u.cgs_value),
                    dimensions=(self.dimensions * u.dimensions),
                    registry=self.registry)

    def __div__(self, u):
        """ Divide Unit by u (Unit object). """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to divide a Unit object by '%s' (type %s). This behavior is undefined." % (u, type(u)) )

        return Unit(self.expr / u.expr,
                    cgs_value=(self.cgs_value / u.cgs_value),
                    dimensions=(self.dimensions / u.dimensions),
                    registry=self.registry)

    def __pow__(self, p):
        """ Take Unit to power p (float). """
        try:
            p = sympify(p)
        except ValueError:
            raise UnitOperationError("Tried to take a Unit object to the power '%s' (type %s). Failed to cast it to a float." % (p, type(p)) )

        return Unit(self.expr**p, cgs_value=(self.cgs_value**p),
                    dimensions=(self.dimensions**p), registry=self.registry)

    def __eq__(self, u):
        """ Test unit equality. """
        if not isinstance(u, Unit):
            raise UnitOperationError("Tried to test equality between a Unit object and '%s' (type %s). This behavior is undefined." % (u, type(u)) )

        return (self.cgs_value == u.cgs_value and self.dimensions == u.dimensions)

    #
    # End unit operations
    #

    def same_dimensions_as(self, other_unit):
        """ Test if dimensions are the same. """
        return (self.dimensions / other_unit.dimensions) == sympy_one

    @property
    def is_dimensionless(self):
        return self.dimensions == sympy_one

    # @todo: might be a simpler/smarter sympy way to do this...
    def get_cgs_equivalent(self):
        """
        Create and return dimensionally-equivalent cgs units.

        """
        cgs_units_string = "g**(%s) * cm**(%s) * s**(%s) * K**(%s)" % \
            (self.dimensions.expand().as_coeff_exponent(mass)[1],
             self.dimensions.expand().as_coeff_exponent(length)[1],
             self.dimensions.expand().as_coeff_exponent(time)[1],
             self.dimensions.expand().as_coeff_exponent(temperature)[1])

        return Unit(cgs_units_string, cgs_value=1.0,
                    dimensions=self.dimensions, registry=self.registry)

    def get_conversion_factor(self, other_units):
        return get_conversion_factor(self, other_units)

#
# Unit manipulation functions
#

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
    # if args are not Unit objects, construct them
    if not isinstance(old_units, Unit):
        old_units = Unit(old_units)
    if not isinstance(new_units, Unit):
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


#
# Helper functions
#

def _make_symbols_positive(expr):
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
def _validate_dimensions(d):
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
        return _validate_dimensions(d.args[0])**_validate_dimensions(d.args[1])

    elif isinstance(d, Mul):
        total_mul = 1
        for arg in d.args:
            total_mul *= _validate_dimensions(arg)
        return total_mul

    # should never get here
    raise UnitParseError("Bad dimensionality expression '%s'." % d)


def _get_unit_data_from_expr(unit_expr, unit_symbol_lut):
    """
    Grabs the total cgs_value and dimensions from a valid unit expression.

    Parameters
    ----------
    unit_expr: Unit object, or sympy Expr object
        The expression containing unit symbols.
    unit_symbol_lut: dict
        Provides the unit data for each valid unit symbol.

    """
    # The simplest case first
    if isinstance(unit_expr, Unit):
        return (unit_expr.cgs_value, unit_expr.dimensions)

    # Now for the sympy possibilities
    if isinstance(unit_expr, Symbol):
        return _lookup_unit_symbol(str(unit_expr), unit_symbol_lut)

    if isinstance(unit_expr, Number):
        # not sure if this should be (1, 1)...
        return (float(unit_expr), sympy_one)

    if isinstance(unit_expr, Pow):
        unit_data = _get_unit_data_from_expr(unit_expr.args[0], unit_symbol_lut)
        power = unit_expr.args[1]
        return (unit_data[0]**power, unit_data[1]**power)

    if isinstance(unit_expr, Mul):
        cgs_value = 1.0
        dimensions = 1
        for expr in unit_expr.args:
            unit_data = _get_unit_data_from_expr(expr, unit_symbol_lut)
            cgs_value *= unit_data[0]
            dimensions *= unit_data[1]

        return (cgs_value, dimensions)

    raise UnitParseError("Cannot parse for unit data from '%s'. Please supply an expression of only Unit, Symbol, Pow, and Mul objects." % str(unit_expr))


def _lookup_unit_symbol(symbol_str, unit_symbol_lut):
    """
    Searches for the unit data typle corresponding to the given symbol.

    Parameters
    ----------
    symbol_str : str
        The unit symbol to look up.
    unit_symbol_lut : dict
        Dictionary with symbols as keys and unit data tuples as values.

    """

    if symbol_str in unit_symbol_lut:
        # lookup successful, return the tuple directly
        return default_unit_symbol_lut[symbol_str]

    # could still be a known symbol with a prefix
    possible_prefix = symbol_str[0]
    if possible_prefix in unit_prefixes:
        # the first character could be a prefix, check the rest of the symbol
        symbol_wo_prefix = symbol_str[1:]

        if symbol_wo_prefix in unit_symbol_lut:
            # lookup successful, it's a symbol with a prefix
            unit_data = unit_symbol_lut[symbol_wo_prefix]
            prefix_value = unit_prefixes[possible_prefix]

            # don't forget to account for the prefix value!
            return (unit_data[0] * prefix_value, unit_data[1])

    # no dice
    raise UnitParseError("Could not find unit symbol '%s' in the provided symbols. Please define this unit symbol." % symbol_str)
