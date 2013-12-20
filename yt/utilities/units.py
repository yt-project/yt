"""
Symbolic unit handling.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from sympy import Expr, Mul, Number, Pow, Rational, Symbol
from sympy import nsimplify, posify, sympify, latex
from sympy.parsing.sympy_parser import parse_expr
from collections import defaultdict
from yt.utilities.physical_ratios import \
    cm_per_pc, cm_per_ly, cm_per_au, cm_per_rsun, \
    mass_sun_grams, sec_per_year, sec_per_day, sec_per_hr, \
    sec_per_min, temp_sun_kelvin, luminosity_sun_ergs_per_sec, \
    metallicity_sun, erg_per_eV, amu_grams, mass_electron_grams, \
    hubble_constant_hertz
import string

# Define a sympy one object.
sympy_one = sympify(1)

#
# Exceptions
#

class SymbolNotFoundError(Exception):
    pass

class UnitParseError(Exception):
    pass

class InvalidUnitOperation(Exception):
    pass

#
# Base dimensions
#

mass = Symbol("(mass)", positive=True)
length = Symbol("(length)", positive=True)
time = Symbol("(time)", positive=True)
temperature = Symbol("(temperature)", positive=True)
metallicity = Symbol("(metallicity)", positive=True)

base_dimensions = [mass, length, time, temperature, metallicity]

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

    # "code" units, default to CGS conversion.
    # These default values are overridden in the code frontends
    "code_length" : (1.0, length),
    "unitary"   : (1.0, length),
    "code_mass" : (1.0, mass),
    "code_time" : (1.0, time),
    "code_velocity" : (1.0, velocity),
    "code_magnetic" : (1.0, magnetic_field),
    "code_temperature" : (1.0, temperature),
    "code_metallicity" : (1.0, metallicity),

    # other cgs
    "dyne": (1.0, force),
    "erg":  (1.0, energy),
    "esu":  (1.0, charge),

    # some SI
    "m": (1.0e2, length),
    "J": (1.0e7, energy),
    "Hz": (1.0, rate),

    # Imperial units
    "ft": (30.48, length),
    "mile": (160934, length),

    # dimensionless stuff
    "h": (1.0, dimensionless), # needs to be added for rho_crit_now

    # times
    "min": (sec_per_min, time),
    "hr":  (sec_per_hr, time),
    "day": (sec_per_day, time),
    "yr":  (sec_per_year, time),

    # Solar units
    "Msun": ( mass_sun_grams, mass),
    "msun": ( mass_sun_grams, mass),
    "Rsun": ( cm_per_rsun, length),
    "rsun": ( cm_per_rsun, length),
    "Lsun": ( luminosity_sun_ergs_per_sec, power),
    "Tsun": ( temp_sun_kelvin, temperature),
    "Zsun": ( metallicity_sun, metallicity),

    # astro distances
    "AU": (cm_per_au, length),
    "au": (cm_per_au, length),
    "ly": (cm_per_ly, length),
    "pc": (cm_per_pc, length),

    # electric stuff
    "gauss": (1.0, magnetic_field),

    # other astro
    "H_0": (hubble_constant_hertz, rate),  # check cf

    # misc
    "eV": (erg_per_eV, energy),
    "amu": (amu_grams, mass),
    "me": (mass_electron_grams, mass),

}

# Add LaTeX representations for units with trivial representations.
latex_symbol_lut = {
    "unitary" : "",
    "code_length" : "\\rm{code}\/\\rm{length}",
    "code_time" : "\\rm{code}\/\\rm{time}",
    "code_mass" : "\\rm{code}\/\\rm{mass}",
    "code_temperature" : "\\rm{code}\/\\rm{temperature}",
    "code_metallicity" : "\\rm{code}\/\\rm{metallicity}",
    "code_velocity" : "\\rm{code}\/\\rm{velocity}",
    "Msun" : "\\rm{M}_\\odot",
    "msun" : "\\rm{M}_\\odot",
    "Rsun" : "\\rm{R}_\\odot",
    "rsun" : "\\rm{R}_\\odot",
    "Lsun" : "\\rm{L}_\\odot",
    "Tsun" : "\\rm{T}_\\odot",
    "Zsun" : "\\rm{Z}_\\odot",
}
for key in default_unit_symbol_lut:
    if key not in latex_symbol_lut:
        latex_symbol_lut[key] = "\\rm{" + key + "}"

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
        self.unit_objs = {}

        if add_default_symbols:
            self.lut.update(default_unit_symbol_lut)

    def __getitem__(self, key):
        return self.lut[key]

    def add(self, symbol, cgs_value, dimensions, tex_repr=None):
        """
        Add a symbol to this registry.

        """
        # Validate
        if not isinstance(cgs_value, float):
            raise UnitParseError("cgs_value must be a float, got a %s." \
                                 % type(cgs_value))

        _validate_dimensions(dimensions)

        # Add to symbol lut
        if tex_repr is None:
            latex_symbol_lut[symbol] = "\\rm{" + symbol + "}"
        else:
            latex_symbol_lut[symbol] = tex_repr

        # Add to lut
        if tex_repr is None: tex_repr = symbol
        self.lut.update( {symbol: (cgs_value, dimensions)} )

    def remove(self, symbol):
        """
        Remove the entry for the unit matching `symbol`.

        """
        if symbol not in self.lut:
            raise SymbolNotFoundError(
                "Tried to remove the symbol '%s', but it does not exist" \
                "in this registry." % symbol)

        del self.lut[symbol]

    def modify(self, symbol, cgs_value):
        """
        Change the cgs value of a dimension.  Useful for adjusting code units
        after parsing parameters."

        """
        if symbol not in self.lut:
            raise SymbolNotFoundError(
                "Tried to remove the symbol '%s', but it does not exist" \
                "in this registry." % symbol)

        if hasattr(cgs_value, "in_cgs"):
            cgs_value = float(cgs_value.in_cgs().value)
        self.lut[symbol] = (cgs_value, self.lut[symbol][1])

    def keys(self):
        """
        Print out the units contained in the lookup table.

        """
        return self.lut.keys()

default_unit_registry = UnitRegistry()
SYMPIFY_ONE = sympify(1)

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

    def __new__(cls, unit_expr=SYMPIFY_ONE, cgs_value=None, dimensions=None,
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
        unit_key = None
        if isinstance(unit_expr, str):
            if registry and unit_expr in registry.unit_objs:
                return registry.unit_objs[unit_expr]
            elif unit_expr in default_unit_registry.unit_objs:
                return default_unit_registry.unit_objs[unit_expr]
            else:
                unit_key = unit_expr
                if not unit_expr:
                    # Bug catch...
                    # if unit_expr is an empty string, parse_expr fails hard...
                    unit_expr = "1"
                unit_expr = parse_expr(unit_expr)
        elif isinstance(unit_expr, Unit):
            # grab the unit object's sympy expression.
            unit_expr = unit_expr.expr
        # Make sure we have an Expr at this point.
        if not isinstance(unit_expr, Expr):
            raise UnitParseError("Unit representation must be a string or " \
                                 "sympy Expr. %s has type %s." \
                                 % (unit_expr, type(unit_expr)))

        if registry is None:
            # Caller did not set the registry, so use the default.
            registry = default_unit_registry

        # done with argument checking...

        # sympify, make positive symbols, and nsimplify the expr
        if unit_expr != SYMPIFY_ONE:
            unit_expr = sympify(unit_expr)
            unit_expr = _make_symbols_positive(unit_expr)
            if any([atom.is_Float for atom in unit_expr.atoms()]):
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
                raise UnitParseError("Could not use cgs_value as a float. " \
                                     "cgs_value is '%s' (type %s)." \
                                     % (cgs_value, type(cgs_value)) )

            # check that dimensions is valid
            _validate_dimensions( sympify(dimensions) )
        else:
            # lookup the unit symbols
            cgs_value, dimensions = _get_unit_data_from_expr(unit_expr, registry.lut)

        # Sympy trick to get dimensions powers as Rationals
        if not dimensions.is_Atom:
            if any([atom.is_Float for atom in dimensions.atoms()]):
                dimensions = nsimplify(dimensions)

        # Create obj with superclass construct.
        obj = Expr.__new__(cls, **assumptions)

        # Attach attributes to obj.
        obj.expr = unit_expr
        obj.is_atomic = is_atomic
        obj.cgs_value = cgs_value
        obj.dimensions = dimensions
        obj.registry = registry

        if unit_key:
            registry.unit_objs[unit_key] = obj

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
            return "dimensionless"
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
            raise InvalidUnitOperation("Tried to multiply a Unit object with " \
                                       "'%s' (type %s). This behavior is " \
                                       "undefined." % (u, type(u)) )

        return Unit(self.expr * u.expr,
                    cgs_value=(self.cgs_value * u.cgs_value),
                    dimensions=(self.dimensions * u.dimensions),
                    registry=self.registry)

    def __div__(self, u):
        """ Divide Unit by u (Unit object). """
        if not isinstance(u, Unit):
            raise InvalidUnitOperation("Tried to divide a Unit object by '%s' "\
                                       "(type %s). This behavior is "
                                       "undefined." % (u, type(u)) )

        return Unit(self.expr / u.expr,
                    cgs_value=(self.cgs_value / u.cgs_value),
                    dimensions=(self.dimensions / u.dimensions),
                    registry=self.registry)

    def __pow__(self, p):
        """ Take Unit to power p (float). """
        try:
            p = sympify(p)
        except ValueError:
            raise InvalidUnitOperation("Tried to take a Unit object to the " \
                                       "power '%s' (type %s). Failed to cast " \
                                       "it to a float." % (p, type(p)) )

        return Unit(self.expr**p, cgs_value=(self.cgs_value**p),
                    dimensions=(self.dimensions**p), registry=self.registry)

    def __eq__(self, u):
        """ Test unit equality. """
        if not isinstance(u, Unit):
            return False
        return \
          (self.cgs_value == u.cgs_value and self.dimensions == u.dimensions)

    def __ne__(self, u):
        """ Test unit inequality. """
        if not isinstance(u, Unit):
            return True
        return \
          (self.cgs_value != u.cgs_value or self.dimensions != u.dimensions)


    #
    # End unit operations
    #

    def same_dimensions_as(self, other_unit):
        """ Test if dimensions are the same. """
        return (self.dimensions / other_unit.dimensions) == sympy_one

    @property
    def is_dimensionless(self):
        return self.dimensions == sympy_one

    @property
    def is_code_unit(self):
        for atom in self.expr.atoms():
            if str(atom).startswith("code") or atom.is_Number:
                pass
            else:
                return False
        return True

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

    def latex_representation(self):
        symbol_table = {}
        for ex in self.expr.free_symbols:
            symbol_table[ex] = latex_symbol_lut[str(ex)]
        return latex(self.expr, symbol_names=symbol_table,
                     fold_frac_powers=True, fold_short_frac=True)
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
        raise InvalidUnitOperation(
            "Cannot convert from %s to %s because the dimensions do not "
            "match: %s and %s" % (old_units, new_units, old_units.dimensions,
                                  new_units.dimensions))

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
        if not s.is_positive:
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
            raise UnitParseError("dimensionality expression contains an "
                                 "unknown symbol '%s'." % d)

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
        conv = float(unit_data[0]**power)
        unit = unit_data[1]**power
        return (conv, unit)

    if isinstance(unit_expr, Mul):
        cgs_value = 1.0
        dimensions = 1
        for expr in unit_expr.args:
            unit_data = _get_unit_data_from_expr(expr, unit_symbol_lut)
            cgs_value *= unit_data[0]
            dimensions *= unit_data[1]

        return (float(cgs_value), dimensions)

    raise UnitParseError("Cannot parse for unit data from '%s'. Please supply" \
                         " an expression of only Unit, Symbol, Pow, and Mul" \
                         "objects." % str(unit_expr))


def _lookup_unit_symbol(symbol_str, unit_symbol_lut):
    """
    Searches for the unit data tuple corresponding to the given symbol.

    Parameters
    ----------
    symbol_str : str
        The unit symbol to look up.
    unit_symbol_lut : dict
        Dictionary with symbols as keys and unit data tuples as values.

    """
    if symbol_str in unit_symbol_lut:
        # lookup successful, return the tuple directly
        return unit_symbol_lut[symbol_str]

    # could still be a known symbol with a prefix
    possible_prefix = symbol_str[0]
    if possible_prefix in unit_prefixes:
        # the first character could be a prefix, check the rest of the symbol
        symbol_wo_prefix = symbol_str[1:]

        if symbol_wo_prefix in unit_symbol_lut:
            # lookup successful, it's a symbol with a prefix
            unit_data = unit_symbol_lut[symbol_wo_prefix]
            prefix_value = unit_prefixes[possible_prefix]

            if symbol_str not in latex_symbol_lut:
                latex_symbol_lut[symbol_str] = \
                    string.replace(latex_symbol_lut[symbol_wo_prefix],
                                   '{'+symbol_wo_prefix+'}', '{'+symbol_str+'}')

            # don't forget to account for the prefix value!
            return (unit_data[0] * prefix_value, unit_data[1])

    # no dice
    raise UnitParseError("Could not find unit symbol '%s' in the provided " \
                         "symbols." % symbol_str)
