"""
A class that represents a unit symbol.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from sympy import \
    Expr, Mul, Add, Number, \
    Pow, Symbol, Integer, \
    Float, Basic, Rational, sqrt
from sympy.core.numbers import One
from sympy import sympify, latex, symbols
from sympy.parsing.sympy_parser import \
    parse_expr, auto_number, rationalize
from keyword import iskeyword
from yt.units.dimensions import \
    base_dimensions, temperature, \
    dimensionless
from yt.units.unit_lookup_table import \
    latex_symbol_lut, unit_prefixes, \
    prefixable_units, cgs_base_units, \
    mks_base_units, latex_prefixes, \
    cgs_conversions
from yt.units.unit_registry import UnitRegistry

import copy
import string
import token

class UnitParseError(Exception):
    pass

class InvalidUnitOperation(Exception):
    pass

default_unit_registry = UnitRegistry()

sympy_one = sympify(1)

global_dict = {
    'Symbol': Symbol,
    'Integer': Integer,
    'Float': Float,
    'Rational': Rational,
    'sqrt': sqrt
}

def auto_positive_symbol(tokens, local_dict, global_dict):
    """
    Inserts calls to ``Symbol`` for undefined variables.
    Passes in positive=True as a keyword argument.
    Adapted from sympy.sympy.parsing.sympy_parser.auto_symbol
    """
    result = []
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == token.NAME:
            name = tokVal

            if (name in ['True', 'False', 'None']
                or iskeyword(name)
                or name in local_dict
                # Don't convert attribute access
                or (prevTok[0] == token.OP and prevTok[1] == '.')
                # Don't convert keyword arguments
                or (prevTok[0] == token.OP and prevTok[1] in ('(', ',')
                    and nextTokNum == token.OP and nextTokVal == '=')):
                result.append((token.NAME, name))
                continue
            elif name in global_dict:
                obj = global_dict[name]
                if isinstance(obj, (Basic, type)) or callable(obj):
                    result.append((token.NAME, name))
                    continue

            result.extend([
                (token.NAME, 'Symbol'),
                (token.OP, '('),
                (token.NAME, repr(str(name))),
                (token.OP, ','),
                (token.NAME, 'positive'),
                (token.OP, '='),
                (token.NAME, 'True'),
                (token.OP, ')'),
            ])
        else:
            result.append((tokNum, tokVal))

        prevTok = (tokNum, tokVal)

    return result

unit_text_transform = (auto_positive_symbol, rationalize, auto_number)

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
    __slots__ = ["expr", "is_atomic", "cgs_value", "cgs_offset", "dimensions",
                 "registry", "cgs_conversion", "is_mks"]

    def __new__(cls, unit_expr=sympy_one, cgs_value=None, cgs_offset=0.0,
                dimensions=None, registry=None, **assumptions):
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
            It must contain only mass, length, time, temperature and angle
            symbols.
        offset : float
            The offset necessary to normalize temperature units to a common
            zero point.
        registry : UnitRegistry object
            The unit registry we use to interpret unit symbols.

        """
        # Simplest case. If user passes a Unit object, just use the expr.
        unit_key = None
        if isinstance(unit_expr, (str, bytes, unicode)):
            if isinstance(unit_expr, bytes):
                unit_expr = unit_expr.decode("utf-8")

            if registry and unit_expr in registry.unit_objs:
                return registry.unit_objs[unit_expr]
            else:
                unit_key = unit_expr
                if not unit_expr:
                    # Bug catch...
                    # if unit_expr is an empty string, parse_expr fails hard...
                    unit_expr = "1"
                unit_expr = parse_expr(unit_expr, global_dict=global_dict,
                                       transformations=unit_text_transform)
        elif isinstance(unit_expr, Unit):
            # grab the unit object's sympy expression.
            unit_expr = unit_expr.expr
        # Make sure we have an Expr at this point.
        if not isinstance(unit_expr, Expr):
            raise UnitParseError("Unit representation must be a string or " \
                                 "sympy Expr. %s has type %s." \
                                 % (unit_expr, type(unit_expr)))

        if unit_expr == sympy_one and dimensions is None:
            dimensions = dimensionless

        if registry is None:
            # Caller did not set the registry, so use the default.
            registry = default_unit_registry

        # done with argument checking...

        # see if the unit is atomic.
        is_atomic = False
        if isinstance(unit_expr, Symbol):
            is_atomic = True

        #
        # check cgs_value and dimensions
        #

        if cgs_value is not None:
            # check that cgs_value is a float or can be converted to one
            try:
                cgs_value = float(cgs_value)
            except ValueError:
                raise UnitParseError("Could not use cgs_value as a float. " \
                                     "cgs_value is '%s' (type %s)." \
                                     % (cgs_value, type(cgs_value)) )

            # check that dimensions is valid
            if dimensions is not None:
                validate_dimensions(dimensions)
        else:
            # lookup the unit symbols
            unit_data = _get_unit_data_from_expr(unit_expr, registry.lut)
            cgs_value = unit_data[0]
            dimensions = unit_data[1]
            if len(unit_data) == 3:
                cgs_offset = unit_data[2]

        # Create obj with superclass construct.
        obj = Expr.__new__(cls, **assumptions)

        # Attach attributes to obj.
        obj.expr = unit_expr
        obj.is_atomic = is_atomic
        obj.cgs_value = cgs_value
        obj.cgs_offset = cgs_offset
        obj.dimensions = dimensions
        obj.registry = registry

        check_atoms = [atom for atom in unit_expr.free_symbols
                       if str(atom) in cgs_conversions]
        if len(check_atoms) > 0:
            conversions = []
            for atom in check_atoms:
                conversions.append((atom,symbols(cgs_conversions[str(atom)])))
            conversion = Unit(unit_expr=unit_expr.subs(conversions),
                              registry=registry)
            is_mks = True
        else:
            conversion = None
            is_mks = False
        obj.cgs_conversion = conversion
        obj.is_mks = is_mks

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
            raise InvalidUnitOperation("Tried to multiply a Unit object with "
                                       "'%s' (type %s). This behavior is "
                                       "undefined." % (u, type(u)))

        cgs_offset = 0.0
        if self.cgs_offset or u.cgs_offset:
            if u.dimensions is temperature and self.is_dimensionless:
                cgs_offset = u.cgs_offset
            elif self.dimensions is temperature and u.is_dimensionless:
                cgs_offset = self.cgs_offset
            else:
                raise InvalidUnitOperation("Quantities with units of Fahrenheit "
                                           "and Celsius cannot be multiplied.")

        return Unit(self.expr * u.expr,
                    cgs_value=(self.cgs_value * u.cgs_value),
                    cgs_offset=cgs_offset,
                    dimensions=(self.dimensions * u.dimensions),
                    registry=self.registry)

    def __div__(self, u):
        """ Divide Unit by u (Unit object). """
        if not isinstance(u, Unit):
            raise InvalidUnitOperation("Tried to divide a Unit object by '%s' "
                                       "(type %s). This behavior is "
                                       "undefined." % (u, type(u)))

        cgs_offset = 0.0
        if self.cgs_offset or u.cgs_offset:
            if u.dimensions is dims.temperature and self.is_dimensionless:
                cgs_offset = u.cgs_offset
            elif self.dimensions is dims.temperature and u.is_dimensionless:
                cgs_offset = self.cgs_offset
            else:
                raise InvalidUnitOperation("Quantities with units of Farhenheit "
                                           "and Celsius cannot be multiplied.")

        return Unit(self.expr / u.expr,
                    cgs_value=(self.cgs_value / u.cgs_value),
                    cgs_offset=cgs_offset,
                    dimensions=(self.dimensions / u.dimensions),
                    registry=self.registry)

    __truediv__ = __div__

    def __pow__(self, p):
        """ Take Unit to power p (float). """
        try:
            p = Rational(str(p)).limit_denominator()
        except ValueError:
            raise InvalidUnitOperation("Tried to take a Unit object to the " \
                                       "power '%s' (type %s). Failed to cast " \
                                       "it to a float." % (p, type(p)) )

        return Unit(self.expr**p, cgs_value=(self.cgs_value**p),
                    dimensions=(self.dimensions**p),
                    registry=self.registry)

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

    def copy(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memodict=None):
        if memodict is None:
            memodict = {}
        expr = str(self.expr)
        cgs_value = copy.deepcopy(self.cgs_value)
        cgs_offset = copy.deepcopy(self.cgs_offset)
        dimensions = copy.deepcopy(self.dimensions)
        lut = copy.deepcopy(self.registry.lut)
        registry = UnitRegistry(lut=lut)
        return Unit(expr, cgs_value, cgs_offset, dimensions, registry)

    #
    # End unit operations
    #

    def same_dimensions_as(self, other_unit):
        """ Test if dimensions are the same. """
        first_check = False
        second_check = False
        if self.cgs_conversion:
            first_check = self.cgs_conversion.dimensions / other_unit.dimensions == sympy_one
        if other_unit.cgs_conversion:
            second_check = other_unit.cgs_conversion.dimensions / self.dimensions == sympy_one
        if first_check or second_check:
            return True
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

    def _get_system_unit_string(self, base_units):
        # The dimensions of a unit object is the product of the base dimensions.
        # Use sympy to factor the dimensions into base CGS unit symbols.
        units = []
        my_dims = self.dimensions.expand()
        for dim in base_units:
            unit_string = base_units[dim]
            power_string = "**(%s)" % my_dims.as_coeff_exponent(dim)[1]
            units.append("".join([unit_string, power_string]))
        return " * ".join(units)


    def get_cgs_equivalent(self):
        """
        Create and return dimensionally-equivalent cgs units.

        """
        if self.cgs_conversion:
            units = self.cgs_conversion
        else:
            units = self
        units_string = units._get_system_unit_string(cgs_base_units)
        return Unit(units_string, cgs_value=1.0,
                    dimensions=units.dimensions, registry=self.registry)

    def get_mks_equivalent(self):
        """
        Create and return dimensionally-equivalent mks units.

        """
        if self.cgs_conversion and not self.is_mks:
            units = self.cgs_conversion
        else:
            units = self
        units_string = units._get_system_unit_string(mks_base_units)
        cgs_value = (get_conversion_factor(units, units.get_cgs_equivalent())[0] /
                     get_conversion_factor(units, Unit(units_string))[0])
        return Unit(units_string, cgs_value=cgs_value,
                    dimensions=units.dimensions, registry=self.registry)

    def get_conversion_factor(self, other_units):
        return get_conversion_factor(self, other_units)

    def latex_representation(self):
        symbol_table = {}
        for ex in self.expr.free_symbols:
            symbol_table[ex] = latex_symbol_lut[str(ex)]
        return latex(self.expr, symbol_names=symbol_table,
                     mul_symbol="dot", fold_frac_powers=True,
                     fold_short_frac=True)
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
    offset : float or None
        Offset between the old unit and new unit.

    """
    ratio = old_units.cgs_value / new_units.cgs_value
    if old_units.cgs_offset == 0 and new_units.cgs_offset == 0:
        return (ratio, None)
    else:
        if old_units.dimensions is temperature:
            return ratio, ratio*old_units.cgs_offset - new_units.cgs_offset
        else:
            raise InvalidUnitOperation(
                "Fahrenheit and Celsius are not absolute temperature scales "
                "and cannot be used in compound unit symbols.")

#
# Helper functions
#

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
        return (float(unit_expr), sympy_one)

    if isinstance(unit_expr, Pow):
        unit_data = _get_unit_data_from_expr(unit_expr.args[0], unit_symbol_lut)
        power = unit_expr.args[1]
        if isinstance(power, Symbol):
            raise UnitParseError("Invalid unit expression '%s'." % unit_expr)
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

        if symbol_wo_prefix in unit_symbol_lut and symbol_wo_prefix in prefixable_units:
            # lookup successful, it's a symbol with a prefix
            unit_data = unit_symbol_lut[symbol_wo_prefix]
            prefix_value = unit_prefixes[possible_prefix]

            if symbol_str not in latex_symbol_lut:
                if possible_prefix in latex_prefixes:
                    sstr = symbol_str.replace(possible_prefix,
                                              '{'+latex_prefixes[possible_prefix]+'}')
                else:
                    sstr = symbol_str
                latex_symbol_lut[symbol_str] = \
                    latex_symbol_lut[symbol_wo_prefix].replace(
                                   '{'+symbol_wo_prefix+'}', '{'+sstr+'}')

            # don't forget to account for the prefix value!
            return (unit_data[0] * prefix_value, unit_data[1])

    # no dice
    raise UnitParseError("Could not find unit symbol '%s' in the provided " \
                         "symbols." % symbol_str)

def validate_dimensions(dimensions):
    if isinstance(dimensions, Mul):
        for dim in dimensions.args:
            validate_dimensions(dim)
    elif isinstance(dimensions, Symbol):
        if dimensions not in base_dimensions:
            raise UnitParseError("Dimensionality expression contains an "
                                 "unknown symbol '%s'." % dimensions)
    elif isinstance(dimensions, Pow):
        if not isinstance(dimensions.args[1], Number):
            raise UnitParseError("Dimensionality expression '%s' contains a "
                                 "unit symbol as a power." % dimensions)
    elif isinstance(dimensions, (Add, Number)):
        if not isinstance(dimensions, One):
            raise UnitParseError("Only dimensions that are instances of Pow, "
                                 "Mul, or symbols in the base dimensions are "
                                 "allowed.  Got dimensions '%s'" % dimensions)
    elif not isinstance(dimensions, Basic):
        raise UnitParseError("Bad dimensionality expression '%s'." % dimensions)
