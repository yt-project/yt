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

from yt.extern.six import text_type
from sympy import \
    Expr, Mul, Add, Number, \
    Pow, Symbol, Integer, \
    Float, Basic, Rational, sqrt
from sympy.core.numbers import One
from sympy import sympify, latex
from sympy.parsing.sympy_parser import \
    parse_expr, auto_number, rationalize
from keyword import iskeyword
from yt.units.dimensions import \
    base_dimensions, temperature, \
    dimensionless, current_mks, \
    angle
from yt.units.equivalencies import \
    equivalence_registry
from yt.units.unit_lookup_table import \
    unit_prefixes, prefixable_units, latex_prefixes, \
    default_base_units
from yt.units.unit_registry import \
    UnitRegistry, \
    UnitParseError
from yt.utilities.exceptions import YTUnitsNotReducible

import copy
import token

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

unit_system_registry = {}

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

def get_latex_representation(expr, registry):
    symbol_table = {}
    for ex in expr.free_symbols:
        try:
            symbol_table[ex] = registry.lut[str(ex)][3]
        except:
            symbol_table[ex] = r"\rm{" + str(ex).replace('_', '\ ') + "}"

    # invert the symbol table dict to look for keys with identical values
    invert_symbols = {}
    for key, value in symbol_table.items():
        if value not in invert_symbols:
            invert_symbols[value] = [key]
        else:
            invert_symbols[value].append(key)

    # if there are any units with identical latex representations, substitute
    # units to avoid  uncanceled terms in the final latex expresion.
    for val in invert_symbols:
        symbols = invert_symbols[val]
        for i in range(1, len(symbols)):
            expr = expr.subs(symbols[i], symbols[0])
    prefix = None
    if isinstance(expr, Mul):
        coeffs = expr.as_coeff_Mul()
        if coeffs[0] == 1 or not isinstance(coeffs[0], Float):
            pass
        else:
            expr = coeffs[1]
            prefix = Float(coeffs[0], 2)
    latex_repr = latex(expr, symbol_names=symbol_table, mul_symbol="dot",
                       fold_frac_powers=True, fold_short_frac=True)

    if prefix is not None:
        latex_repr = latex(prefix, mul_symbol="times") + '\\ ' + latex_repr

    if latex_repr == '1':
        return ''
    else:
        return latex_repr

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
    __slots__ = ["expr", "is_atomic", "base_value", "base_offset", "dimensions",
                 "registry", "_latex_repr"]

    def __new__(cls, unit_expr=sympy_one, base_value=None, base_offset=0.0,
                dimensions=None, registry=None, latex_repr=None, **assumptions):
        """
        Create a new unit. May be an atomic unit (like a gram) or combinations
        of atomic units (like g / cm**3).

        Parameters
        ----------
        unit_expr : Unit object, sympy.core.expr.Expr object, or str
            The symbolic unit expression.
        base_value : float
            The unit's value in yt's base units.
        base_offset : float
            The offset necessary to normalize temperature units to a common
            zero point.
        dimensions : sympy.core.expr.Expr
            A sympy expression representing the dimensionality of this unit.
            It must contain only mass, length, time, temperature and angle
            symbols.
        registry : UnitRegistry object
            The unit registry we use to interpret unit symbols.
        latex_repr : string
            A string to render the unit as LaTeX

        Additional keyword arguments are passed as assumptions to the Sympy Expr
        initializer

        """
        # Simplest case. If user passes a Unit object, just use the expr.
        unit_key = None
        if isinstance(unit_expr, (str, bytes, text_type)):
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
                try:
                    unit_expr = parse_expr(unit_expr, global_dict=global_dict,
                                           transformations=unit_text_transform)
                except SyntaxError as e:
                    msg = ("Unit expression %s raised an error "
                           "during parsing:\n%s" % (unit_expr, repr(e)))
                    raise UnitParseError(msg)
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
        # check base_value and dimensions
        #

        if base_value is not None:
            # check that base_value is a float or can be converted to one
            try:
                base_value = float(base_value)
            except ValueError:
                raise UnitParseError("Could not use base_value as a float. " \
                                     "base_value is '%s' (type %s)." \
                                     % (base_value, type(base_value)) )

            # check that dimensions is valid
            if dimensions is not None:
                validate_dimensions(dimensions)
        else:
            # lookup the unit symbols
            unit_data = _get_unit_data_from_expr(unit_expr, registry.lut)
            base_value = unit_data[0]
            dimensions = unit_data[1]
            if len(unit_data) > 2:
                base_offset = unit_data[2]
                latex_repr = unit_data[3]
            else:
                base_offset = 0.0

        # Create obj with superclass construct.
        obj = Expr.__new__(cls, **assumptions)

        # Attach attributes to obj.
        obj.expr = unit_expr
        obj.is_atomic = is_atomic
        obj.base_value = base_value
        obj.base_offset = base_offset
        obj.dimensions = dimensions
        obj._latex_repr = latex_repr
        obj.registry = registry

        if unit_key is not None:
            registry.unit_objs[unit_key] = obj

        # Return `obj` so __init__ can handle it.

        return obj

    _latex_expr = None
    @property
    def latex_repr(self):
        if self._latex_repr is not None:
            return self._latex_repr
        if self.expr.is_Atom:
            expr = self.expr
        else:
            expr = self.expr.copy()
        self._latex_repr = get_latex_representation(expr, self.registry)
        return self._latex_repr

    ### Some sympy conventions
    def __getnewargs__(self):
        return (self.expr, self.is_atomic, self.base_value, self.dimensions,
                self.registry)

    def __hash__(self):
        return super(Unit, self).__hash__()

    def _hashable_content(self):
        return (self.expr, self.is_atomic, self.base_value, self.dimensions,
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

        base_offset = 0.0
        if self.base_offset or u.base_offset:
            if u.dimensions in (temperature, angle) and self.is_dimensionless:
                base_offset = u.base_offset
            elif self.dimensions in (temperature, angle) and u.is_dimensionless:
                base_offset = self.base_offset
            else:
                raise InvalidUnitOperation("Quantities with units of Fahrenheit "
                                           "and Celsius or angles cannot be multiplied.")

        return Unit(self.expr * u.expr,
                    base_value=(self.base_value * u.base_value),
                    base_offset=base_offset,
                    dimensions=(self.dimensions * u.dimensions),
                    registry=self.registry)

    def __div__(self, u):
        """ Divide Unit by u (Unit object). """
        if not isinstance(u, Unit):
            raise InvalidUnitOperation("Tried to divide a Unit object by '%s' "
                                       "(type %s). This behavior is "
                                       "undefined." % (u, type(u)))

        base_offset = 0.0
        if self.base_offset or u.base_offset:
            if u.dimensions in (temperature, angle) and self.is_dimensionless:
                base_offset = u.base_offset
            elif self.dimensions in (temperature, angle) and u.is_dimensionless:
                base_offset = self.base_offset
            else:
                raise InvalidUnitOperation("Quantities with units of Farhenheit "
                                           "and Celsius cannot be multiplied.")

        return Unit(self.expr / u.expr,
                    base_value=(self.base_value / u.base_value),
                    base_offset=base_offset,
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

        return Unit(self.expr**p, base_value=(self.base_value**p),
                    dimensions=(self.dimensions**p),
                    registry=self.registry)

    def __eq__(self, u):
        """ Test unit equality. """
        if not isinstance(u, Unit):
            return False
        return \
          (self.base_value == u.base_value and self.dimensions == u.dimensions)

    def __ne__(self, u):
        """ Test unit inequality. """
        if not isinstance(u, Unit):
            return True
        if self.base_value != u.base_value:
            return True
        # use 'is' comparison dimensions to avoid expensive sympy operation
        if self.dimensions is u.dimensions:
            return False
        return self.dimensions != u.dimensions

    def copy(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memodict=None):
        if memodict is None:
            memodict = {}
        expr = str(self.expr)
        base_value = copy.deepcopy(self.base_value)
        base_offset = copy.deepcopy(self.base_offset)
        dimensions = copy.deepcopy(self.dimensions)
        lut = copy.deepcopy(self.registry.lut)
        registry = UnitRegistry(lut=lut)
        return Unit(expr, base_value, base_offset, dimensions, registry)

    #
    # End unit operations
    #

    def same_dimensions_as(self, other_unit):
        """ Test if dimensions are the same. """
        # test first for 'is' equality to avoid expensive sympy operation
        if self.dimensions is other_unit.dimensions:
            return True
        return (self.dimensions / other_unit.dimensions) == sympy_one

    @property
    def is_dimensionless(self):
        return self.dimensions is sympy_one

    @property
    def is_code_unit(self):
        for atom in self.expr.atoms():
            if str(atom).startswith("code") or atom.is_Number:
                pass
            else:
                return False
        return True

    def list_equivalencies(self):
        """
        Lists the possible equivalencies associated with this unit object
        """
        for k, v in equivalence_registry.items():
            if self.has_equivalent(k):
                print(v())

    def has_equivalent(self, equiv):
        """
        Check to see if this unit object as an equivalent unit in *equiv*.
        """
        try:
            this_equiv = equivalence_registry[equiv]()
        except KeyError:
            raise KeyError("No such equivalence \"%s\"." % equiv)
        old_dims = self.dimensions
        return old_dims in this_equiv.dims

    def get_base_equivalent(self, unit_system="cgs"):
        """
        Create and return dimensionally-equivalent units in a specified base.
        """
        yt_base_unit_string = _get_system_unit_string(self.dimensions, default_base_units)
        yt_base_unit = Unit(yt_base_unit_string, base_value=1.0,
                            dimensions=self.dimensions, registry=self.registry)
        if unit_system == "cgs":
            if current_mks in self.dimensions.free_symbols:
                raise YTUnitsNotReducible(self, "cgs")
            return yt_base_unit
        else:
            if unit_system == "code":
                raise RuntimeError(r'You must refer to a dataset instance to convert to a '
                                   r'code unit system. Try again with unit_system=ds instead, '
                                   r'where \'ds\' is your dataset.')
            unit_system = unit_system_registry[str(unit_system)]
            units_string = _get_system_unit_string(self.dimensions, unit_system.base_units)
            u = Unit(units_string, registry=self.registry)
            base_value = get_conversion_factor(self, yt_base_unit)[0]
            base_value /= get_conversion_factor(self, u)[0]
            return Unit(units_string, base_value=base_value,
                        dimensions=self.dimensions, registry=self.registry)

    def get_cgs_equivalent(self):
        """
        Create and return dimensionally-equivalent cgs units.
        """
        return self.get_base_equivalent(unit_system="cgs")

    def get_mks_equivalent(self):
        """
        Create and return dimensionally-equivalent mks units.
        """
        return self.get_base_equivalent(unit_system="mks")

    def get_conversion_factor(self, other_units):
        return get_conversion_factor(self, other_units)

    def latex_representation(self):
        return self.latex_repr

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
    ratio = old_units.base_value / new_units.base_value
    if old_units.base_offset == 0 and new_units.base_offset == 0:
        return (ratio, None)
    else:
        if old_units.dimensions in (temperature, angle):
            return ratio, ratio*old_units.base_offset - new_units.base_offset
        else:
            raise InvalidUnitOperation(
                "Fahrenheit and Celsius are not absolute temperature scales "
                "and cannot be used in compound unit symbols.")

#
# Helper functions
#

def _get_unit_data_from_expr(unit_expr, unit_symbol_lut):
    """
    Grabs the total base_value and dimensions from a valid unit expression.

    Parameters
    ----------
    unit_expr: Unit object, or sympy Expr object
        The expression containing unit symbols.
    unit_symbol_lut: dict
        Provides the unit data for each valid unit symbol.

    """
    # The simplest case first
    if isinstance(unit_expr, Unit):
        return (unit_expr.base_value, unit_expr.dimensions)

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
        base_value = 1.0
        dimensions = 1
        for expr in unit_expr.args:
            unit_data = _get_unit_data_from_expr(expr, unit_symbol_lut)
            base_value *= unit_data[0]
            dimensions *= unit_data[1]

        return (float(base_value), dimensions)

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

        unit_is_si_prefixable = (symbol_wo_prefix in unit_symbol_lut and
                                 symbol_wo_prefix in prefixable_units)

        if unit_is_si_prefixable is True:
            # lookup successful, it's a symbol with a prefix
            unit_data = unit_symbol_lut[symbol_wo_prefix]
            prefix_value = unit_prefixes[possible_prefix]

            if possible_prefix in latex_prefixes:
                latex_repr = symbol_str.replace(
                    possible_prefix, '{'+latex_prefixes[possible_prefix]+'}')
            else:
                # Need to add some special handling for comoving units
                # this is fine for now, but it wouldn't work for a general
                # unit that has an arbitrary LaTeX representation
                if symbol_wo_prefix != 'cm' and symbol_wo_prefix.endswith('cm'):
                    sub_symbol_wo_prefix = symbol_wo_prefix[:-2]
                    sub_symbol_str = symbol_str[:-2]
                else:
                    sub_symbol_wo_prefix = symbol_wo_prefix
                    sub_symbol_str = symbol_str
                latex_repr = unit_data[3].replace(
                    '{' + sub_symbol_wo_prefix + '}', '{' + sub_symbol_str + '}')

            # Leave offset and dimensions the same, but adjust scale factor and
            # LaTeX representation
            ret = (unit_data[0] * prefix_value, unit_data[1], unit_data[2],
                   latex_repr)

            unit_symbol_lut[symbol_str] = ret

            return ret

    # no dice
    if symbol_str.startswith('code_'):
        raise UnitParseError(
            "Code units have not been defined. \n"
            "Try creating the array or quantity using ds.arr or ds.quan instead.")
    else:
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

def _get_system_unit_string(dimensions, base_units):
    # The dimensions of a unit object is the product of the base dimensions.
    # Use sympy to factor the dimensions into base CGS unit symbols.
    units = []
    my_dims = dimensions.expand()
    if my_dims is dimensionless:
        return ""
    if my_dims in base_units:
        return base_units[my_dims]
    for factor in my_dims.as_ordered_factors():
        dim = list(factor.free_symbols)[0]
        unit_string = str(base_units[dim])
        if factor.is_Pow:
            power_string = "**(%s)" % factor.as_base_exp()[1]
        else:
            power_string = ""
        units.append("(%s)%s" % (unit_string, power_string))
    return " * ".join(units)
