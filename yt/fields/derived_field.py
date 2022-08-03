import contextlib
import inspect
import re
from typing import Optional, Tuple, Union

from more_itertools import always_iterable

import yt.units.dimensions as ytdims
from yt._maintenance.deprecation import issue_deprecation_warning
from yt.funcs import iter_fields, validate_field_key
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.logger import ytLogger as mylog

from .field_detector import FieldDetector
from .field_exceptions import (
    FieldUnitsError,
    NeedsDataField,
    NeedsGridType,
    NeedsOriginalGrid,
    NeedsParameter,
    NeedsProperty,
)


def TranslationFunc(field_name):
    def _TranslationFunc(field, data):
        # We do a bunch of in-place modifications, so we will copy this.
        return data[field_name].copy()

    _TranslationFunc.alias_name = field_name
    return _TranslationFunc


def NullFunc(field, data):
    raise YTFieldNotFound(field.name)


def DeprecatedFieldFunc(ret_field, func, since, removal):
    def _DeprecatedFieldFunc(field, data):
        # Only log a warning if we've already done
        # field detection
        if data.ds.fields_detected:
            args = [field.name, since]
            msg = "The Derived Field %s is deprecated as of yt v%s "
            if removal is not None:
                msg += "and will be removed in yt v%s "
                args.append(removal)
            if ret_field != field.name:
                msg += ", use %s instead"
                args.append(ret_field)
            mylog.warning(msg, *args)
        return func(field, data)

    return _DeprecatedFieldFunc


class DerivedField:
    """
    This is the base class used to describe a cell-by-cell derived field.

    Parameters
    ----------

    name : str
       is the name of the field.
    function : callable
       A function handle that defines the field.  Should accept
       arguments (field, data)
    units : str
       A plain text string encoding the unit, or a query to a unit system of
       a dataset. Powers must be in Python syntax (** instead of ^). If set
       to 'auto' or None (default), units will be inferred from the return value
       of the field function.
    take_log : bool
       Describes whether the field should be logged
    validators : list
       A list of :class:`FieldValidator` objects
    sampling_type : string, default = "cell"
        How is the field sampled?  This can be one of the following options at
        present: "cell" (cell-centered), "discrete" (or "particle") for
        discretely sampled data.
    vector_field : bool
       Describes the dimensionality of the field.  Currently unused.
    display_field : bool
       Governs its appearance in the dropdowns in Reason
    not_in_all : bool
       Used for baryon fields from the data that are not in all the grids
    display_name : str
       A name used in the plots
    output_units : str
       For fields that exist on disk, which we may want to convert to other
       fields or that get aliased to themselves, we can specify a different
       desired output unit than the unit found on disk.
    dimensions : str or object from yt.units.dimensions
       The dimensions of the field, only used for error checking with units='auto'.
    nodal_flag : array-like with three components
       This describes how the field is centered within a cell. If nodal_flag
       is [0, 0, 0], then the field is cell-centered. If any of the components
       of nodal_flag are 1, then the field is nodal in that direction, meaning
       it is defined at the lo and hi sides of the cell rather than at the center.
       For example, a field with nodal_flag = [1, 0, 0] would be defined at the
       middle of the 2 x-faces of each cell. nodal_flag = [0, 1, 1] would mean the
       that the field defined at the centers of the 4 edges that are normal to the
       x axis, while nodal_flag = [1, 1, 1] would be defined at the 8 cell corners.
    """

    _inherited_particle_filter = False

    def __init__(
        self,
        name: Tuple[str, str],
        sampling_type,
        function,
        units: Optional[Union[str, bytes, Unit]] = None,
        take_log=True,
        validators=None,
        vector_field=False,
        display_field=True,
        not_in_all=False,
        display_name=None,
        output_units=None,
        dimensions=None,
        ds=None,
        nodal_flag=None,
        *,
        alias: Optional["DerivedField"] = None,
        particle_type=None,
    ):
        validate_field_key(name)
        self.name = name
        self.take_log = take_log
        self.display_name = display_name
        self.not_in_all = not_in_all
        self.display_field = display_field
        if particle_type is not None:
            issue_deprecation_warning(
                "The 'particle_type' keyword argument is deprecated. "
                "Please use sampling_type='particle' instead.",
                since="3.2",
                removal="4.2",
            )
            sampling_type = "particle"
        self.sampling_type = sampling_type
        self.vector_field = vector_field
        self.ds = ds

        if self.ds is not None:
            self._ionization_label_format = self.ds._ionization_label_format
        else:
            self._ionization_label_format = "roman_numeral"

        if nodal_flag is None:
            self.nodal_flag = [0, 0, 0]
        else:
            self.nodal_flag = nodal_flag

        self._function = function

        self.validators = list(always_iterable(validators))

        # handle units
        self.units: Optional[Union[str, bytes, Unit]]
        if units in (None, "auto"):
            self.units = None
        elif isinstance(units, str):
            self.units = units
        elif isinstance(units, Unit):
            self.units = str(units)
        elif isinstance(units, bytes):
            self.units = units.decode("utf-8")
        else:
            raise FieldUnitsError(
                "Cannot handle units '%s' (type %s). "
                "Please provide a string or Unit "
                "object." % (units, type(units))
            )
        if output_units is None:
            output_units = self.units
        self.output_units = output_units

        if isinstance(dimensions, str):
            dimensions = getattr(ytdims, dimensions)
        self.dimensions = dimensions

        if alias is None:
            self._shared_aliases_list = [self]
        else:
            self._shared_aliases_list = alias._shared_aliases_list
            self._shared_aliases_list.append(self)

    def _copy_def(self):
        dd = {}
        dd["name"] = self.name
        dd["units"] = self.units
        dd["take_log"] = self.take_log
        dd["validators"] = list(self.validators)
        dd["sampling_type"] = self.sampling_type
        dd["vector_field"] = self.vector_field
        dd["display_field"] = True
        dd["not_in_all"] = self.not_in_all
        dd["display_name"] = self.display_name
        return dd

    @property
    def is_sph_field(self):
        if self.sampling_type == "cell":
            return False
        is_sph_field = False
        if self.is_alias:
            name = self.alias_name
        else:
            name = self.name
        if hasattr(self.ds, "_sph_ptypes"):
            is_sph_field |= name[0] in (self.ds._sph_ptypes + ("gas",))
        return is_sph_field

    @property
    def local_sampling(self):
        return self.sampling_type in ("discrete", "particle", "local")

    def get_units(self):
        if self.ds is not None:
            u = Unit(self.units, registry=self.ds.unit_registry)
        else:
            u = Unit(self.units)
        return u.latex_representation()

    def get_projected_units(self):
        if self.ds is not None:
            u = Unit(self.units, registry=self.ds.unit_registry)
        else:
            u = Unit(self.units)
        return (u * Unit("cm")).latex_representation()

    def check_available(self, data):
        """
        This raises an exception of the appropriate type if the set of
        validation mechanisms are not met, and otherwise returns True.
        """
        for validator in self.validators:
            validator(data)
        # If we don't get an exception, we're good to go
        return True

    def get_dependencies(self, *args, **kwargs):
        """
        This returns a list of names of fields that this field depends on.
        """
        e = FieldDetector(*args, **kwargs)
        e[self.name]
        return e

    def _get_needed_parameters(self, fd):
        params = []
        values = []
        permute_params = {}
        vals = [v for v in self.validators if isinstance(v, ValidateParameter)]
        for val in vals:
            if val.parameter_values is not None:
                permute_params.update(val.parameter_values)
            else:
                params.extend(val.parameters)
                values.extend([fd.get_field_parameter(fp) for fp in val.parameters])
        return dict(zip(params, values)), permute_params

    _unit_registry = None

    @contextlib.contextmanager
    def unit_registry(self, data):
        old_registry = self._unit_registry
        if hasattr(data, "unit_registry"):
            ur = data.unit_registry
        elif hasattr(data, "ds"):
            ur = data.ds.unit_registry
        else:
            ur = None
        self._unit_registry = ur
        yield
        self._unit_registry = old_registry

    def __call__(self, data):
        """Return the value of the field in a given *data* object."""
        self.check_available(data)
        original_fields = data.keys()  # Copy
        if self._function is NullFunc:
            raise RuntimeError(
                "Something has gone terribly wrong, _function is NullFunc "
                + f"for {self.name}"
            )
        with self.unit_registry(data):
            dd = self._function(self, data)
        for field_name in data.keys():
            if field_name not in original_fields:
                del data[field_name]
        return dd

    def get_source(self):
        """
        Return a string containing the source of the function (if possible.)
        """
        return inspect.getsource(self._function)

    def get_label(self, projected=False):
        """
        Return a data label for the given field, including units.
        """
        name = self.name[1]
        if self.display_name is not None:
            name = self.display_name

        # Start with the field name
        data_label = r"$\rm{%s}" % name

        # Grab the correct units
        if projected:
            raise NotImplementedError
        else:
            if self.ds is not None:
                units = Unit(self.units, registry=self.ds.unit_registry)
            else:
                units = Unit(self.units)
        # Add unit label
        if not units.is_dimensionless:
            data_label += r"\ \ \left(%s\right)" % (units.latex_representation())

        data_label += r"$"
        return data_label

    @property
    def alias_field(self) -> bool:
        issue_deprecation_warning(
            "DerivedField.alias_field is a deprecated equivalent to DerivedField.is_alias ",
            since="4.1.0",
        )
        return self.is_alias

    @property
    def is_alias(self) -> bool:
        return self._shared_aliases_list.index(self) > 0

    def is_alias_to(self, other: "DerivedField") -> bool:
        return self._shared_aliases_list is other._shared_aliases_list

    @property
    def alias_name(self) -> Optional[Tuple[str, str]]:
        if self.is_alias:
            return self._shared_aliases_list[0].name
        return None

    def __repr__(self):
        if self._function is NullFunc:
            s = "On-Disk Field "
        elif self.is_alias:
            s = f"Alias Field for {self.alias_name!r} "
        else:
            s = "Derived Field "
        s += f"{self.name!r}: (units: {self.units!r}"
        if self.display_name is not None:
            s += f", display_name: {self.display_name!r}"
        if self.sampling_type == "particle":
            s += ", particle field"
        s += ")"
        return s

    def _is_ion(self):
        p = re.compile("_p[0-9]+_")
        result = False
        if p.search(self.name[1]) is not None:
            result = True
        return result

    def _ion_to_label(self):
        # check to see if the output format has changed
        if self.ds is not None:
            self._ionization_label_format = self.ds._ionization_label_format

        pnum2rom = {
            "0": "I",
            "1": "II",
            "2": "III",
            "3": "IV",
            "4": "V",
            "5": "VI",
            "6": "VII",
            "7": "VIII",
            "8": "IX",
            "9": "X",
            "10": "XI",
            "11": "XII",
            "12": "XIII",
            "13": "XIV",
            "14": "XV",
            "15": "XVI",
            "16": "XVII",
            "17": "XVIII",
            "18": "XIX",
            "19": "XX",
        }

        # first look for charge to decide if it is an ion
        p = re.compile("_p[0-9]+_")
        m = p.search(self.name[1])
        if m is not None:

            # Find the ionization state
            pstr = m.string[m.start() + 1 : m.end() - 1]
            segments = self.name[1].split("_")

            # find the ionization index
            for i, s in enumerate(segments):
                if s == pstr:
                    ipstr = i

            for i, s in enumerate(segments):
                # If its the species we don't want to change the capitalization
                if i == ipstr - 1:
                    continue
                segments[i] = s.capitalize()

            species = segments[ipstr - 1]

            # If there is a number in the species part of the label
            # that indicates part of a molecule
            symbols = []
            for symb in species:
                # can't just use underscore b/c gets replaced later with space
                if symb.isdigit():
                    symbols.append("latexsub{" + symb + "}")
                else:
                    symbols.append(symb)
            species_label = "".join(symbols)

            # Use roman numerals for ionization

            if self._ionization_label_format == "roman_numeral":
                roman = pnum2rom[pstr[1:]]
                label = (
                    species_label
                    + r"\ "
                    + roman
                    + r"\ "
                    + r"\ ".join(segments[ipstr + 1 :])
                )

            # use +/- for ionization
            else:
                sign = "+" * int(pstr[1:])
                label = (
                    "{"
                    + species_label
                    + "}"
                    + "^{"
                    + sign
                    + "}"
                    + r"\ "
                    + r"\ ".join(segments[ipstr + 1 :])
                )

        else:
            label = self.name[1]
        return label

    def get_latex_display_name(self):
        label = self.display_name
        if label is None:
            if self._is_ion():
                fname = self._ion_to_label()
                label = r"$\rm{" + fname.replace("_", r"\ ") + r"}$"
                label = label.replace("latexsub", "_")
            else:
                label = r"$\rm{" + self.name[1].replace("_", r"\ ").title() + r"}$"
        elif label.find("$") == -1:
            label = label.replace(" ", r"\ ")
            label = r"$\rm{" + label + r"}$"
        return label

    def __copy__(self):
        # a shallow copy doesn't copy the _shared_alias_list attr
        # This method is implemented in support to ParticleFilter.wrap_func
        return type(self)(
            name=self.name,
            sampling_type=self.sampling_type,
            function=self._function,
            units=self.units,
            take_log=self.take_log,
            validators=self.validators,
            vector_field=self.vector_field,
            display_field=self.display_field,
            not_in_all=self.not_in_all,
            display_name=self.display_name,
            output_units=self.output_units,
            dimensions=self.dimensions,
            ds=self.ds,
            nodal_flag=self.nodal_flag,
        )


class FieldValidator:
    pass


class ValidateParameter(FieldValidator):
    def __init__(self, parameters, parameter_values=None):
        """
        This validator ensures that the dataset has a given parameter.

        If *parameter_values* is supplied, this will also ensure that the field
        is available for all permutations of the field parameter.
        """
        FieldValidator.__init__(self)
        self.parameters = list(always_iterable(parameters))
        self.parameter_values = parameter_values

    def __call__(self, data):
        doesnt_have = []
        for p in self.parameters:
            if not data.has_field_parameter(p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsParameter(doesnt_have)
        return True


class ValidateDataField(FieldValidator):
    def __init__(self, field):
        """
        This validator ensures that the output file has a given data field stored
        in it.
        """
        FieldValidator.__init__(self)
        self.fields = list(iter_fields(field))

    def __call__(self, data):
        doesnt_have = []
        if isinstance(data, FieldDetector):
            return True
        for f in self.fields:
            if f not in data.index.field_list:
                doesnt_have.append(f)
        if len(doesnt_have) > 0:
            raise NeedsDataField(doesnt_have)
        return True


class ValidateProperty(FieldValidator):
    def __init__(self, prop):
        """
        This validator ensures that the data object has a given python attribute.
        """
        FieldValidator.__init__(self)
        self.prop = list(always_iterable(prop))

    def __call__(self, data):
        doesnt_have = []
        for p in self.prop:
            if not hasattr(data, p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsProperty(doesnt_have)
        return True


class ValidateSpatial(FieldValidator):
    def __init__(self, ghost_zones=0, fields=None):
        """
        This validator ensures that the data handed to the field is of spatial
        nature -- that is to say, 3-D.
        """
        FieldValidator.__init__(self)
        self.ghost_zones = ghost_zones
        self.fields = fields

    def __call__(self, data):
        # When we say spatial information, we really mean
        # that it has a three-dimensional data structure
        # if isinstance(data, FieldDetector): return True
        if not getattr(data, "_spatial", False):
            raise NeedsGridType(self.ghost_zones, self.fields)
        if self.ghost_zones <= data._num_ghost_zones:
            return True
        raise NeedsGridType(self.ghost_zones, self.fields)


class ValidateGridType(FieldValidator):
    def __init__(self):
        """
        This validator ensures that the data handed to the field is an actual
        grid patch, not a covering grid of any kind.
        """
        FieldValidator.__init__(self)

    def __call__(self, data):
        # We need to make sure that it's an actual AMR grid
        if isinstance(data, FieldDetector):
            return True
        if getattr(data, "_type_name", None) == "grid":
            return True
        raise NeedsOriginalGrid()
