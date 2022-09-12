from collections import UserDict
from collections.abc import Callable
from numbers import Number as numeric_type
from typing import Optional, Tuple

import numpy as np
from unyt.exceptions import UnitConversionError

from yt._typing import KnownFieldsT
from yt.config import ytcfg
from yt.fields.field_exceptions import NeedsConfiguration
from yt.funcs import mylog, obj_length, only_on_root
from yt.geometry.geometry_handler import is_curvilinear
from yt.units.dimensions import dimensionless  # type: ignore
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.exceptions import (
    YTCoordinateNotImplemented,
    YTDomainOverflow,
    YTFieldNotFound,
)

from .derived_field import DeprecatedFieldFunc, DerivedField, NullFunc, TranslationFunc
from .field_plugin_registry import field_plugins
from .particle_fields import (
    add_union_field,
    particle_deposition_functions,
    particle_scalar_functions,
    particle_vector_functions,
    sph_whitelist_fields,
    standard_particle_fields,
)


class FieldInfoContainer(UserDict):
    """
    This is a generic field container.  It contains a list of potential derived
    fields, all of which know how to act on a data object and return a value.
    This object handles converting units as well as validating the availability
    of a given field.

    """

    fallback = None
    known_other_fields: KnownFieldsT = ()
    known_particle_fields: KnownFieldsT = ()
    extra_union_fields: Tuple[Tuple[str, str], ...] = ()

    def __init__(self, ds, field_list, slice_info=None):
        super().__init__()
        self._show_field_errors = []
        self.ds = ds
        # Now we start setting things up.
        self.field_list = field_list
        self.slice_info = slice_info
        self.field_aliases = {}
        self.species_names = []
        if ds is not None and is_curvilinear(ds.geometry):
            self.curvilinear = True
        else:
            self.curvilinear = False
        self.setup_fluid_aliases()

    def setup_fluid_fields(self):
        pass

    def setup_fluid_index_fields(self):
        # Now we get all our index types and set up aliases to them
        if self.ds is None:
            return
        index_fields = {f for _, f in self if _ == "index"}
        for ftype in self.ds.fluid_types:
            if ftype in ("index", "deposit"):
                continue
            for f in index_fields:
                if (ftype, f) in self:
                    continue
                self.alias((ftype, f), ("index", f))

    def setup_particle_fields(self, ptype, ftype="gas", num_neighbors=64):
        skip_output_units = ("code_length",)
        for f, (units, aliases, dn) in sorted(self.known_particle_fields):
            units = self.ds.field_units.get((ptype, f), units)
            output_units = units
            if (
                f in aliases or ptype not in self.ds.particle_types_raw
            ) and units not in skip_output_units:
                u = Unit(units, registry=self.ds.unit_registry)
                if u.dimensions is not dimensionless:
                    output_units = str(self.ds.unit_system[u.dimensions])
            if (ptype, f) not in self.field_list:
                continue
            self.add_output_field(
                (ptype, f),
                sampling_type="particle",
                units=units,
                display_name=dn,
                output_units=output_units,
            )
            for alias in aliases:
                self.alias((ptype, alias), (ptype, f), units=output_units)

        # We'll either have particle_position or particle_position_[xyz]
        if (ptype, "particle_position") in self.field_list or (
            ptype,
            "particle_position",
        ) in self.field_aliases:
            particle_scalar_functions(
                ptype, "particle_position", "particle_velocity", self
            )
        else:
            # We need to check to make sure that there's a "known field" that
            # overlaps with one of the vector fields.  For instance, if we are
            # in the Stream frontend, and we have a set of scalar position
            # fields, they will overlap with -- and be overridden by -- the
            # "known" vector field that the frontend creates.  So the easiest
            # thing to do is to simply remove the on-disk field (which doesn't
            # exist) and replace it with a derived field.
            if (ptype, "particle_position") in self and self[
                ptype, "particle_position"
            ]._function == NullFunc:
                self.pop((ptype, "particle_position"))
            particle_vector_functions(
                ptype,
                [f"particle_position_{ax}" for ax in "xyz"],
                [f"particle_velocity_{ax}" for ax in "xyz"],
                self,
            )
        particle_deposition_functions(ptype, "particle_position", "particle_mass", self)
        standard_particle_fields(self, ptype)
        # Now we check for any leftover particle fields
        for field in sorted(self.field_list):
            if field in self:
                continue
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] not in self.ds.particle_types:
                continue
            units = self.ds.field_units.get(field, None)
            if units is None:
                try:
                    units = ytcfg.get("fields", *field, "units")
                except KeyError:
                    units = ""
            self.add_output_field(
                field,
                sampling_type="particle",
                units=units,
            )
        self.setup_smoothed_fields(ptype, num_neighbors=num_neighbors, ftype=ftype)

    def setup_extra_union_fields(self, ptype="all"):
        if ptype != "all":
            raise RuntimeError(
                "setup_extra_union_fields is currently"
                + 'only enabled for particle type "all".'
            )
        for units, field in self.extra_union_fields:
            add_union_field(self, ptype, field, units)

    def setup_smoothed_fields(self, ptype, num_neighbors=64, ftype="gas"):
        # We can in principle compute this, but it is not yet implemented.
        if (ptype, "density") not in self or not hasattr(self.ds, "_sph_ptypes"):
            return
        new_aliases = []
        for ptype2, alias_name in list(self):
            if ptype2 != ptype:
                continue
            if alias_name not in sph_whitelist_fields:
                if alias_name.startswith("particle_"):
                    pass
                else:
                    continue
            uni_alias_name = alias_name
            if "particle_position_" in alias_name:
                uni_alias_name = alias_name.replace("particle_position_", "")
            elif "particle_" in alias_name:
                uni_alias_name = alias_name.replace("particle_", "")
            new_aliases.append(
                (
                    (ftype, uni_alias_name),
                    (ptype, alias_name),
                )
            )
            new_aliases.append(
                (
                    (ptype, uni_alias_name),
                    (ptype, alias_name),
                )
            )
            for alias, source in new_aliases:
                self.alias(alias, source)

    # Collect the names for all aliases if geometry is curvilinear
    def get_aliases_gallery(self):
        aliases_gallery = []
        known_other_fields = dict(self.known_other_fields)
        if self.curvilinear:
            for field in sorted(self.field_list):
                if field[0] in self.ds.particle_types:
                    continue
                args = known_other_fields.get(field[1], ("", [], None))
                units, aliases, display_name = args
                for alias in aliases:
                    aliases_gallery.append(alias)
        return aliases_gallery

    def setup_fluid_aliases(self, ftype="gas"):
        known_other_fields = dict(self.known_other_fields)

        # For non-Cartesian geometry, convert alias of vector fields to
        # curvilinear coordinates
        aliases_gallery = self.get_aliases_gallery()

        for field in sorted(self.field_list):
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] in self.ds.particle_types:
                continue
            args = known_other_fields.get(field[1], None)
            if args is not None:
                units, aliases, display_name = args
            else:
                try:
                    node = ytcfg.get("fields", *field).as_dict()
                except KeyError:
                    node = dict()

                units = node.get("units", "")
                aliases = node.get("aliases", [])
                display_name = node.get("display_name", None)

            # We allow field_units to override this.  First we check if the
            # field *name* is in there, then the field *tuple*.
            units = self.ds.field_units.get(field[1], units)
            units = self.ds.field_units.get(field, units)
            if not isinstance(units, str) and args[0] != "":
                units = f"(({args[0]})*{units})"
            if (
                isinstance(units, (numeric_type, np.number, np.ndarray))
                and args[0] == ""
                and units != 1.0
            ):
                mylog.warning(
                    "Cannot interpret units: %s * %s, setting to dimensionless.",
                    units,
                    args[0],
                )
                units = ""
            elif units == 1.0:
                units = ""
            self.add_output_field(
                field, sampling_type="cell", units=units, display_name=display_name
            )
            axis_names = self.ds.coordinates.axis_order
            for alias in aliases:
                if (
                    self.curvilinear
                ):  # For non-Cartesian geometry, convert vector aliases

                    if alias[-2:] not in ["_x", "_y", "_z"]:
                        to_convert = False
                    else:
                        for suffix in ["x", "y", "z"]:
                            if f"{alias[:-2]}_{suffix}" not in aliases_gallery:
                                to_convert = False
                                break
                        to_convert = True
                    if to_convert:
                        if alias[-2:] == "_x":
                            alias = f"{alias[:-2]}_{axis_names[0]}"
                        elif alias[-2:] == "_y":
                            alias = f"{alias[:-2]}_{axis_names[1]}"
                        elif alias[-2:] == "_z":
                            alias = f"{alias[:-2]}_{axis_names[2]}"
                self.alias((ftype, alias), field)

    @staticmethod
    def _sanitize_sampling_type(sampling_type: str) -> str:
        """Detect conflicts between deprecated and new parameters to specify the
        sampling type in a new field.

        This is a helper function to add_field methods.

        Parameters
        ----------
        sampling_type : str
            One of "cell", "particle" or "local" (case insensitive)

        Raises
        ------
        ValueError
            For unsupported values in sampling_type
        """
        if not isinstance(sampling_type, str):
            raise TypeError("sampling_type should be a string.")

        sampling_type = sampling_type.lower()
        acceptable_samplings = ("cell", "particle", "local")
        if sampling_type not in acceptable_samplings:
            raise ValueError(
                f"Received invalid sampling type {sampling_type!r}. "
                f"Expected any of {acceptable_samplings}"
            )
        return sampling_type

    def add_field(
        self,
        name: Tuple[str, str],
        function: Callable,
        sampling_type: str,
        *,
        alias: Optional[DerivedField] = None,
        force_override: bool = False,
        **kwargs,
    ) -> None:
        """
        Add a new field, along with supplemental metadata, to the list of
        available fields.  This respects a number of arguments, all of which
        are passed on to the constructor for
        :class:`~yt.data_objects.api.DerivedField`.

        Parameters
        ----------

        name : tuple[str, str]
           field (or particle) type, field name
        function : callable
           A function handle that defines the field.  Should accept
           arguments (field, data)
        sampling_type: str
           "cell" or "particle" or "local"
        force_override: bool
           If False (default), an error will be raised if a field of the same name already exists.
        alias: DerivedField (optional):
           existing field to be aliased
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^). If set to "auto" the units
           will be inferred from the return value of the field function.
        take_log : bool
           Describes whether the field should be logged
        validators : list
           A list of :class:`FieldValidator` objects
        vector_field : bool
           Describes the dimensionality of the field.  Currently unused.
        display_name : str
           A name used in the plots

        """
        # Handle the case where the field has already been added.
        if not force_override and name in self:
            return

        kwargs.setdefault("ds", self.ds)

        sampling_type = self._sanitize_sampling_type(sampling_type)

        if (
            not isinstance(name, str)
            and obj_length(name) == 2
            and all(isinstance(e, str) for e in name)
        ):
            self[name] = DerivedField(
                name, sampling_type, function, alias=alias, **kwargs
            )
        else:
            raise ValueError(f"Expected name to be a tuple[str, str], got {name}")

    def load_all_plugins(self, ftype: Optional[str] = "gas"):
        if ftype is None:
            return
        mylog.debug("Loading field plugins for field type: %s.", ftype)
        loaded = []
        for n in sorted(field_plugins):
            loaded += self.load_plugin(n, ftype)
            only_on_root(mylog.debug, "Loaded %s (%s new fields)", n, len(loaded))
        self.find_dependencies(loaded)

    def load_plugin(self, plugin_name, ftype="gas", skip_check=False):
        if callable(plugin_name):
            f = plugin_name
        else:
            f = field_plugins[plugin_name]
        orig = set(self.items())
        f(self, ftype, slice_info=self.slice_info)
        loaded = [n for n, v in set(self.items()).difference(orig)]
        return loaded

    def find_dependencies(self, loaded):
        deps, unavailable = self.check_derived_fields(loaded)
        self.ds.field_dependencies.update(deps)
        # Note we may have duplicated
        dfl = set(self.ds.derived_field_list).union(deps.keys())
        self.ds.derived_field_list = sorted(dfl)
        return loaded, unavailable

    def add_output_field(self, name, sampling_type, **kwargs):
        if name[1] == "density":
            if name in self:
                # this should not happen, but it does
                # it'd be best to raise an error here but
                # it may take a while to cleanup internal issues
                return
        kwargs.setdefault("ds", self.ds)
        self[name] = DerivedField(name, sampling_type, NullFunc, **kwargs)

    def alias(
        self,
        alias_name: Tuple[str, str],
        original_name: Tuple[str, str],
        units: Optional[str] = None,
        deprecate: Optional[Tuple[str, Optional[str]]] = None,
    ):
        """
        Alias one field to another field.

        Parameters
        ----------
        alias_name : Tuple[str, str]
            The new field name.
        original_name : Tuple[str, str]
            The field to be aliased.
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^). If set to "auto" the units
           will be inferred from the return value of the field function.
        deprecate : tuple[str, str | None] | None
            If this is set, then the tuple contains two string version
            numbers: the first marking the version when the field was
            deprecated, and the second marking when the field will be
            removed.
        """
        if original_name not in self:
            return
        if units is None:
            # We default to CGS here, but in principle, this can be pluggable
            # as well.

            # self[original_name].units may be set to `None` at this point
            # to signal that units should be autoset later
            oru = self[original_name].units
            if oru is None:
                units = None
            else:
                u = Unit(oru, registry=self.ds.unit_registry)
                if u.dimensions is not dimensionless:
                    units = str(self.ds.unit_system[u.dimensions])
                else:
                    units = oru

        self.field_aliases[alias_name] = original_name
        function = TranslationFunc(original_name)
        if deprecate is not None:
            self.add_deprecated_field(
                alias_name,
                function=function,
                sampling_type=self[original_name].sampling_type,
                display_name=self[original_name].display_name,
                units=units,
                since=deprecate[0],
                removal=deprecate[1],
                ret_name=original_name,
            )
        else:
            self.add_field(
                alias_name,
                function=function,
                sampling_type=self[original_name].sampling_type,
                display_name=self[original_name].display_name,
                units=units,
                alias=self[original_name],
            )

    def add_deprecated_field(
        self,
        name,
        function,
        sampling_type,
        since,
        removal=None,
        ret_name=None,
        **kwargs,
    ):
        """
        Add a new field which is deprecated, along with supplemental metadata,
        to the list of available fields.  This respects a number of arguments,
        all of which are passed on to the constructor for
        :class:`~yt.data_objects.api.DerivedField`.

        Parameters
        ----------
        name : str
           is the name of the field.
        function : callable
           A function handle that defines the field.  Should accept
           arguments (field, data)
        sampling_type : str
           "cell" or "particle" or "local"
        since : str
            The version string marking when this field was deprecated.
        removal : str
            The version string marking when this field will be removed.
        ret_name : str
            The name of the field which will actually be returned, used
            only by :meth:`~yt.fields.field_info_container.FieldInfoContainer.alias`.
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^). If set to "auto" the units
           will be inferred from the return value of the field function.
        take_log : bool
           Describes whether the field should be logged
        validators : list
           A list of :class:`FieldValidator` objects
        vector_field : bool
           Describes the dimensionality of the field.  Currently unused.
        display_name : str
           A name used in the plots
        """
        if ret_name is None:
            ret_name = name
        self.add_field(
            name,
            function=DeprecatedFieldFunc(ret_name, function, since, removal),
            sampling_type=sampling_type,
            **kwargs,
        )

    def has_key(self, key):
        # This gets used a lot
        if key in self:
            return True
        if self.fallback is None:
            return False
        return key in self.fallback

    def __missing__(self, key):
        if self.fallback is None:
            raise KeyError(f"No field named {key}")
        return self.fallback[key]

    @classmethod
    def create_with_fallback(cls, fallback, name=""):
        obj = cls()
        obj.fallback = fallback
        obj.name = name
        return obj

    def __contains__(self, key):
        if super().__contains__(key):
            return True
        if self.fallback is None:
            return False
        return key in self.fallback

    def __iter__(self):
        yield from super().__iter__()
        if self.fallback is not None:
            yield from self.fallback

    def keys(self):
        keys = super().keys()
        if self.fallback:
            keys += list(self.fallback.keys())
        return keys

    def check_derived_fields(self, fields_to_check=None):

        # The following exceptions lists were obtained by expanding an
        # all-catching `except Exception`.
        # We define
        # - a blacklist (exceptions that we know should be caught)
        # - a whitelist (exceptions that should be handled)
        # - a greylist (exceptions that may be covering bugs but should be checked)
        # See https://github.com/yt-project/yt/issues/2853
        # in the long run, the greylist should be removed
        blacklist = ()
        whitelist = (NotImplementedError,)
        greylist = (
            YTFieldNotFound,
            YTDomainOverflow,
            YTCoordinateNotImplemented,
            NeedsConfiguration,
            TypeError,
            ValueError,
            IndexError,
            AttributeError,
            KeyError,
            # code smells -> those are very likely bugs
            UnitConversionError,  # solved in GH PR 2897 ?
            # RecursionError is clearly a bug, and was already solved once
            # in GH PR 2851
            RecursionError,
        )

        deps = {}
        unavailable = []
        fields_to_check = fields_to_check or list(self.keys())
        for field in fields_to_check:
            fi = self[field]
            try:
                # fd: field detector
                fd = fi.get_dependencies(ds=self.ds)
            except blacklist as err:
                print(f"{err.__class__} raised for field {field}")
                raise SystemExit(1) from err
            except (*whitelist, *greylist) as e:
                if field in self._show_field_errors:
                    raise
                if not isinstance(e, YTFieldNotFound):
                    # if we're doing field tests, raise an error
                    # see yt.fields.tests.test_fields
                    if hasattr(self.ds, "_field_test_dataset"):
                        raise
                    mylog.debug(
                        "Raises %s during field %s detection.", str(type(e)), field
                    )
                self.pop(field)
                continue
            # This next bit checks that we can't somehow generate everything.
            # We also manually update the 'requested' attribute
            missing = not all(f in self.field_list for f in fd.requested)
            if missing:
                self.pop(field)
                unavailable.append(field)
                continue
            fd.requested = set(fd.requested)
            deps[field] = fd
            mylog.debug("Succeeded with %s (needs %s)", field, fd.requested)

        # now populate the derived field list with results
        # this violates isolation principles and should be refactored
        dfl = set(self.ds.derived_field_list).union(deps.keys())
        dfl = sorted(dfl)

        if not hasattr(self.ds.index, "meshes"):
            # the meshes attribute characterizes an unstructured-mesh data structure

            # ideally this filtering should not be required
            # and this could maybe be handled in fi.get_dependencies
            # but it's a lot easier to do here

            filtered_dfl = []
            for field in dfl:
                try:
                    ftype, fname = field
                    if "vertex" in fname:
                        continue
                except ValueError:
                    # in very rare cases, there can a field represented by a single
                    # string, like "emissivity"
                    # this try block _should_ be removed and the error fixed upstream
                    # for reference, a test that would break is
                    # yt/data_objects/tests/test_fluxes.py::ExporterTests
                    pass
                filtered_dfl.append(field)
            dfl = filtered_dfl

        self.ds.derived_field_list = dfl
        self._set_linear_fields()
        return deps, unavailable

    def _set_linear_fields(self):
        """
        Sets which fields use linear as their default scaling in Profiles and
        PhasePlots. Default for all fields is set to log, so this sets which
        are linear.  For now, set linear to geometric fields: position and
        velocity coordinates.
        """
        non_log_prefixes = ("", "velocity_", "particle_position_", "particle_velocity_")
        coords = ("x", "y", "z")
        non_log_fields = [
            prefix + coord for prefix in non_log_prefixes for coord in coords
        ]
        for field in self.ds.derived_field_list:
            if field[1] in non_log_fields:
                self[field].take_log = False
