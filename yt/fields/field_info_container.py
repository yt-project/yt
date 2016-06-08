"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from numbers import Number as numeric_type

from yt.extern.six import string_types
from yt.funcs import mylog, only_on_root
from yt.units.unit_object import Unit
from .derived_field import \
    DerivedField, \
    NullFunc, \
    TranslationFunc
from yt.utilities.exceptions import \
    YTFieldNotFound
from .field_plugin_registry import \
    field_plugins
from .particle_fields import \
    add_union_field, \
    particle_deposition_functions, \
    particle_vector_functions, \
    particle_scalar_functions, \
    standard_particle_fields, \
    add_volume_weighted_smoothed_field, \
    sph_whitelist_fields

def tupleize(inp):
    if isinstance(inp, tuple):
        return inp
    # prepending with a '?' ensures that the sort order is the same in py2 and
    # py3, since names of field types shouldn't begin with punctuation
    return ('?', inp, )

class FieldInfoContainer(dict):
    """
    This is a generic field container.  It contains a list of potential derived
    fields, all of which know how to act on a data object and return a value.
    This object handles converting units as well as validating the availability
    of a given field.

    """
    fallback = None
    known_other_fields = ()
    known_particle_fields = ()
    extra_union_fields = ()

    def __init__(self, ds, field_list, slice_info = None):
        self._show_field_errors = []
        self.ds = ds
        # Now we start setting things up.
        self.field_list = field_list
        self.slice_info = slice_info
        self.field_aliases = {}
        self.species_names = []
        self.setup_fluid_aliases()

    def setup_fluid_fields(self):
        pass

    def setup_fluid_index_fields(self):
        # Now we get all our index types and set up aliases to them
        if self.ds is None: return
        index_fields = set([f for _, f in self if _ == "index"])
        for ftype in self.ds.fluid_types:
            if ftype in ("index", "deposit"): continue
            for f in index_fields:
                if (ftype, f) in self: continue
                self.alias((ftype, f), ("index", f))
                # Different field types have different default units.
                # We want to make sure the aliased field will have
                # the same units as the "index" field.
                self[(ftype, f)].units = self["index", f].units

    def setup_particle_fields(self, ptype, ftype='gas', num_neighbors=64 ):
        skip_output_units = ("code_length",)
        for f, (units, aliases, dn) in sorted(self.known_particle_fields):
            units = self.ds.field_units.get((ptype, f), units)
            if (f in aliases or ptype not in self.ds.particle_types_raw) and \
                units not in skip_output_units:
                u = Unit(units, registry = self.ds.unit_registry)
                output_units = str(self.ds.unit_system[u.dimensions])
            else:
                output_units = units
            if (ptype, f) not in self.field_list:
                continue
            self.add_output_field((ptype, f),
                units = units, particle_type = True,
                display_name = dn, output_units = output_units)
            for alias in aliases:
                self.alias((ptype, alias), (ptype, f), units = output_units)

        # We'll either have particle_position or particle_position_[xyz]
        if (ptype, "particle_position") in self.field_list or \
           (ptype, "particle_position") in self.field_aliases:
            particle_scalar_functions(ptype,
                   "particle_position", "particle_velocity",
                   self)
        else:
            # We need to check to make sure that there's a "known field" that
            # overlaps with one of the vector fields.  For instance, if we are
            # in the Stream frontend, and we have a set of scalar position
            # fields, they will overlap with -- and be overridden by -- the
            # "known" vector field that the frontend creates.  So the easiest
            # thing to do is to simply remove the on-disk field (which doesn't
            # exist) and replace it with a derived field.
            if (ptype, "particle_position") in self and \
                 self[ptype, "particle_position"]._function == NullFunc:
                self.pop((ptype, "particle_position"))
            particle_vector_functions(ptype,
                    ["particle_position_%s" % ax for ax in 'xyz'],
                    ["particle_velocity_%s" % ax for ax in 'xyz'],
                    self)
        particle_deposition_functions(ptype, "particle_position",
            "particle_mass", self)
        standard_particle_fields(self, ptype)
        # Now we check for any leftover particle fields
        for field in sorted(self.field_list):
            if field in self: continue
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] not in self.ds.particle_types:
                continue
            self.add_output_field(field, 
                                  units = self.ds.field_units.get(field, ""),
                                  particle_type = True)
        self.setup_smoothed_fields(ptype, 
                                   num_neighbors=num_neighbors,
                                   ftype=ftype)

    def setup_extra_union_fields(self, ptype="all"):
        if ptype != "all":
            raise RuntimeError("setup_extra_union_fields is currently" + 
                               "only enabled for particle type \"all\".")
        for units, field in self.extra_union_fields:
            add_union_field(self, ptype, field, units)

    def setup_smoothed_fields(self, ptype, num_neighbors = 64, ftype = "gas"):
        # We can in principle compute this, but it is not yet implemented.
        if (ptype, "density") not in self:
            return
        if (ptype, "smoothing_length") in self:
            sml_name = "smoothing_length"
        else:
            sml_name = None
        new_aliases = []
        for ptype2, alias_name in list(self):
            if ptype2 != ptype:
                continue
            if alias_name not in sph_whitelist_fields:
                continue
            fn = add_volume_weighted_smoothed_field(
                ptype, "particle_position", "particle_mass",
                sml_name, "density", alias_name, self,
                num_neighbors)
            if 'particle_' in alias_name:
                alias_name = alias_name.replace('particle_', '')
            new_aliases.append(((ftype, alias_name), fn[0]))
        for alias, source in new_aliases:
            #print "Aliasing %s => %s" % (alias, source)
            self.alias(alias, source)

    def setup_fluid_aliases(self):
        known_other_fields = dict(self.known_other_fields)
        for field in sorted(self.field_list):
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] in self.ds.particle_types:
                continue
            args = known_other_fields.get(
                field[1], ("", [], None))
            units, aliases, display_name = args
            # We allow field_units to override this.  First we check if the
            # field *name* is in there, then the field *tuple*.
            units = self.ds.field_units.get(field[1], units)
            units = self.ds.field_units.get(field, units)
            if not isinstance(units, string_types) and args[0] != "":
                units = "((%s)*%s)" % (args[0], units)
            if isinstance(units, (numeric_type, np.number, np.ndarray)) and \
                args[0] == "" and units != 1.0:
                mylog.warning("Cannot interpret units: %s * %s, " +
                              "setting to dimensionless.", units, args[0])
                units = ""
            elif units == 1.0:
                units = ""
            self.add_output_field(field, units = units,
                                  display_name = display_name)
            for alias in aliases:
                self.alias(("gas", alias), field)

    def add_field(self, name, function=None, **kwargs):
        """
        Add a new field, along with supplemental metadata, to the list of
        available fields.  This respects a number of arguments, all of which
        are passed on to the constructor for
        :class:`~yt.data_objects.api.DerivedField`.

        Parameters
        ----------

        name : str
           is the name of the field.
        function : callable
           A function handle that defines the field.  Should accept
           arguments (field, data)
        units : str
           A plain text string encoding the unit.  Powers must be in
           python syntax (** instead of ^). If set to "auto" the units
           will be inferred from the return value of the field function.
        take_log : bool
           Describes whether the field should be logged
        validators : list
           A list of :class:`FieldValidator` objects
        particle_type : bool
           Is this a particle (1D) field?
        vector_field : bool
           Describes the dimensionality of the field.  Currently unused.
        display_name : str
           A name used in the plots

        """
        override = kwargs.pop("force_override", False)
        # Handle the case where the field has already been added.
        if not override and name in self:
            # See below.
            if function is None:
                def create_function(f):
                    return f
                return create_function
            return
        # add_field can be used in two different ways: it can be called
        # directly, or used as a decorator. If called directly, the
        # function will be passed in as an argument, and we simply create
        # the derived field and exit. If used as a decorator, function will
        # be None. In that case, we return a function that will be applied
        # to the function that the decorator is applied to.
        kwargs.setdefault('ds', self.ds)
        if function is None:
            def create_function(f):
                self[name] = DerivedField(name, f, **kwargs)
                return f
            return create_function

        if isinstance(name, tuple):
            self[name] = DerivedField(name, function, **kwargs)
            return

        if kwargs.get("particle_type", False):
            ftype = 'all'
        else:
            ftype = self.ds.default_fluid_type

        if (ftype, name) not in self:
            tuple_name = (ftype, name)
            self[tuple_name] = DerivedField(tuple_name, function, **kwargs)
            self.alias(name, tuple_name)
        else:
            self[name] = DerivedField(name, function, **kwargs)

    def load_all_plugins(self, ftype="gas"):
        loaded = []
        for n in sorted(field_plugins):
            loaded += self.load_plugin(n, ftype)
            only_on_root(mylog.debug, "Loaded %s (%s new fields)",
                         n, len(loaded))
        self.find_dependencies(loaded)

    def load_plugin(self, plugin_name, ftype = "gas", skip_check = False):
        if callable(plugin_name):
            f = plugin_name
        else:
            f = field_plugins[plugin_name]
        orig = set(self.items())
        f(self, ftype, slice_info = self.slice_info)
        loaded = [n for n, v in set(self.items()).difference(orig)]
        return loaded

    def find_dependencies(self, loaded):
        deps, unavailable = self.check_derived_fields(loaded)
        self.ds.field_dependencies.update(deps)
        # Note we may have duplicated
        dfl = set(self.ds.derived_field_list).union(deps.keys())
        self.ds.derived_field_list = list(sorted(dfl, key=tupleize))
        return loaded, unavailable

    def add_output_field(self, name, **kwargs):
        kwargs.setdefault('ds', self.ds)
        self[name] = DerivedField(name, NullFunc, **kwargs)

    def alias(self, alias_name, original_name, units = None):
        if original_name not in self: return
        if units is None:
            # We default to CGS here, but in principle, this can be pluggable
            # as well.
            u = Unit(self[original_name].units,
                      registry = self.ds.unit_registry)
            units = str(self.ds.unit_system[u.dimensions])
        self.field_aliases[alias_name] = original_name
        self.add_field(alias_name,
            function = TranslationFunc(original_name),
            particle_type = self[original_name].particle_type,
            display_name = self[original_name].display_name,
            units = units)

    def has_key(self, key):
        # This gets used a lot
        if key in self: return True
        if self.fallback is None: return False
        return key in self.fallback

    def __missing__(self, key):
        if self.fallback is None:
            raise KeyError("No field named %s" % (key,))
        return self.fallback[key]

    @classmethod
    def create_with_fallback(cls, fallback, name = ""):
        obj = cls()
        obj.fallback = fallback
        obj.name = name
        return obj

    def __contains__(self, key):
        if dict.__contains__(self, key): return True
        if self.fallback is None: return False
        return key in self.fallback

    def __iter__(self):
        for f in dict.__iter__(self):
            yield f
        if self.fallback is not None:
            for f in self.fallback: yield f

    def keys(self):
        keys = dict.keys(self)
        if self.fallback:
            keys += list(self.fallback.keys())
        return keys

    def check_derived_fields(self, fields_to_check = None):
        deps = {}
        unavailable = []
        fields_to_check = fields_to_check or list(self.keys())
        for field in fields_to_check:
            mylog.debug("Checking %s", field)
            if field not in self: raise RuntimeError
            fi = self[field]
            try:
                fd = fi.get_dependencies(ds = self.ds)
            except Exception as e:
                if field in self._show_field_errors:
                    raise
                if type(e) != YTFieldNotFound:
                    mylog.debug("Raises %s during field %s detection.",
                                str(type(e)), field)
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
        dfl = set(self.ds.derived_field_list).union(deps.keys())
        self.ds.derived_field_list = list(sorted(dfl, key=tupleize))
        return deps, unavailable
