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

from yt.funcs import mylog
from yt.utilities.units import Unit
from yt.data_objects.yt_array import YTArray
from .derived_field import \
    DerivedField, \
    NullFunc, \
    TranslationFunc
from .field_detector import \
    FieldDetector
from yt.utilities.exceptions import \
    YTFieldNotFound
from .field_plugin_registry import \
    field_plugins
from yt.utilities.units import \
    Unit

class FieldInfoContainer(dict): # Resistance has utility
    """
    This is a generic field container.  It contains a list of potential derived
    fields, all of which know how to act on a data object and return a value.
    This object handles converting units as well as validating the availability
    of a given field.

    """
    fallback = None

    def __init__(self, pf, field_list, slice_info = None):
        self.pf = pf
        # Now we start setting things up.
        self.field_list = field_list
        self.slice_info = slice_info

    def add_field(self, name, function=None, **kwargs):
        """
        Add a new field, along with supplemental metadata, to the list of
        available fields.  This respects a number of arguments, all of which
        are passed on to the constructor for
        :class:`~yt.data_objects.api.DerivedField`.

        """
        if function == None:
            def create_function(function):
                self[name] = DerivedField(name, function, **kwargs)
                return function
            return create_function
        self[name] = DerivedField(name, function, **kwargs)

    def load_all_plugins(self, ftype="gas"):
        for n in sorted(field_plugins):
            loaded, not_loaded = self.load_plugin(n, ftype)
            mylog.info("Loaded %s (%s new fields)",
                n, len(loaded))

    def load_plugin(self, plugin_name, ftype = "gas"):
        orig = set(self.items())
        f = field_plugins[plugin_name]
        f(self, ftype, slice_info = self.slice_info)
        loaded = [n for n, v in set(self.items()).difference(orig)]
        deps, unavailable = self.check_derived_fields(loaded)
        self.pf.field_dependencies.update(deps)
        # Note we may have duplicated
        dfl = set(self.pf.h.derived_field_list).union(deps.keys())
        self.pf.h.derived_field_list = list(sorted(dfl))
        return loaded, unavailable

    def add_output_field(self, name, **kwargs):
        self[name] = DerivedField(name, NullFunc, **kwargs)

    def alias(self, alias_name, original_name, units = None):
        if units is None:
            # We default to CGS here, but in principle, this can be pluggable
            # as well.
            u = Unit(self[original_name].units)
            units = str(u.get_cgs_equivalent())
        self.add_field(alias_name,
            function = TranslationFunc(original_name),
            units = units)

    def add_grad(self, field, **kwargs):
        """
        Creates the partial derivative of a given field. This function will
        autogenerate the names of the gradient fields.

        """
        sl = slice(2,None,None)
        sr = slice(None,-2,None)

        def _gradx(f, data):
            grad = data[field][sl,1:-1,1:-1] - data[field][sr,1:-1,1:-1]
            grad /= 2.0*data["dx"].flat[0]
            g = np.zeros(data[field].shape, dtype='float64')
            g[1:-1,1:-1,1:-1] = grad
            return g

        def _grady(f, data):
            grad = data[field][1:-1,sl,1:-1] - data[field][1:-1,sr,1:-1]
            grad /= 2.0*data["dy"].flat[0]
            g = np.zeros(data[field].shape, dtype='float64')
            g[1:-1,1:-1,1:-1] = grad
            return g

        def _gradz(f, data):
            grad = data[field][1:-1,1:-1,sl] - data[field][1:-1,1:-1,sr]
            grad /= 2.0*data["dz"].flat[0]
            g = np.zeros(data[field].shape, dtype='float64')
            g[1:-1,1:-1,1:-1] = grad
            return g

        d_kwargs = kwargs.copy()
        if "display_name" in kwargs: del d_kwargs["display_name"]

        for ax in "xyz":
            if "display_name" in kwargs:
                disp_name = r"%s\_%s" % (kwargs["display_name"], ax)
            else:
                disp_name = r"\partial %s/\partial %s" % (field, ax)
            name = "Grad_%s_%s" % (field, ax)
            self[name] = DerivedField(name, function=eval('_grad%s' % ax),
                         take_log=False, validators=[ValidateSpatial(1,[field])],
                         display_name = disp_name, **d_kwargs)

        def _grad(f, data) :
            a = np.power(data["Grad_%s_x" % field],2)
            b = np.power(data["Grad_%s_y" % field],2)
            c = np.power(data["Grad_%s_z" % field],2)
            norm = np.sqrt(a+b+c)
            return norm

        if "display_name" in kwargs:
            disp_name = kwargs["display_name"]
        else:
            disp_name = r"\Vert\nabla %s\Vert" % (field)
        name = "Grad_%s" % field
        self[name] = DerivedField(name, function=_grad, take_log=False,
                                  display_name = disp_name, **d_kwargs)
        mylog.info("Added new fields: Grad_%s_x, Grad_%s_y, Grad_%s_z, Grad_%s" \
                   % (field, field, field, field))

    def has_key(self, key):
        # This gets used a lot
        if key in self: return True
        if self.fallback is None: return False
        return self.fallback.has_key(key)

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
            keys += self.fallback.keys()
        return keys

    def check_derived_fields(self, fields_to_check = None):
        deps = {}
        unavailable = []
        fields_to_check = fields_to_check or self.keys()
        for field in fields_to_check:
            mylog.debug("Checking %s", field)
            if field not in self: raise RuntimeError
            fi = self[field]
            try:
                fd = fi.get_dependencies(pf = self.pf)
            except Exception as e:
                if "mach_number" in field:
                    raise
                if type(e) != YTFieldNotFound:
                    mylog.debug("Raises %s during field %s detection.",
                                str(type(e)), field)
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
        dfl = set(self.pf.h.derived_field_list).union(deps.keys())
        self.pf.h.derived_field_list = list(sorted(dfl))
        return deps, unavailable
