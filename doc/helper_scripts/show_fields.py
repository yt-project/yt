import inspect
from yt.mods import *
from yt.testing import *
import numpy as np
from yt.utilities.cosmology import \
     Cosmology
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.frontends.stream.fields import \
    StreamFieldInfo
from yt.fields.derived_field import NullFunc
from yt.units.yt_array import YTArray, Unit

fields, units = [], []

for fname, (code_units, aliases, dn) in StreamFieldInfo.known_other_fields:
    fields.append(("gas", fname))
    units.append(code_units)
base_ds = fake_random_ds(4, fields = fields, units = units)
base_ds.index
base_ds.cosmological_simulation = 1
base_ds.cosmology = Cosmology()
from yt.config import ytcfg
ytcfg["yt","__withintesting"] = "True"
np.seterr(all = 'ignore')

def _strip_ftype(field):
    if not isinstance(field, tuple):
        return field
    elif field[0] == "all":
        return field
    return field[1]

np.random.seed(int(0x4d3d3d3))
units = [base_ds._get_field_info(*f).units for f in fields]
fields = [_strip_ftype(f) for f in fields]
ds = fake_random_ds(16, fields = fields, units = units)
ds.parameters["HydroMethod"] = "streaming"
ds.parameters["EOSType"] = 1.0
ds.parameters["EOSSoundSpeed"] = 1.0
ds.conversion_factors["Time"] = 1.0
ds.conversion_factors.update( dict((f, 1.0) for f in fields) )
ds.gamma = 5.0/3.0
ds.current_redshift = 0.0001
ds.cosmological_simulation = 1
ds.hubble_constant = 0.7
ds.omega_matter = 0.27
ds.omega_lambda = 0.73
ds.cosmology = Cosmology(hubble_constant=ds.hubble_constant,
                         omega_matter=ds.omega_matter,
                         omega_lambda=ds.omega_lambda,
                         unit_registry=ds.unit_registry)

header = r"""
.. _field-list:

Field List
==========

This is a list of many of the fields available in ``yt``.  We have attempted to
include most of the fields that are accessible through the plugin system,
however it is possible to generate many more permutations, particularly through
vector operations.

Try using the ``ds.field_list`` and ``ds.derived_field_list`` to view the
native and derived fields available for your dataset respectively. For example
to display the native fields in alphabetical order:

.. notebook-cell::

  from yt.mods import *
  ds = load("Enzo_64/DD0043/data0043")
  for i in sorted(ds.field_list):
    print i

"""

print header

seen = []

def print_all_fields(fl):
    for fn in sorted(fl):
        df = fl[fn]
        f = df._function
        s = "%s" % (df.name,)
        print s
        print "+" * len(s)
        print
        if len(df.units) > 0:
            u = Unit(df.units, registry = ds.unit_registry)
            print "   * Units: :math:`%s`" % (u.latex_representation())
        print "   * Particle Type: %s" % (df.particle_type)
        print
        print "**Field Source**"
        print
        if f == NullFunc:
            print "No source available."
            print
            continue
        else:
            print ".. code-block:: python"
            print
            for line in inspect.getsource(f).split("\n"):
                print "  " + line
            print

ds.index
print_all_fields(ds.field_info)
