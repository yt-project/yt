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
from yt.frontends.api import _frontends
from yt.fields.derived_field import NullFunc
import yt.frontends as frontends_module
from yt.units.yt_array import YTArray, Unit
from yt.units import dimensions

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
for my_unit in ["m", "pc", "AU", "au"]:
    new_unit = "%scm" % my_unit
    ds.unit_registry.add(new_unit, base_ds.unit_registry.lut[my_unit][0],
                         dimensions.length, "\\rm{%s}/(1+z)" % my_unit)



header = r"""
.. _field-list:

Field List
==========

This is a list of many of the fields available in yt.  We have attempted to
include most of the fields that are accessible through the plugin system, as 
well as the fields that are known by the frontends, however it is possible to 
generate many more permutations, particularly through vector operations. For 
more information about the fields framework, see :ref:`fields`.

Some fields are recognized by specific frontends only. These are typically 
fields like density and temperature that have their own names and units in 
the different frontend datasets. Often, these fields are aliased to their 
yt-named counterpart fields (typically 'gas' fieldtypes). For example, in 
the ``FLASH`` frontend, the ``dens`` field (i.e. ``(flash, dens)``) is aliased 
to the gas field density (i.e. ``(gas, density)``), similarly ``(flash, velx)`` 
is aliased to ``(gas, velocity_x)``, and so on. In what follows, if a field 
is aliased it will be noted.

Try using the ``ds.field_list`` and ``ds.derived_field_list`` to view the
native and derived fields available for your dataset respectively. For example
to display the native fields in alphabetical order:

.. notebook-cell::

  import yt
  ds = yt.load("Enzo_64/DD0043/data0043")
  for i in sorted(ds.field_list):
    print i

To figure out out what all of the field types here mean, see
:ref:`known-field-types`.

.. contents:: Table of Contents
   :depth: 1
   :local:
   :backlinks: none

.. _yt-fields:

Universal Fields
----------------
"""

footer = """

Index of Fields
---------------

.. contents:: 
   :depth: 3
   :backlinks: none

"""
print header

seen = []

def fix_units(units, in_cgs=False):
    unit_object = Unit(units, registry=ds.unit_registry)
    if in_cgs:
        unit_object = unit_object.get_cgs_equivalent()
    latex = unit_object.latex_representation()
    return latex.replace('\/','~')

def print_all_fields(fl):
    for fn in sorted(fl):
        df = fl[fn]
        f = df._function
        s = "%s" % (df.name,)
        print s
        print "^" * len(s)
        print
        if len(df.units) > 0:
            # Most universal fields are in CGS except for these special fields
            if df.name[1] in ['particle_position', 'particle_position_x', \
                         'particle_position_y', 'particle_position_z', \
                         'entropy', 'kT', 'metallicity', 'dx', 'dy', 'dz',\
                         'cell_volume', 'x', 'y', 'z']:
                print "   * Units: :math:`%s`" % fix_units(df.units)
            else:
                print "   * Units: :math:`%s`" % fix_units(df.units, in_cgs=True)
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

def print_frontend_field(ftype, field, ptype):
    name = field[0]
    units = field[1][0]
    aliases = ["``%s``" % f for f in field[1][1]]
    if ftype is not "particle_type":
        ftype = "'"+ftype+"'"
    s = "(%s, '%s')" % (ftype, name)
    print s
    print "^" * len(s)
    print
    if len(units) > 0:
        print "   * Units: :math:`\mathrm{%s}`" % fix_units(units)
    if len(aliases) > 0:
        print "   * Aliased to: %s" % " ".join(aliases)
    print "   * Particle Type: %s" % (ptype)
    print

current_frontends = [f for f in _frontends if f not in ["stream"]]

for frontend in current_frontends:
    this_f = getattr(frontends_module, frontend)
    field_info_names = [fi for fi in dir(this_f) if "FieldInfo" in fi]
    dataset_names = [dset for dset in dir(this_f) if "Dataset" in dset]

    if frontend == "sph":
        field_info_names = \
          ['TipsyFieldInfo' if 'Tipsy' in d else 'SPHFieldInfo' for d in dataset_names]
    elif frontend == "boxlib":
        field_info_names = []
        for d in dataset_names:
            if "Maestro" in d:  
                field_info_names.append("MaestroFieldInfo")
            elif "Castro" in d: 
                field_info_names.append("CastroFieldInfo")
            else: 
                field_info_names.append("BoxlibFieldInfo")

    for dset_name, fi_name in zip(dataset_names, field_info_names):
        fi = getattr(this_f, fi_name)
        nfields = 0
        if hasattr(fi, "known_other_fields"):
            known_other_fields = fi.known_other_fields
            nfields += len(known_other_fields)
        else:
            known_other_fields = []
        if hasattr(fi, "known_particle_fields"):
            known_particle_fields = fi.known_particle_fields
            if 'Tipsy' in fi_name:
                known_particle_fields += tuple(fi.aux_particle_fields.values())
            nfields += len(known_particle_fields)
        else:
            known_particle_fields = []
        if nfields > 0:
            print ".. _%s_specific_fields:\n" % dset_name.replace("Dataset", "")
            h = "%s-Specific Fields" % dset_name.replace("Dataset", "")
            print h
            print "-" * len(h) + "\n"
            for field in known_other_fields:
                print_frontend_field(frontend, field, False)
            for field in known_particle_fields:
                if frontend in ["sph", "halo_catalogs", "sdf"]:
                    print_frontend_field("particle_type", field, True)
                else:
                    print_frontend_field("io", field, True)

print footer
