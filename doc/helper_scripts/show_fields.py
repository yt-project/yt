import inspect

import numpy as np

import yt.frontends as frontends_module
from yt.config import ytcfg
from yt.fields.derived_field import NullFunc
from yt.frontends.api import _frontends
from yt.frontends.stream.fields import StreamFieldInfo
from yt.funcs import obj_length
from yt.testing import fake_random_ds
from yt.units import dimensions
from yt.units.yt_array import Unit
from yt.utilities.cosmology import Cosmology

fields, units = [], []

for fname, (code_units, _aliases, _dn) in StreamFieldInfo.known_other_fields:
    fields.append(("gas", fname))
    units.append(code_units)
base_ds = fake_random_ds(4, fields=fields, units=units)
base_ds.index
base_ds.cosmological_simulation = 1
base_ds.cosmology = Cosmology()


ytcfg["yt", "internals", "within_testing"] = True
np.seterr(all="ignore")


def _strip_ftype(field):
    if not isinstance(field, tuple):
        return field
    elif field[0] == "all":
        return field
    return field[1]


np.random.seed(int(0x4D3D3D3))
units = [base_ds._get_field_info(*f).units for f in fields]
fields = [_strip_ftype(f) for f in fields]
ds = fake_random_ds(16, fields=fields, units=units, particles=1)
ds.parameters["HydroMethod"] = "streaming"
ds.parameters["EOSType"] = 1.0
ds.parameters["EOSSoundSpeed"] = 1.0
ds.conversion_factors["Time"] = 1.0
ds.conversion_factors.update({f: 1.0 for f in fields})
ds.gamma = 5.0 / 3.0
ds.current_redshift = 0.0001
ds.cosmological_simulation = 1
ds.hubble_constant = 0.7
ds.omega_matter = 0.27
ds.omega_lambda = 0.73
ds.cosmology = Cosmology(
    hubble_constant=ds.hubble_constant,
    omega_matter=ds.omega_matter,
    omega_lambda=ds.omega_lambda,
    unit_registry=ds.unit_registry,
)
for my_unit in ["m", "pc", "AU", "au"]:
    new_unit = f"{my_unit}cm"
    my_u = Unit(my_unit, registry=ds.unit_registry)
    ds.unit_registry.add(
        new_unit,
        my_u.base_value,
        dimensions.length,
        "\\rm{%s}/(1+z)" % my_unit,
        prefixable=True,
    )


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
    print(i)

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
print(header)

seen = []


def fix_units(units, in_cgs=False):
    unit_object = Unit(units, registry=ds.unit_registry)
    if in_cgs:
        unit_object = unit_object.get_cgs_equivalent()
    latex = unit_object.latex_representation()
    return latex.replace(r"\ ", "~")


def print_all_fields(fl):
    for fn in sorted(fl):
        df = fl[fn]
        f = df._function
        s = f"{df.name}"
        print(s)
        print("^" * len(s))
        print()
        if obj_length(df.units) > 0:
            # Most universal fields are in CGS except for these special fields
            if df.name[1] in [
                "particle_position",
                "particle_position_x",
                "particle_position_y",
                "particle_position_z",
                "entropy",
                "kT",
                "metallicity",
                "dx",
                "dy",
                "dz",
                "cell_volume",
                "x",
                "y",
                "z",
            ]:
                print(f"   * Units: :math:`{fix_units(df.units)}`")
            else:
                print(f"   * Units: :math:`{fix_units(df.units, in_cgs=True)}`")
        print(f"   * Sampling Method: {df.sampling_type}")
        print()
        print("**Field Source**")
        print()
        if f == NullFunc:
            print("No source available.")
            print()
            continue
        else:
            print(".. code-block:: python")
            print()
            for line in inspect.getsource(f).split("\n"):
                print("  " + line)
            print()


ds.index
print_all_fields(ds.field_info)


class FieldInfo:
    """a simple container to hold the information about fields"""

    def __init__(self, ftype, field, ptype):
        name = field[0]
        self.units = ""
        u = field[1][0]
        if len(u) > 0:
            self.units = r":math:`\mathrm{%s}`" % fix_units(u)
        a = [f"``{f}``" for f in field[1][1] if f]
        self.aliases = " ".join(a)
        self.dname = ""
        if field[1][2] is not None:
            self.dname = f":math:`{field[1][2]}`"

        if ftype != "particle_type":
            ftype = f"'{ftype}'"
        self.name = f"({ftype}, '{name}')"
        self.ptype = ptype


current_frontends = [f for f in _frontends if f not in ["stream"]]

for frontend in current_frontends:
    this_f = getattr(frontends_module, frontend)
    field_info_names = [fi for fi in dir(this_f) if "FieldInfo" in fi]
    dataset_names = [dset for dset in dir(this_f) if "Dataset" in dset]

    if frontend == "gadget":
        # Drop duplicate entry for GadgetHDF5, add special case for FieldInfo
        # entry
        dataset_names = ["GadgetDataset"]
        field_info_names = ["GadgetFieldInfo"]
    elif frontend == "boxlib":
        field_info_names = []
        for d in dataset_names:
            if "Maestro" in d:
                field_info_names.append("MaestroFieldInfo")
            elif "Castro" in d:
                field_info_names.append("CastroFieldInfo")
            else:
                field_info_names.append("BoxlibFieldInfo")
    elif frontend == "chombo":
        # remove low dimensional field info containers for ChomboPIC
        field_info_names = [
            f for f in field_info_names if "1D" not in f and "2D" not in f
        ]

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
            if "Tipsy" in fi_name:
                known_particle_fields += tuple(fi.aux_particle_fields.values())
            nfields += len(known_particle_fields)
        else:
            known_particle_fields = []
        if nfields > 0:
            print(f".. _{dset_name.replace('Dataset', '')}_specific_fields:\n")
            h = f"{dset_name.replace('Dataset', '')}-Specific Fields"
            print(h)
            print("-" * len(h) + "\n")

            field_stuff = []
            for field in known_other_fields:
                field_stuff.append(FieldInfo(frontend, field, False))
            for field in known_particle_fields:
                if frontend in ["sph", "halo_catalogs", "sdf"]:
                    field_stuff.append(FieldInfo("particle_type", field, True))
                else:
                    field_stuff.append(FieldInfo("io", field, True))

            # output
            len_name = 10
            len_units = 5
            len_aliases = 7
            len_part = 9
            len_disp = 12
            for f in field_stuff:
                len_name = max(len_name, len(f.name))
                len_aliases = max(len_aliases, len(f.aliases))
                len_units = max(len_units, len(f.units))
                len_disp = max(len_disp, len(f.dname))

            fstr = "{nm:{nw}}  {un:{uw}}  {al:{aw}}  {pt:{pw}}  {dp:{dw}}"
            header = fstr.format(
                nm="field name",
                nw=len_name,
                un="units",
                uw=len_units,
                al="aliases",
                aw=len_aliases,
                pt="particle?",
                pw=len_part,
                dp="display name",
                dw=len_disp,
            )

            div = fstr.format(
                nm="=" * len_name,
                nw=len_name,
                un="=" * len_units,
                uw=len_units,
                al="=" * len_aliases,
                aw=len_aliases,
                pt="=" * len_part,
                pw=len_part,
                dp="=" * len_disp,
                dw=len_disp,
            )
            print(div)
            print(header)
            print(div)

            for f in field_stuff:
                print(
                    fstr.format(
                        nm=f.name,
                        nw=len_name,
                        un=f.units,
                        uw=len_units,
                        al=f.aliases,
                        aw=len_aliases,
                        pt=f.ptype,
                        pw=len_part,
                        dp=f.dname,
                        dw=len_disp,
                    )
                )

            print(div)
            print("")

print(footer)
