import inspect
from yt.mods import *


def islambda(f):
    return inspect.isfunction(f) and \
           f.__name__ == (lambda: True).__name__

header = r"""
.. _field-list:

Field List
==========

This is a list of all fields available in ``yt``.  It has been organized by the
type of code that each field is supported by.  "Universal" fields are available
everywhere, "Enzo" fields in Enzo datasets, "Orion" fields in Orion datasets,
and so on.

Try using the ``pf.field_list`` and ``pf.derived_field_list`` to view the
native and derived fields available for your dataset respectively. For example
to display the native fields in alphabetical order:

.. notebook-cell::

  from yt.mods import *
  pf = load("Enzo_64/DD0043/data0043")
  for i in sorted(pf.field_list):
    print i

.. note:: Universal fields will be overridden by a code-specific field.

.. rubric:: Table of Contents

.. contents::
   :depth: 2
   :local:
   :backlinks: none
"""

print header

seen = []


def print_all_fields(fl):
    for fn in sorted(fl):
        df = fl[fn]
        f = df._function
        cv = df._convert_function
        if [f, cv] in seen:
            continue
        seen.append([f, cv])
        print "%s" % (df.name)
        print "+" * len(df.name)
        print
        if len(df._units) > 0:
            print "   * Units: :math:`%s`" % (df._units)
        if len(df._projected_units) > 0:
            print "   * Projected Units: :math:`%s`" % (df._projected_units)
        print "   * Particle Type: %s" % (df.particle_type)
        print
        print "**Field Source**"
        print
        if islambda(f):
            print "No source available."
            print
            continue
        else:
            print ".. code-block:: python"
            print
            for line in inspect.getsource(f).split("\n"):
                print "  " + line
            print
        print "**Convert Function Source**"
        print
        if islambda(cv):
            print "No source available."
            print
            continue
        else:
            print ".. code-block:: python"
            print
            for line in inspect.getsource(cv).split("\n"):
                print "  " + line
            print


print "Universal Field List"
print "--------------------"
print
print_all_fields(FieldInfo)

print "Enzo-Specific Field List"
print "------------------------"
print
print_all_fields(EnzoFieldInfo)

print "Orion-Specific Field List"
print "-------------------------"
print
print_all_fields(OrionFieldInfo)

print "FLASH-Specific Field List"
print "-------------------------"
print
print_all_fields(FLASHFieldInfo)

print "Athena-Specific Field List"
print "--------------------------"
print
print_all_fields(AthenaFieldInfo)

print "Nyx-Specific Field List"
print "-----------------------"
print
print_all_fields(NyxFieldInfo)

print "Chombo-Specific Field List"
print "--------------------------"
print
print_all_fields(ChomboFieldInfo)

print "Pluto-Specific Field List"
print "--------------------------"
print
print_all_fields(PlutoFieldInfo)

print "Grid-Data-Format-Specific Field List"
print "------------------------------------"
print
print_all_fields(GDFFieldInfo)

print "Generic-Format (Stream) Field List"
print "----------------------------------"
print
print_all_fields(StreamFieldInfo)
