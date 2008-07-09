Derived Fields
==============

Okay, now we have a region.  We can do stuff with this -- in fact, it has the
same data-access protocol as a slice or a grid!  Plus, math on these things
is pretty fast:::

   >>> electron_fraction = my_region["Electron_Density"]/my_region["Density"]

But the thing is, that one happens to be a derived field already:::

   >>> electron_fraction = my_region["Electron_Fraction"]

In fact, every "_Density" field is automatically set up as a "_Fraction" field.
But what other fields are available to us?::

   >>> print lagos.fieldInfo.keys()

Defining Simple Fields
----------------------

There are **lots**.  Not only that, but you can set up your own fields, too,
which then become available to all data-types. Here's a very simple example,
where we just cube Density.::

   >>> def _DensityCubed(field, data):
   ...     return data["Density"]**3.0
   ...
   >>> def _ConvertDensityCubed(data):
   ...     return 1.0
   ...
   >>> # Now we have our field functions, and we just need to tell the fieldinfo
   >>> # object about it.
   >>> lagos.add_field("DensityCubed", function=_DensityCubed,
   ...     convert_function= _ConvertDensityCubed, units=r"\rm{g}^3/\rm{cm}^9")

This last part might need some explanation.  Derived fields are composed of
three main parts.  The first is the function that defines the field.  It always
accepts two variables -- field, and data.  'Field' is the DerivedField and
'data' is the data object.  Note that we stick to a pretty strict means of
returning fields -- *everything is returned in cgs*!  This means that when we
request "Density" from 'data' it is returned to us *already* in CGS, so we don't
need to convert the units. This may not always be the case for derived fields,
so the second function defines the necessary unit conversion.

More Complicated Fields
-----------------------

The last part, the call to add_field, is the most crucial.  It creates a
DerivedField object which will control the behavior of the derived field in a
couple crucial ways -- what you see here includes the string used for the
units, the string used for the name of the field, the functions that govern the
creation and conversion of the fiels.  However, you can also use 'validators'
-- in fact, for some fields, you *should* use validators.  Validators are the
most powerful part of the derived field process, because they can allow you to
require ghost zones, which will be automatically generated.  Check the
"AveragedDensity" and "DivV" fields in the source code for an example of this.
