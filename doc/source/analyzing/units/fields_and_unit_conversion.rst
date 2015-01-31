.. _data_selection_and_fields:

Fields and Unit Conversion
==========================

.. notebook:: 2)_Fields_and_unit_conversion.ipynb

Derived Fields
--------------

.. This needs to be added outside the notebook since user-defined derived fields
   require a 'fresh' kernel.

The following example creates a derived field for the square root of the cell
volume.

.. notebook-cell::

   import yt
   import numpy as np

   # Function defining the derived field
   def root_cell_volume(field, data):
      return np.sqrt(data['cell_volume'])

   # Load the dataset
   ds = yt.load('HiresIsolatedGalaxy/DD0044/DD0044')

   # Add the field to the dataset, linking to the derived field function and 
   # units of the field
   ds.add_field(("gas", "root_cell_volume"), units="cm**(3/2)", function=root_cell_volume)

   # Access the derived field like any other field
   ad = ds.all_data()
   ad['root_cell_volume']

No special unit logic needs to happen inside of the function - `np.sqrt` will
convert the units of the `density` field appropriately:

.. notebook-cell::
   :skip_exceptions:

   import yt
   import numpy as np

   ds = yt.load('HiresIsolatedGalaxy/DD0044/DD0044')
   ad = ds.all_data()

   print ad['cell_volume'].in_cgs()
   print np.sqrt(ad['cell_volume'].in_cgs())

That said, it is necessary to specify the units in the call to the
:code:`add_field` function.  Not only does this ensure the returned units
will be exactly what you expect, it also allows an in-place conversion of units,
just in case the function returns a field with dimensionally equivalent units.

For example, let's redo the above example but ask for units of
:code:`Mpc**(3/2)`:

.. notebook-cell::

   import yt
   import numpy as np

   def root_cell_volume(field, data):
      return np.sqrt(data['cell_volume'])

   ds = yt.load('HiresIsolatedGalaxy/DD0044/DD0044')

   # Here we set the default units to Mpc^(3/2)
   ds.add_field(("gas", "root_cell_volume"), units="Mpc**(3/2)", function=root_cell_volume)

   ad = ds.all_data()
   ad['root_cell_volume']
