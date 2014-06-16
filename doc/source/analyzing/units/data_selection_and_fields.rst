.. _data_selection_and_fields:

Data selection and fields
=========================

.. notebook:: 2)_Data_Selection_and_fields.ipynb

Derived Fields
--------------

.. This needs to be added outside the notebook since user-defined derived fields
   require a 'fresh' kernel.

.. warning:: Note: derived field definitions need to happen *before* a dataset
             is loaded.  This means changes to the following cells will only be
             picked up on a fresh kernel.  Select Kernel -> Restart on the
             IPython menu bar to restart the kernel.

New derived fields can be added just like in old vesions of yt.  The most
straightforward way to do this is to apply the `derived_field` decorator on a
function that defines a field.

The following example creates a derived field for the square root of the cell
volume.

.. notebook-cell::

   from yt.mods import *
   import numpy as np

   @derived_field(name='root_cell_volume', units='cm**(3/2)')
   def root_cell_volume(field, data):
     return np.sqrt(data['cell_volume'])

   ds = load('HiresIsolatedGalaxy/DD0044/DD0044')

   dd = ds.all_data()
   dd['root_cell_volume']

No special unit logic needs to happen inside of the function - `np.sqrt` will
convert the units of the `density` field appropriately:

.. notebook-cell::
   :skip_exceptions:

   from yt.mods import *
   import numpy as np

   ds = load('HiresIsolatedGalaxy/DD0044/DD0044')
   dd = ds.all_data()

   print dd['cell_volume'].in_cgs()
   print np.sqrt(dd['cell_volume'].in_cgs())

That said, it is necessary to specify the units in the call to the
:code:`@derived_field` decorator.  Not only does this ensure the returned units
will be exactly what you expect, it also allows an in-place conversion of units,
just in case the function returns a field with dimensionally equivalent units.

For example, let's redo the above example but ask for units of
:code:`Mpc**(3/2)`:

.. notebook-cell::

   from yt.mods import *

   @derived_field(name='root_cell_volume', units='Mpc**(3/2)')
   def root_cell_volume(field, data):
     return np.sqrt(data['cell_volume'])

   ds = load('HiresIsolatedGalaxy/DD0044/DD0044')

   dd = ds.all_data()
   dd['root_cell_volume']
