Constructing Data Objects
-------------------------

These recipes demonstrate a few uncommon methods of constructing data objects
from a simulation.

.. _cookbook-find_clumps:

Identifying Clumps
~~~~~~~~~~~~~~~~~~

This is a recipe to show how to find topologically connected sets of cells
inside a dataset.  It returns these clumps and they can be inspected or
visualized as would any other data object.  More detail on this method can be
found in `astro-ph/0806.1653`.

.. yt_cookbook:: find_clumps.py

Boolean Data Objects
~~~~~~~~~~~~~~~~~~~~

Below shows the creation of a number of boolean data objects, which are built
upon previously-defined data objects. The boolean data objects can be used like
any other, except for a few cases.  Please see :ref:`boolean_data_objects` for
more information.

.. yt_cookbook:: boolean_data_objects.py

Extracting Fixed Resolution Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a recipe to show how to open a dataset and extract it to a file at a
fixed resolution with no interpolation or smoothing.  Additionally, this recipe
shows how to insert a dataset into an external HDF5 file using h5py.

.. yt_cookbook:: extract_fixed_resolution_data.py
