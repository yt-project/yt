Constructing Data Objects
-------------------------

These recipes demonstrate a few uncommon methods of constructing data objects
from a simulation.

Creating Particle Filters
~~~~~~~~~~~~~~~~~~~~~~~~~

Create particle filters based on the age of star particles in an isolated
disk galaxy simulation.  Determine the total mass of each stellar age bin
in the simulation.  Generate projections for each of the stellar age bins.

.. yt_cookbook:: particle_filter.py

.. _cookbook-find_clumps:

Identifying Clumps
~~~~~~~~~~~~~~~~~~

This is a recipe to show how to find topologically connected sets of cells
inside a dataset.  It returns these clumps and they can be inspected or
visualized as would any other data object.  More detail on this method can be
found in `Smith et al. 2009
<https://ui.adsabs.harvard.edu/abs/2009ApJ...691..441S>`_.

.. yt_cookbook:: find_clumps.py

.. _extract_frb:

Extracting Fixed Resolution Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a recipe to show how to open a dataset and extract it to a file at a
fixed resolution with no interpolation or smoothing.  Additionally, this recipe
shows how to insert a dataset into an external HDF5 file using h5py.

.. yt_cookbook:: extract_fixed_resolution_data.py
