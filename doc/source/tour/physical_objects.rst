Physical Objects
================

Assuming the simulation ran correctly, and debugging doesn't need to occur, the
primary means of interaction with simulation output should be through
physically meaningful objects.  Specifically, rather than viewing a collection
of points in a single grid patch, one should be able to obtain and manipulate
collections of points inside spheres, disks, 'clumps' and so forth.  To
facilitate that, ``yt`` provides means of obtaining objects.

Creating Objects
----------------

These objects are all created from the hierarchy object; for instance, to
obtain a sphere centered at ``(0.3, 0.5, 0.8)`` with a radius of one
megaparsec, you could do:

.. code-block:: python

   sp = pf.h.sphere([0.3,0.5,0.8], 1.0/pf['mpc'])

The variable ``sp`` is now an object that can be examined, manipulated and
analyzed.  There are several different objects available -- for more
information, see :ref:`data_objects`.

Means of Interactions
---------------------

These objects have several means of interactions, but the most primitive is
simply to obtain all of the data that is contained in a given object.

.. code-block:: python

   print sp["CellMassMsun"].sum()

This accesses every cell located within the sphere, calculates its mass in
solar masses, and prints the sum.  The actual array of data is a NumPy array;
these have a number of methods, of which ``sum`` is but one -- for more
information, see the NumPy documentation.  (Additionally, you can just call
``help`` on the array: ``help(sp["CellMassMsun"])`` .)
