.. _demeshening:

How Particles are Indexed
=========================

With yt-4.0, the method by which particles are indexed changed considerably.
Whereas in previous versions, particles were indexed based on their position in
an octree (the structure of which was determined by particle number density),
in yt-4.0 this system was overhauled to utilize a `bitmap
index<https://en.wikipedia.org/wiki/Bitmap_index>`_ based on a space-filling
curve, using a `enhanced word-aligned
hybrid<https://github.com/lemire/ewahboolarray>` boolean array as their
backend.

.. note::

   You may see scattered references to "the demeshening" (including in the
   filename for this document!)  This was a humorous name used in the yt
   development process to refer to removing a global (octree) mesh for
   particle codes.

Effectively, what this does is allow yt to treat the particles as discrete
objects (or with an area of influence) and use their positions in a multi-level
index to optimize and minimize the disk operations necessary to load only those
particles it needs.

This is described in some detail in the `yt 4.0
paper<https://yt-project.github.io/yt-4.0-paper/>`_ in the section entitled
`Indexing Discrete-Point
Datasets<https://yt-project.github.io/yt-4.0-paper/#sec:point_indexing>`_.

In brief, however, what this relies on is two numbers, ``index_order1`` and
``index_order2``.  These control the "coarse" and "refined" sets of indices,
and they are supplied to any particle dataset ``load()`` in the form of a tuple
as the argument ``index_order``.  By default these are set to 5 and 7,
respectively, but it is entirely likely that a different set of values will
work better for your purposes.

For example, if you were to use the sample Gadget-3 dataset, you could override
the default values and use values of 5 and 5 by specifying this argument to the
``load_sample`` function; this works with ``load`` as well.

.. code-block:: python

   ds = yt.load_sample("Gadget3-snap-format2", index_order=(5, 5))

Index Order and Efficiency
--------------------------



Index Caching
-------------
