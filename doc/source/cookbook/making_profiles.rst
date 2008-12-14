.. _making_profiles:

Making Profiles
===============

Making 1D Profiles
------------------

This is an example of one of the most common tasks to perform on a simulation
-- a radial profile.  I would like to note, however, that the binning field
(that is set to "Radius" here) can be set to anything.

Additionally, the 'accumulation' parameter can be useful in one-dimensional
profiles, as it sums up from low-to-high along a given axis.  Here we use it to
construct a mass-enclosed field.

.. literalinclude:: ../../../examples/cookbook_making_1d_profiles.py
   :language: python
   :linenos:


Making 2D Profiles
------------------

In this next example we construct a two-dimensional, weighted-average binning
of a sphere centered at the point of maximum density in an output.  This sets
up the basic :class:`BinnedProfile2D` object, which we then add fields to.
Note that we allow the default weighting to be applied to the x-velocity
field, which will give us an average value, weighted by CellMassMsun.  For the
second field we add, we explicitly ask to have the values returned to us
unweighted, which means we are given the distribution of mass along the dimensions
of NumberDensity and Temperature.

.. literalinclude:: ../../../examples/cookbook_making_2d_profiles.py
   :language: python
   :linenos:
