Reduced Objects
===============

Once a physically meaningful collection of data has been selected, often
analysis needs to be conducted that reduces that object to some smaller set of
numbers, often numbers that can be plotted or included in a paper.  As such,
``yt`` includes facilities for doing such things!

Derived Quantities
------------------

Derived quantities are the 'simple numbers' you can get out of objects -- this
includes things like (baryon) spin parameter, the angular momentum vector, the
gravitational boundedness of a given physical object and so forth.  More
information, including a list of quantities, can be found in
:ref:`derived_quantities`.  These get called in an admittedly slightly awkward
way, but that is because they can also accept parameters:

.. code-block:: python

   is_bound = my_sphere.quantities["IsBound"]()
   L_vec = my_sphere.quantities["AngularMomentumVector"]()

Some return booleans, others return arrays, and some just return numbers.

Profiles
--------

The idea of a profile is rather more like an n-dimensional histogram with an
arbitrary weighting.  There are only a few things to keep in mind here -- the
profiles will truncate any points that do not fall within their boundaries (if
you do not specify boundaries, they will auto-select them to include all
points) and the value in a given bin is governed by the value of the 'weight.'
If the 'weight' is set to ``None``, then the values in a given bin are simply
added; typically it is set to ``CellMassMsun`` which will multiply each value
by the mass in its cell, sum, and then divide by the total mass in that bin.

Typically 1-D profiles are generated either as 'radial profiles' or
'probability distribution functions' and 2-D profiles are thought of as phase
plots.  A common usage pattern is to generate a phase plot of mass in a bin in
the plane of two quantities, then overplot the 1-D profile with a weighted
average.  In this way, the spread in a quantity as well as its average value
can be displayed simultaneously.

You can generate these profiles yourself, but because their API is completely
exposed by the :class:`~yt.raven.PlotCollection` object, you rarely need to.  But, if you
are determined to, see :ref:`making_profiles`.  Otherwise, check out
:ref:`profiles_and_phase_plots`.

Projections
-----------

Projections are, in a sense, a means of casting rays through the simulational
volume.  
