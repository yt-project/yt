Plotting Profiles
=================

One-Dimensional Profiles
------------------------

This example first plots a very simple one-dimensional profile  of density and
temperature.  Note that the PlotCollection object deals with finding the
location of the most dense point, and the profile is automatically centered
there.  The next stage of the script is a bit more complicated, as it extracts
a region and uses that as an input to another diagram.  Finally, we save it.
(``cookbook_plotting_1dprofiles.py``)

.. literalinclude:: ../../../examples/cookbook_plotting_1dprofiles.py
   :language: python
   :linenos:

Two-Dimensional Profiles
------------------------

This example is almost identical to the one above, except that we add an
additional field to the specification.  Note that in the first call to
:meth:`add_phase_sphere` we don't specify a weight; it defaults to displaying
the ``CellMassMsun``-weighted average in each bin; in the second example, we
specify ``weight = None`` which means that it will do no averaging, and instead
simply plot the total (sum) in each bin.  This allows us to see
mass-distribution.  (``cookbook_plotting_2dprofiles.py``)

.. literalinclude:: ../../../examples/cookbook_plotting_2dprofiles.py
   :language: python
   :linenos:
