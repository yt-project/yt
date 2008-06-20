Basic Analysis Tasks
====================

Because yt provides you with an low-cost entry-point for examining both
particle and baryon data, it can be used for very quickly-written analysis
tasks.  To that end, some of the simplest ones are just looking at the
distribution of mass in a given simulation.

Summing Up Mass
---------------

The first, most basic thing to try is to sum up all the mass in a given region.
Here we simply sum up over the entire domain.
(``cookbook_mass_sum.py``)

.. literalinclude:: ../../../examples/cookbook_mass_sum.py
   :language: python
   :linenos:

Summing Up Halo Masses
----------------------

With the inclusion of HOP in the code, we can now identify halos and use those
as inputs into our calculations of things like the total mass in a region.  For
instance, one can iterate over a set of HOP-centers and examine (as above) the
mass as divided into stars, dark matter, and baryons.
(``cookbook_hop_mass_sum.py``)

.. literalinclude:: ../../../examples/cookbook_hop_mass_sum.py
   :language: python
   :linenos:

