Profiles, Contours and Halo Finding
===================================

Profiling
---------

Profiling in yt is a means of generating arbitrary histograms -- for instance,
phase diagrams, radial profiles, and even more complicated 3D examinations.

.. currentmodule:: yt.lagos

.. autoclass:: yt.lagos.BinnedProfile1D
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.BinnedProfile2D
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.BinnedProfile3D
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.StoredBinnedProfile3D
   :members:
   :show-inheritance:

Contour Finding
---------------

Typically this is done via the :meth:`extract_connected_sets` on a data object.
However, you can call it manually, as is done in the clump finding scripts.

.. autofunction:: yt.lagos.identify_contours

.. autoclass:: yt.lagos.GridConsiderationQueue
   :members:

Halo Finding
------------

yt now includes the HOP algorithm and implementation from the Enzo source
distribution, with some modifications by both Stephen Skory and myself.

.. autoclass:: yt.lagos.hop.HopList
   :members:

.. autoclass:: yt.lagos.hop.HopGroup
   :members:
