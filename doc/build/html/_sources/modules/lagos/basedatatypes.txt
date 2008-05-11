Data Containers
===============

yt provides a number of data containers, defined such that they satisfy a
logical need.  Each of these provides only the finest-resolution cells, unless
an option is available to restrict the levels from which they draw, as is the
case with :class:`~yt.lagos.EnzoCoveringGrid`.

Each of these implements the same primary protocol - all data values can be
access dictionary-style:

.. code-block:: python

   >>> object["Density"]
   >>> object["Density"].max()

And so on.  For more information, see the
:ref:`the tutorial topic <basic_objects>` on basic objects.

.. currentmodule:: yt.lagos

Base Containers
---------------

.. autoclass:: yt.lagos.FakeGridForParticles
   :members:

.. autoclass:: yt.lagos.EnzoData
   :members:

.. autoclass:: yt.lagos.Enzo1DData
   :members:
   :inherited-members:

.. autoclass:: yt.lagos.Enzo2DData
   :members:
   :inherited-members:

.. autoclass:: yt.lagos.Enzo3DData
   :members:
   :inherited-members:

1D Data Containers
------------------

Each of these inherits from :class:`~yt.lagos.Enzo1DData`, and so has all the
member functions defined there.

.. autoclass:: yt.lagos.EnzoOrthoRayBase
   :members:

2D Data Containers
------------------

.. autoclass:: yt.lagos.EnzoSliceBase
   :members:

.. autoclass:: yt.lagos.EnzoCuttingPlaneBase
   :members:

.. autoclass:: yt.lagos.EnzoProjBase
   :members:

3D Data Containers
------------------

.. autoclass:: yt.lagos.EnzoSphereBase
   :members:

.. autoclass:: yt.lagos.EnzoRegionBase
   :members:

.. autoclass:: yt.lagos.EnzoCylinderBase
   :members:

.. autoclass:: yt.lagos.EnzoGridCollection
   :members:

.. autoclass:: yt.lagos.EnzoCoveringGrid
   :members:

.. autoclass:: yt.lagos.ExtractedRegionBase
   :members:

