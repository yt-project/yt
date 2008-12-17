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

.. currentmodule:: yt.lagos

Base Classes
------------

.. autoclass:: yt.lagos.EnzoData
   :members:

.. autoclass:: yt.lagos.Enzo1DData
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.Enzo2DData
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.Enzo3DData
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: yt.lagos.FakeGridForParticles
   :members:

.. _data_objects:

Physical Objects
================

1D Data Containers
------------------

Each of these inherits from :class:`~yt.lagos.Enzo1DData`, and so has all the
member functions defined there.

.. autoclass:: yt.lagos.EnzoOrthoRayBase
   :members:
   :show-inheritance:

2D Data Containers
------------------

.. autoclass:: yt.lagos.EnzoSliceBase
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoCuttingPlaneBase
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoProjBase
   :members:
   :show-inheritance:

3D Data Containers
------------------

.. autoclass:: yt.lagos.EnzoSphereBase
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoRegionBase
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoCylinderBase
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoGridCollection
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoCoveringGrid
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.EnzoSmoothedCoveringGrid
   :members:
   :show-inheritance:

.. autoclass:: yt.lagos.ExtractedRegionBase
   :members:
   :show-inheritance:

.. _derived_quantities:

Derived Quantities
------------------

All of these are accessed via the ``.quantities[]`` object, and feeding it the
function name without the leading underscore.  For instance:

.. code-block:: python

   my_sphere.quantities["TotalMass"]()

They all accept the ``lazy_reader`` option, which governs whether the
calculation is performed out of core or not.

.. currentmodule:: yt.lagos.DerivedQuantities

.. autofunction:: yt.lagos.DerivedQuantities._TotalMass

.. autofunction:: yt.lagos.DerivedQuantities._CenterOfMass

.. autofunction:: yt.lagos.DerivedQuantities._WeightedAverageQuantity

.. autofunction:: yt.lagos.DerivedQuantities._BulkVelocity

.. autofunction:: yt.lagos.DerivedQuantities._AngularMomentumVector

.. autofunction:: yt.lagos.DerivedQuantities._BaryonSpinParameter

.. autofunction:: yt.lagos.DerivedQuantities._IsBound

.. autofunction:: yt.lagos.DerivedQuantities._Extrema

