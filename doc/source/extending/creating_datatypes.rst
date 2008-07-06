Creating 3D Datatypes
=====================

The three-dimensional datatypes in yt follow a fairly simple protocol.  The
basic principle is that if you want to define a region in space, that region
must be identifable from some sort of cut applied against the cells --
typically, in yt, this is done by examining the geomery.  (The
:class:`yt.lagos.ExtractedRegionBase` type is a notable exception to this, as
it is defined as a subset of an existing data object.)

In principle, you can define any number of 3D data objects, as long as the
following methods are implemented to protocol specifications.

.. function:: __init__(self, args, kwargs)

   This function can accept any number of arguments but must eventually call
   Enzo3DData.__init__.  It is used to set up the various parameters that
   define the object.

.. function:: _get_list_of_grids(self)

   This function must set the property _grids to be a list of the grids
   that should be considered to be a part of the data object.  Each of these
   will be partly or completely contained within the object.

.. function:: _is_fully_enclosed(self, grid)

   This function returns true if the entire grid is part of the data object
   and false if it is only partly enclosed.

.. function:: _get_cut_mask(self, grid)

   This function returns a boolean mask in the shape of the grid.  All of the
   cells set to 'True' will be included in the data object and all of those set
   to 'False' will be excluded.  Typically this is done via some logical
   operation.

For a good example of how to do this, see the
:class:`yt.lagos.EnzoCylinderBase` source code.
