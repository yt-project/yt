.. _basic_objects:

Basic Objects
=============

Derived Fields and Indices
--------------------------

Let's try out some of the derived fields.::

   >>> print my_grid["Entropy"]

Entropy is a derived field.  But how does it work?  Well, we have a built-in
source-code inspector.  ::

   >>> print lagos.fieldInfo["Entropy"].get_source()

But what if we only want a couple values?  We can directly grab indices, in C
order.  ::

   >>> print my_grid["Density"][0,10,4] # 0-indexed!

We can also slice, using the slice operator: ':' ::

   >>> print my_grid["Density"][0,:,4]
   >>> print my_grid["Density"][0,2:4,4]

The first will return to us everything where the x index is zero and the z
index is four.  The second will return the same, except restrict the return to
only between indices 2 and 4.

Okay, so now you've seen how to grab fields from a grid.  The thing about
getting fields from a grid, though, is that it works the exact same way as
all the other data types!  (Well, maybe not projection.)  So that raises the
question, how do we get more data objects?  Well, the hierarchy is our friend
here.  We can get slices, projections, spheres, rectangular prisms,
fixed-resolution grids and cutting planes from the herarchy.

.. index:: Slices, 2ddata;Slices

Slices
------

Let's start with a slice.  A slice needs some info, but how do we know what
type?  Well, let's use 'help'!::

   >>> help(a.h.slice)

Looks like we need to feed it an axis, a coordinate, and maybe some fields
and a center.  Axes are always done as numbers, with 0 = x, 1 = y, 2 = z.::

   >>> my_slice = a.h.slice(0, 0.5, "Density")  # We don't really need to pass the field!

Now let's test out getting new fields.::

   >>> print my_slice["Pressure"], my_slice["Density"]

Okay, so enough of the slices.  (Although they can do some more stuff, like
shift around inside the domain, and so on and so forth.)  Slices are typically
used implicity with plotting (more on that later) and not very often from
the interpreter.

Boxes
-----

Let's do a region.  A region is a box with a left edge and a right edge, and a
center.  (You need a center for some other calculations, not for getting the
region.)  Let's get started by grabbing the ENTIRE domain.::

   >>> my_region = a.h.region([0.5,0.5,0.5], [0.0,0.0,0.0], [1.0,1.0,1.0])

This data object will respect a fixed protocol, as will all the 3D objects.
You ask it for a field, and it returns that field: ::

   >>> my_region["Density"]

We'll talk more about which fields are available in the next section.
