The Quick Guide to yt
=====================

If you're impatient, like me, you probably just want to pull up some data and
take a look at it.  This guide will help you out!

Starting Python
---------------

If you've used the installation script that comes with yt, you should
have an isolated environment containing Python 2.5, Matplotlib, wxPython, and
yt.  Be sure to finish up the instructions by *prepending* the
``LD_LIBRARY_PATH``, ``PATH`` and ``PYTHONPATH`` environment variables with the
output of the script!

If you've done that, go ahead and start up yt:

.. code-block:: bash

   $ yt

It should start you up in an interpreter, and the namespace will be populated
with the stuff you need.  Really, the command ``yt`` just opens up Python and
loads up yt -- nothing too fancy!

You're all set, so let's move on to the next step -- actually opening up your
data!

Opening Your Data File
----------------------

You'll need to know the location of the parameter file from the output you want
to look at.  Let's pretend, for the sake of argument, it's
``/scratch/mturk/DataDump0010.dir/DataDump0010`` and that we have all the right
permissions.  So let's open it, and see what the maximum density is.

.. code-block:: python

   >>> pf = EnzoStaticOutput("/scratch/mturk/DataDump0010.dir/DataDump0010")
   >>> v, c = pf.h.find_max("Density")

And then in the variable ``v`` we have the value of the most dense cell, and in
``c`` we have the location of that point.

Making Plots
------------

But hey, what good is the data if we can't see it?  So let's make some plots!
First we need to get a :class:`PlotCollection` object, and then we'll add some
slices and projections to it.  Note that we use 0, 1, 2 to refer to 'x', 'y', 'z'
axes.

.. code-block:: python

   >>> pc = PlotCollection(pf)
   >>> pc.add_slice("Temperature", 0)
   >>> pc.add_projection("Density", 2)

It makes these plots all off-screen.  (If you had used the
:class:`PlotCollectionInteractive` object, they'd be there, displayed, as soon
as you added them.)

We can also adjust the width of the plots very easily:

.. code-block:: python

   >>> pc.set_width(100, 'kpc')

The center is set to the most dense location by default.  (For more
information, see the documentation for :class:`PlotCollection`.)

Saving Plots
------------

Because all of these plots are off-screen, we save to the file system before we
can see them.

.. code-block:: python

   >>> pc.save("hi")

And that's it!  The plots get saved out, and it returns to you a list of their
filenames.

A Few More Plots
----------------

You can also add profiles -- radial or otherwise -- and phase diagrams very
easily.

.. code-block:: python

   >>> pc.add_profile_sphere(100.0, 'kpc', ['Density', 'Temperature'])
   >>> pc.add_phase_sphere(10.0, 'pc', ['Density', 'Temperature', 
   ...                                  'H2I_Fraction'])

But again, you have to save these out before you can view them.  Note that the
phase plots default to showing a weighted-average in each bin -- weighted by
the cell mass in solar masses.  If you want to see a distribution of mass,
you'll need to specify you don't want an average:

.. code-block:: python

   >>> pc.add_phase_sphere(10.0, 'pc', ['Density', 'Temperature', 
   ...                                  'CellMassMsun'], weight=None)

