The Quick Guide to yt
=====================

If you're impatient, like me, you probably just want to pull up some data and
take a look at it.  This guide will help you out!

Starting IPython
----------------

If you've used the installation script that comes with yt, you should have an
isolated environment containing Python 2.5, Matplotlib, yt, IPython, and maybe
wxPython.  Be sure to finish up the instructions by *prepending* the
``LD_LIBRARY_PATH``, ``PATH`` and ``PYTHONPATH`` environment variables with the
output of the script.

If you've done that, go ahead and start up our interactive yt environment:

.. code-block:: bash

   $ iyt

It should start you up in an interpreter, and the namespace will be populated
with the stuff you need.  Really, the command ``iyt`` just opens up IPython and
loads up yt, with some special commands available for you.

You're all set, so let's move on to the next step -- actually opening up your
data!

Opening Your Data File
----------------------

You'll need to know the location of the parameter file from the output you want
to look at.  Let's pretend, for the sake of argument, it's
``/home/mturk/data/galaxy1200.dir/galaxy1200`` and that we have all the right
permissions.  So let's open it, and see what the maximum density is.

.. note:: In IPython, you get filename completion!  So hit tab and it'll guess
   at what you want to open.

.. sourcecode:: ipython

  In [1]: pf = EnzoStaticOutput("/home/mturk/data/galaxy1200.dir/galaxy1200")

  In [2]: v, c = pf.h.find_max("Density")

And then in the variable ``v`` we have the value of the most dense cell, and in
``c`` we have the location of that point.

.. _quick_making_plots:

Making Plots
------------

But hey, what good is the data if we can't see it?  So let's make some plots!
First we need to get a :class:`PlotCollectionInteractive` object, and then
we'll add some slices and projections to it.  Note that we use 0, 1, 2 to refer
to 'x', 'y', 'z' axes.

.. sourcecode:: ipython

  In [3]: pc = PlotCollectionInteractive(pf)
  In [4]: pc.add_slice("Temperature", 0)
  yt.raven   INFO       2008-10-25 11:42:58,429 Added slice of Temperature at x = 0.953125 with 'center' = [0.953125, 0.8046875, 0.6171875]
  Out[4]: <yt.raven.PlotTypes.SlicePlot instance at 0x9882cec>

  In [5]: pc.add_slice("Density", 0)
  yt.raven   INFO       2008-10-25 11:43:45,608 Added slice of Density at x = 0.953125 with 'center' = [0.953125, 0.8046875, 0.6171875]
  Out[5]: <yt.raven.PlotTypes.SlicePlot instance at 0xab83eec>


A window should now pop up for each of these plots.  One will be a line
integral through the simulation, and the other will be a slice.  (If you had
used the :class:`PlotCollection` object, they'd be created off-screen -- this
is the right way to make plots programmatically in scripts.)

We can also adjust the width of the plots very easily:

.. sourcecode:: ipython

  In [6]: pc.set_width(100, 'kpc')

The center is set to the most dense location by default.  (For more
information, see the documentation for :class:`PlotCollection`.)

Saving Plots
------------

Even though the windows are open, we can save these to the file system at high
resolution.

.. sourcecode:: ipython

  In [7]: pc.save("hi")
  Out[7]: ['hi_Slice_x_Temperature.png', 'hi_Slice_x_Density.png']

And that's it!  The plots get saved out, and it returns to you a list of their
filenames.

.. note:: The *save* command will add some data to the end of the filename --
   this helps to keep track of what each saved file is.

A Few More Plots
----------------

You can also add profiles -- radial or otherwise -- and phase diagrams very
easily.

.. sourcecode:: ipython

  In [8]: pc.add_profile_sphere(100.0, 'kpc', ["Density", "Temperature"])
  Out[8]: <yt.raven.PlotTypes.Profile1DPlot instance at 0xada03ec>

  In [9]: pc.add_phase_sphere(10.0, 'pc', ['Density', 'Temperature', 
     ...:                                  'H2I_Fraction'])
  Out[9]: <yt.raven.PlotTypes.PhasePlot instance at 0xada91ef>

Note that the phase plots default to showing a weighted-average in each bin --
weighted by the cell mass in solar masses.  If you want to see a distribution
of mass, you'll need to specify you don't want an average:

.. code-block:: python

  In [10]: pc.add_phase_sphere(10.0, 'pc', ['Density', 'Temperature', 
      ...:                                   'CellMassMsun'], weight=None)

  Out[10]: <yt.raven.PlotTypes.PhasePlot instance at 0xada91ef>
