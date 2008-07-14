.. index::
   single: GUI
   module: reason
   module: yt.reason

Using Reason
============

Reason is the GUI for yt.  It's developed in wxPython, and all attempts have
been made to expose a susbtantial API to enable interactive plot addition, as
well as a command-line for examining raw data.

Starting Reason
---------------

If you've downloaded the binary package for OSX, you can simply double-click
the icon.  Otherwise, if you have `wxPython <http://www.wxpython.org/>`_
installed, and yt is installed globally, you should be able to simply type
``reason`` to initiate it.

When reason starts up, it will execute all the code in your plugin file, so any
fields you've defined there will be defined in the dialog boxes.

Opening a Dataset
-----------------

.. index:: projection, slice

Making Slices and Projections
-----------------------------

.. index:: sphere, EnzoSphere

Obtaining Spheres
-----------------

Making Cutting Planes
---------------------

.. index:: phase plot, cutting plane

Making Phase Plots
------------------

Using the Shell
---------------

