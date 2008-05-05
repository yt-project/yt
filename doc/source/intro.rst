Introduction
============

History
-------

My name is Matthew Turk, and I am the primary author of yt.  I designed and
implemented it during the course of my graduate studies working with Prof. Tom
Abel at Stanford University, under the auspices of the Department of Energy
through the Stanford Linear Accelerator Center and Los Alamos National Lab.  It
has evolved from a simple data-reader and exporter into what I believe is a
fully-featured toolkit.

yt was designed to be a completely Free (as in beer *and* as in freedom)
user-extensible framework for analyzing and visualizing adaptive mesh
refinement data.  It relies on no proprietary software -- although it can be
and has been extended to interface with proprietary software and libraries --
and has been designed from the ground up to enable users to be as immersed in
the data as they desire.

Originally, yt was written as a frontend to HippoDraw, an extensible and
comprehensive plotting package built at SLAC by Paul Kunz.  Over time, however,
it has been generalized to rely on mostly commodity Python packages, and its
dependencies reduced and ease of installation increased.

What yt is and is not
---------------------

In some sense, yt is also designed to be rather utilitarian.  By virtue of the
fact that it has been written in an interpreted language, it can be somewhat
slower than purely C-based analysis codes, although I believe that to be
mitigated by a cleverness of algorithms and a substantially improved
development time for the user in the case of a desire to expand the
functionality.

```EXPAND```

Prerequisites for yt
--------------------

A driving factor in the development of yt over the months leading to release
0.3 has been the reduction of dependencies.  To that extent, only a few
packages are required for the base usage, and a GUI toolkit if you are going to use
the graphical user interface, Reason.

 * `Python <http://python.org/>`_, at least version 2.4, but preferably 2.5.
 * `HDF5 <http://www.hdfgroup.org/>`_, the data storage backend used by Enzo
   and yt (if you can run Enzo, this is already installed!)
 * `NumPy <http://numpy.scipy.org/>`_, the fast numerical backend for Python
 * `MatPlotLib <http://matplotlib.sf.net/>`_, the plotting package
 * `wxPython <http://www.wxpython.org/>`_, the GUI toolkit

(If you are only interested in manipulating data without any graphical plotting
or interfaces, you only need to install NumPy and Python!)

Installing the Necessary Packages
---------------------------------

Installing Python itself is usually quite simple, and often very fast.  Because
we're setting up a small system of packages, even if you have a system-wide
install of python2.5 it can be easier in some cases to create a local directory
structure:

.. code-block:: bash

   $ tar xvf Python-2.5.2.tar.gz
   $ cd Python-2.5.2
   $ ./configure --prefix=$HOME/local/
   $ make install

This will create (if necessary) a directory named local in your home directory,
along with the necessary subdirectories.  When the executable
``$HOME/lcoal/bin/python2.4`` is used to install a package, it will install it
to the ``$HOME/local/`` directory structure.

The python packages are fairly straightforward to install.  The process of
installing packages in python has been greatly simplified over the last few
years with the addition of setuptools, but for these particular packages I
typically recommend installing from source, which for Python packages consists
of changing to the source directory and issues the command:

.. code-block:: bash

   $ tar xvfz $PKGNAME.tar.gz
   $ cd $PKGNAME
   $ python2.5 setup.py install

This method works for NumPy, Matplotlib and yt itself, but for wxPython, I
**strongly** suggest you seek binaries for your platform.  If they are not
available, I recommend you read the ``INSTALL`` file and follow its directions
closely.

