===============
Getting Started
===============

Maintained Installations
========================

As of right now, I maintain an installation of yt on several machines.  I try
to keep them up to date with the stable branch.

DataStar (SDSC)
---------------

To use yt on 
`DataStar <http://www.sdsc.edu/us/resources/datastar/>`_,
you need to ensure you are using the correct
installation of Python 2.4 and the correct set of python packages.  In my
.bashrc I have:

.. code-block:: bash

   export PATH=/usr/local/apps/python_2.4.2_64/bin:$PATH
   export PYTHONPATH=/users/stanford/mturk/local/lib/python2.4/site-packages/

which ensures that you are using my installation of yt and the correct, 64-bit
global installation of Python2.4.

Orange and Red (SLAC)
---------------------

To use yt on `Orange <http://kipac.stanford.edu/collab/computing/hardware/orange>`_
or `Red <http://www.sgi.com/company_info/newsroom/press_releases/2005/april/space_sciences.html>`_
you have to set your pythonpath appropriately.  I have created a script called
'pywrapper.sh' that sets up your PATH, PYTHONPATH and LD_LIBRARY_PATH.

.. code-block:: bash

   $ export ARCH_PATH="/u/ki/mturk/ki12/`uname -p`_local/"
   $ $ARCH_PATH/bin/pywrapper.sh my_script.py

where my_script.py is the script you wish to run.

Binary Packages
===============

Currently binary packages are only supplied for OSX.  See the download page on
the Wiki for up-to-date links.

Installing From Source
======================

Using the Installation Script
-----------------------------

.. note:: The installation script is now the preferred means of installing a
   full set of packages -- but if you are comfortable with python, feel free to 

In the ``doc/`` directory in the yt source distribution, there is a script,
``install_script.sh``, designed to set up a full installation of yt, along with
all the necessary dependencies.  You can run this script from within a checkout
of yt or an expanded tarball.

.. note:: For convenience, yt will be installed in 'develop' mode, which means
   any changes in the source directory will be included the next time you
   import yt!

There are several variabels you can set inside this script.

   ``DEST_DIR``
     This is the location to which all source code will be downloaded and
     resulting built libraries installed.
   ``HDF5_DIR``
     If you wish to link against existing HDF5 (*shared*) libraries, put the
     root path to the installation here.
   ``INST_WXPYTHON``
     This is a boolean, set to 0 or 1, that governs whether or not wxPython
     should be installed.
   ``INST_ZLIB``
     This is a boolean, set to 0 or 1, that governs whether or not zlib
     should be installed.
   ``YT_DIR``
     If you've got a source checkout of YT somewhere else, point to it with
     this!

.. warning:: If you run into problems, particularly anything involving
   ``-fPIC``, it is likely that there's a problem with static libraries.
   Try asking the installer script to install HDF5 and ZLIB.

Installing by Hand
------------------

If you've ever installed a python package by hand before, YT should be easy to
install.  You will need to install the prequisites first.  A driving factor in
the development of yt over the months leading to release 1.0 has been the
reduction of dependencies.  To that extent, only a few packages are required
for the base usage, and a GUI toolkit if you are going to use the graphical
user interface, Reason.

 * `Python <http://python.org/>`_, at least version 2.4, but preferably 2.5 or
   2.6.
 * `HDF5 <http://www.hdfgroup.org/>`_, the data storage backend used by Enzo
   and yt (if you can run Enzo, this is already installed!)
 * `NumPy <http://numpy.scipy.org/>`_, the fast numerical backend for Python
 * `MatPlotLib <http://matplotlib.sf.net/>`_, the plotting package
 * `wxPython <http://www.wxpython.org/>`_, the GUI toolkit (optional)

(If you are only interested in manipulating data without any graphical plotting
or interfaces, you only need to install HDF5, NumPy, and Python!)

Instructions for installing these packages is, unfortunately, beyond the scope
of this document.  However, there are copious directions on how to do so
elsewhere.  You may also consider installing the
`Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_,
which includes all of the necessary packages.

Once these dependencies have been met, YT can be installed in the standard
manner:

.. code-block:: bash

   cd yt/
   python2.5 setup.py install --prefix=/some/where/
