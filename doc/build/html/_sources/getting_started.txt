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

Installing From Source
======================

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
or interfaces, you only need to install HDF5, NumPy, and Python!)

Installing the Necessary Packages
---------------------------------

.. note:: 
   In the ``doc/`` directory in the yt source distribution, there is a script,
   ``install_script.sh``, that I have used in the past to set up a full
   installation of yt.  It may need tweaking or modification, but it gives a
   good idea of the roadmap to installation.

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

