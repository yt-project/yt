.. _getting-and-installing-yt:

Getting and Installing yt
=========================

.. _getting-yt:

Getting yt
----------

yt is a Python package (with some components written in C), using NumPy as a
computation engine, Matplotlib for some visualization tasks and Mercurial for
version control.  Because installation of all of these interlocking parts can 
be time-consuming, yt provides an installation script which downloads and builds
a fully-isolated Python + NumPy + Matplotlib + HDF5 + Mercurial installation.  
yt supports Linux and OSX deployment, with the possibility of deployment on 
other Unix-like systems (XSEDE resources, clusters, etc.).  Windows is not 
supported.

Since the install is fully-isolated, if you get tired of having yt on your 
system, you can just delete its directory, and yt and all of its dependencies
will be removed from your system (no scattered files remaining throughout 
your system).  

To get the installation script, download it from:

.. code-block:: bash

  http://hg.yt-project.org/yt/raw/stable/doc/install_script.sh

.. _installing-yt:

Installing yt
-------------

By default, the bash script will install an array of items, but there are 
additional packages that can be downloaded and installed (e.g. SciPy, enzo, 
etc.). The script has all of these options at the top of the file. You should 
be able to open it and edit it without any knowledge of bash syntax.  
To execute it, run:

.. code-block:: bash

  $ bash install_script.sh

Because the installer is downloading and building a variety of packages from
source, this will likely take a while (e.g. 20 minutes), but you will get 
updates of its status at the command line throughout.

If you receive errors during this process, the installer will provide you 
with a large amount of information to assist in debugging your problems.  The 
file ``yt_install.log`` will contain all of the ``STDOUT`` and ``STDERR`` from 
the entire installation process, so it is usually quite cumbersome.  By looking 
at the last few hundred lines (i.e. ``tail -500 yt_install.log``), you can 
potentially figure out what went wrong.  If you have problems, though, do not 
hesitate to :ref:`contact us <asking-for-help>` for assistance.

.. _activating-yt:

Activating Your Installation
----------------------------

Once the installation has completed, there will be instructions on how to set up 
your shell environment to use yt by executing the activate script.  You must 
run this script in order to have yt properly recognized by your system.  You can 
either add it to your login script, or you must execute it in each shell session 
prior to working with yt.

.. code-block:: bash

  $ source <yt installation directory>/bin/activate

If you use csh or tcsh as your shell, activate that version of the script:

.. code-block:: bash

  $ source <yt installation directory>/bin/activate.csh

If you don't like executing outside scripts on your computer, you can set 
the shell variables manually.  ``YT_DEST`` needs to point to the root of the
directory containing the install. By default, this will be ``yt-<arch>``, where
``<arch>`` is your machine's architecture (usually ``x86_64`` or ``i386``). You 
will also need to set ``LD_LIBRARY_PATH`` and ``PYTHONPATH`` to contain 
``$YT_DEST/lib`` and ``$YT_DEST/python2.7/site-packages``, respectively.

Alternative Installation Methods
--------------------------------

If you want to forego the use of the install script, you need to make sure you
have yt's dependencies installed on your system.  These include: a C compiler,
``HDF5``, ``Freetype``, ``libpng``, ``python``, ``cython``, ``NumPy``, and
``matplotlib``.  From here, you can use ``pip`` (which comes with ``Python``) to
install yt as:

.. code-block:: bash

  $ pip install yt

The source code for yt may be found at the Bitbucket project site and can also be
utilized for installation. If you prefer to use it instead of relying on external
tools, you will need ``mercurial`` to clone the official repo:

.. code-block:: bash

  $ hg clone https://bitbucket.org/yt_analysis/yt
  $ cd yt
  $ hg update yt
  $ python setup.py install --user

It will install yt into ``$HOME/.local/lib64/python2.7/site-packages``. 
Please refer to ``setuptools`` documentation for the additional options.

Provided that the required dependencies are in a predictable location, yt should
be able to find them automatically. However, you can manually specify prefix used
for installation of ``HDF5``, ``Freetype`` and ``libpng`` by using ``hdf5.cfg``,
``freetype.cfg``, ``png.cfg`` or setting ``HDF5_DIR``, ``FTYPE_DIR``, ``PNG_DIR``
environmental variables respectively, e.g.

.. code-block:: bash

  $ echo '/usr/local' > hdf5.cfg
  $ export FTYPE_DIR=/opt/freetype

If you choose this installation method, you do not need to run the activation
script as it is unnecessary.

.. _testing-installation:

Testing Your Installation
-------------------------

To test to make sure everything is installed properly, try running yt at 
the command line:

.. code-block:: bash

  $ yt --help

If this works, you should get a list of the various command-line options for
yt, which means you have successfully installed yt.  Congratulations!  

If you get an error, follow the instructions it gives you to debug the problem.  
Do not hesitate to :ref:`contact us <asking-for-help>` so we can help you 
figure it out.

.. _updating-yt:

Updating yt and its dependencies
--------------------------------

With many active developers, code development sometimes occurs at a furious 
pace in yt.  To make sure you're using the latest version of the code, run
this command at a command-line:

.. code-block:: bash

  $ yt update

Additionally, if you want to make sure you have the latest dependencies 
associated with yt and update the codebase simultaneously, type this:

.. code-block:: bash

  $ yt update --all

.. _removing-yt:

Removing yt and its dependencies
--------------------------------

Because yt and its dependencies are installed in an isolated directory when
you use the script installer, you can easily remove yt and all of its 
dependencies cleanly.  Simply remove the install directory and its 
subdirectories and you're done.  If you *really* had problems with the
code, this is a last defense for solving: remove and then fully
:ref:`re-install <installing-yt>` from the install script again.
