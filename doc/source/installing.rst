.. _getting-and-installing-yt:

Getting and Installing yt
=========================

.. contents::
   :depth: 2
   :local:
   :backlinks: none

.. _getting-yt:

Getting yt
----------

In this document we describe several methods for installing yt. The method that
will work best for you depends on your precise situation:

* If you do not have root access on your computer, are not comfortable managing
  python packages, or are working on a supercomputer or cluster computer, you
  will probably want to use the bash all-in-one installation script.  This builds 
  python, numpy, matplotlib, and yt from source to set up an isolated scientific 
  python environment inside of a single folder in your home directory. See
  :ref:`install-script` for more details.

* If you use the `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ python
  distribution see :ref:`anaconda-installation` for details on how to install
  yt using the ``conda`` package manager.  Source-based installation from the
  mercurial repository or via ``pip`` should also work under Anaconda. Note that
  this is currently the only supported installation mechanism on Windows.

* If you already have a scientific python software stack installed on your
  computer and are comfortable installing python packages,
  :ref:`source-installation` will probably be the best choice. If you have set
  up python using a source-based package manager like `Homebrew
  <http://brew.sh>`_ or `MacPorts <http://www.macports.org/>`_ this choice will
  let you install yt using the python installed by the package manager. Similarly
  for python environments set up via linux package managers so long as you
  have the the necessary compilers installed (e.g. the ``build-essentials``
  package on debian and ubuntu).

.. _install-script:

All-in-One Installation Script
++++++++++++++++++++++++++++++

Because installation of all of the interlocking parts necessary to install yt
itself can be time-consuming, yt provides an all-in-one installation script
which downloads and builds a fully-isolated Python + NumPy + Matplotlib + HDF5 +
Mercurial installation. Since the install script compiles yt's dependencies from
source, you must have C, C++, and optionally Fortran compilers installed.

The install script supports UNIX-like systems, including Linux, OS X, and most
supercomputer and cluster environments. It is particularly suited for deployment
in environments where users do not have root access and can only install
software into their home directory.

Since the install is fully-isolated in a single directory, if you get tired of
having yt on your system, you can just delete the directory and yt and all of
its dependencies will be removed from your system (no scattered files remaining
throughout your system).

Running the Install Script
^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the installation script, download it from:

.. code-block:: bash

  wget http://hg.yt-project.org/yt/raw/stable/doc/install_script.sh

.. _installing-yt:

By default, the bash install script will install an array of items, but there
are additional packages that can be downloaded and installed (e.g. SciPy, enzo,
etc.). The script has all of these options at the top of the file. You should be
able to open it and edit it without any knowledge of bash syntax.  To execute
it, run:

.. code-block:: bash

  bash install_script.sh

Because the installer is downloading and building a variety of packages from
source, this will likely take a while (e.g. 20 minutes), but you will get 
updates of its status at the command line throughout.

If you receive errors during this process, the installer will provide you 
with a large amount of information to assist in debugging your problems.  The 
file ``yt_install.log`` will contain all of the ``stdout`` and ``stderr`` from 
the entire installation process, so it is usually quite cumbersome.  By looking 
at the last few hundred lines (i.e. ``tail -500 yt_install.log``), you can 
potentially figure out what went wrong.  If you have problems, though, do not 
hesitate to :ref:`contact us <asking-for-help>` for assistance.

.. _activating-yt:

Activating Your Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the installation has completed, there will be instructions on how to set up 
your shell environment to use yt by executing the activate script.  You must 
run this script in order to have yt properly recognized by your system.  You can 
either add it to your login script, or you must execute it in each shell session 
prior to working with yt.

.. code-block:: bash

  source <yt installation directory>/bin/activate

If you use csh or tcsh as your shell, activate that version of the script:

.. code-block:: bash

  source <yt installation directory>/bin/activate.csh

If you don't like executing outside scripts on your computer, you can set 
the shell variables manually.  ``YT_DEST`` needs to point to the root of the
directory containing the install. By default, this will be ``yt-<arch>``, where
``<arch>`` is your machine's architecture (usually ``x86_64`` or ``i386``). You 
will also need to set ``LD_LIBRARY_PATH`` and ``PYTHONPATH`` to contain 
``$YT_DEST/lib`` and ``$YT_DEST/python2.7/site-packages``, respectively.

.. _updating-yt:

Updating yt and Its Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With many active developers, code development sometimes occurs at a furious
pace in yt.  To make sure you're using the latest version of the code, run
this command at a command-line:

.. code-block:: bash

  yt update

Additionally, if you want to make sure you have the latest dependencies
associated with yt and update the codebase simultaneously, type this:

.. code-block:: bash

  yt update --all

.. _removing-yt:

Removing yt and Its Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because yt and its dependencies are installed in an isolated directory when
you use the script installer, you can easily remove yt and all of its
dependencies cleanly.  Simply remove the install directory and its
subdirectories and you're done.  If you *really* had problems with the
code, this is a last defense for solving: remove and then fully
:ref:`re-install <installing-yt>` from the install script again.

.. _anaconda-installation:

Installing yt Using Anaconda
++++++++++++++++++++++++++++

Perhaps the quickest way to get yt up and running is to install it using the
`Anaconda Python Distribution <https://store.continuum.io/cshop/anaconda/>`_,
which will provide you with a easy-to-use environment for installing Python
packages.

If you do not want to install the full anaconda python distribution, you can
install a bare-bones Python installation using miniconda.  To install miniconda,
visit http://repo.continuum.io/miniconda/ and download a recent version of the
``Miniconda-x.y.z`` script (corresponding to Python 2.7) for your platform and
system architecture. Next, run the script, e.g.:

.. code-block:: bash

  bash Miniconda-3.3.0-Linux-x86_64.sh

Make sure that the Anaconda ``bin`` directory is in your path, and then issue:

.. code-block:: bash

  conda install yt

which will install yt along with all of its dependencies.

Recipes to build conda packages for yt are available at
https://github.com/conda/conda-recipes.  To build the yt conda recipe, first
clone the conda-recipes repository

.. code-block:: bash

  git clone https://github.com/conda/conda-recipes

Then navigate to the repository root and invoke `conda build`:

.. code-block:: bash

  cd conda-recipes
  conda build ./yt/

Note that building a yt conda package requires a C compiler.

.. _windows-installation:

Installing yt on Windows
^^^^^^^^^^^^^^^^^^^^^^^^

Installation on 64-bit Microsoft Windows platforms is supported using Anaconda (see
:ref:`anaconda-installation`). Also see :ref:`windows-developing` for details on how to build yt
from source in Windows.

.. _source-installation:

Installing yt Using pip or from Source
++++++++++++++++++++++++++++++++++++++

To install yt from source, you must make sure you have yt's dependencies
installed on your system.  These include: a C compiler, ``HDF5``, ``python``,
``Cython``, ``NumPy``, ``matplotlib``, ``sympy``, and ``h5py``. From here, you
can use ``pip`` (which comes with ``Python``) to install the latest stable
version of yt:

.. code-block:: bash

  $ pip install yt

The source code for yt may be found at the Bitbucket project site and can also
be utilized for installation. If you prefer to install the development version
of yt instead of the latest stable release, you will need ``mercurial`` to clone
the official repo:

.. code-block:: bash

  hg clone https://bitbucket.org/yt_analysis/yt
  cd yt
  hg update yt
  python setup.py install --user --prefix=

.. note::

  If you maintain your own user-level python installation separate from the OS-level python
  installation, you can leave off ``--user --prefix=``, although you might need
  ``sudo`` depending on where python is installed. See `This StackOverflow
  discussion
  <http://stackoverflow.com/questions/4495120/combine-user-with-prefix-error-with-setup-py-install>`_
  if you are curious why ``--prefix=`` is neccessary on systems.

This will install yt into a folder in your home directory
(``$HOME/.local/lib64/python2.7/site-packages`` on Linux,
``$HOME/Library/Python/2.7/lib/python/site-packages/`` on OSX) Please refer to
the ``setuptools`` documentation for the additional options.

If you choose this installation method, you do not need to run any activation
script since this will install yt into your global python environment.

If you will be modifying yt, you can also make the clone of the yt mercurial
repository the "active" installed copy:

.. code-block:: bash

  hg clone https://bitbucket.org/yt_analysis/yt
  cd yt
  hg update yt
  python setup.py develop --user --prefix=

As above, you can leave off ``--user --prefix=`` if you want to install yt into the default
package install path.  If you do not have write access for this location, you
might need to use ``sudo``.

Keeping yt Updated via Mercurial
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to maintain your yt installation via updates straight from the
Bitbucket repository or if you want to do some development on your own, we
suggest you check out some of the :ref:`development docs <contributing-code>`,
especially the sections on :ref:`Mercurial <mercurial-with-yt>` and
:ref:`building yt from source <building-yt>`.

You can also make use of the following command to keep yt up to date from the
command line:

.. code-block:: bash

  yt update

This will detect that you have installed yt from the mercurial repository, pull
any changes from bitbucket, and then recompile yt if necessary.

.. _testing-installation:

Testing Your Installation
-------------------------

To test to make sure everything is installed properly, try running yt at
the command line:

.. code-block:: bash

  yt --help

If this works, you should get a list of the various command-line options for
yt, which means you have successfully installed yt.  Congratulations!

If you get an error, follow the instructions it gives you to debug the problem.
Do not hesitate to :ref:`contact us <asking-for-help>` so we can help you
figure it out.  There is also information at :ref:`update-errors`.

If you like, this might be a good time to run the test suite, see :ref:`testing`
for more details.

.. _switching-between-yt-versions:

Switching between yt-2.x and yt-3.x
-----------------------------------

With the release of version 3.0 of yt, development of the legacy yt 2.x series
has been relegated to bugfixes.  That said, we will continue supporting the 2.x
series for the forseeable future.  This makes it easy to use scripts written
for older versions of yt without substantially updating them to support the
new field naming or unit systems in yt version 3.

Currently, the yt-2.x codebase is contained in a named branch in the yt
mercurial repository.  Thus, depending on the method you used to install
yt, there are different instructions for switching versions.

If You Installed yt Using the Installer Script
++++++++++++++++++++++++++++++++++++++++++++++

You already have the mercurial repository, so you simply need to switch
which version you're using.  Navigate to the root of the yt mercurial
repository, update to the desired version, and rebuild the source (some of the
c code requires a compilation step for big changes like this):

.. code-block:: bash

  cd yt-<machine>/src/yt-hg
  hg update <desired-version>
  python setup.py develop

Valid versions to jump to are:

* ``yt`` -- The latest *dev* changes in yt-3.x (can be unstable)
* ``stable`` -- The latest stable release of yt-3.x
* ``yt-2.x`` -- The latest stable release of yt-2.x
    
You can check which version of yt you have installed by invoking ``yt version``
at the command line.  If you encounter problems, see :ref:`update-errors`.

If You Installed yt Using from Source or Using pip
++++++++++++++++++++++++++++++++++++++++++++++++++

If you have installed python via ``pip``, remove 
any extant installations of yt on your system and clone the source mercurial 
repository of yt as described in :ref:`source-installation`.

.. code-block:: bash

  pip uninstall yt
  hg clone https://bitbucket.org/yt_analysis/yt

Now, to switch between versions, you need to navigate to the root of
the mercurial yt repository. Use mercurial to
update to the appropriate version and recompile.  

.. code-block:: bash

  cd yt
  hg update <desired-version>
  python setup.py install --user --prefix=

Valid versions to jump to are:

* ``yt`` -- The latest *dev* changes in yt-3.x (can be unstable)
* ``stable`` -- The latest stable release of yt-3.x
* ``yt-2.x`` -- The latest stable release of yt-2.x
    
You can check which version of yt you have installed by invoking ``yt version``
at the command line.  If you encounter problems, see :ref:`update-errors`.
