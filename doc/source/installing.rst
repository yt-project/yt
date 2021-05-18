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
  python packages, or are working on a machine where you are not allowed to, you
  will probably want to use the bash all-in-one installation script.  This
  creates a python environment using the `miniconda python
  distribution <https://docs.conda.io/en/latest/miniconda.html>`_ and the
  `conda <https://conda.io/en/latest/>`_ package manager inside of a single
  folder in your home directory. See :ref:`install-script` for more details.

* If you use the `Anaconda <https://www.anaconda.com/distribution/>`_ python
  distribution and already have ``conda`` installed, see
  :ref:`anaconda-installation` for details on how to install yt using the
  ``conda`` package manager. Note that this is currently the only supported
  installation mechanism on Windows.

* A viable alternative to the installation based on Anaconda is the use of the
  `Intel Distribution for Python <https://software.intel.com/en-us/distribution-for-python>`_.
  For `Parallel Computation <http://yt-project.org/docs/dev/analyzing/parallel_computation.html>`_
  on Intel architectures, especially on supercomputers, a large
  `performance and scalability improvement <https://arxiv.org/abs/1910.07855>`_
  over several common tasks has been demonstrated.
  Detailed installation instructions are provided below as well, see :ref:`conda-intel-python`.
  No change in the way yt is managed by ``conda`` is required.

* Some operating systems have ``yt`` pre-built packages that can be
  installed with the system package manager. Note that the packages in some of
  these distributions may not be the most recent release. See :ref:`distro-packages`
  for a list of available packages.  You can always get the current stable
  version of ``yt`` via ``conda`` as described in :ref:`anaconda-installation`
  or via ``pip`` as described in :ref:`source-installation`.

* If you want to build a development version of yt or are comfortable with
  compilers and know your way around python packaging,
  :ref:`source-installation` will probably be the best choice. If you have set
  up python using a source-based package manager like `Homebrew
  <https://brew.sh>`_ or `MacPorts <https://www.macports.org/>`_ this choice will
  let you install yt using the python installed by the package
  manager. Similarly, this will also work for python environments set up via
  Linux package managers so long as you have the necessary compilers installed
  (e.g. the ``build-essentials`` package on Debian and Ubuntu).

.. note::
  See `Parallel Computation
  <http://yt-project.org/docs/dev/analyzing/parallel_computation.html>`_
  for a discussion on using yt in parallel.

.. _branches-of-yt:

Branches of yt: ``main``, ``stable``, and ``yt-2.x``
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Before you install yt, you must decide which branch (i.e. version) of the code
you prefer to use:

* ``main`` -- The most up-to-date *development* version with the most current
  features but sometimes unstable (the development version of the next release).
* ``stable`` -- The latest stable release of ``yt-3.x``.
* ``yt-2.x`` -- The last stable release of ``yt-2.x``.

If this is your first time using the code, we recommend using ``stable``, unless
you specifically need some piece of brand-new functionality only available in
``main`` or need to run an old script developed for ``yt-2.x``.  There were major
API and functionality changes made in yt for version 3.0.  For a detailed
description of the changes between versions 2.x (e.g. branch ``yt-2.x``) and 3.x
(e.g. branches ``main`` and ``stable``) see :ref:`yt3differences`.  Lastly, don't
feel like you're locked into one branch when you install yt, because you can
easily change the active branch by following the instructions in
:ref:`switching-between-yt-versions`.

.. _install-script:

All-in-One Installation Script
++++++++++++++++++++++++++++++

Because installation of all of the interlocking parts necessary to install yt
itself can be time-consuming, yt provides an all-in-one installation script
which downloads and builds a fully-isolated installation of Python that includes
NumPy, Matplotlib, H5py, git, and yt.

The install script supports UNIX-like systems, including Linux, OS X, and most
supercomputer and cluster environments. It is particularly suited for deployment
in environments where users do not have root access and can only install
software into their home directory.

Since the install is fully-isolated in a single directory, if you get tired of
having yt on your system, you can just delete the directory and yt and all of
its dependencies will be removed from your system (no scattered files remaining
throughout your system).

.. _installing-yt:

Running the Install Script
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can download the installation script with the following command:

.. code-block:: bash

  $ wget https://raw.githubusercontent.com/yt-project/yt/main/doc/install_script.sh

If you do not have ``wget``, the following should also work:

.. code-block:: bash

  $ curl -OL https://raw.githubusercontent.com/yt-project/yt/main/doc/install_script.sh

By default, the bash install script will create a python environment based on
the `miniconda python distribution <https://docs.conda.io/en/latest/miniconda.html>`_,
and will install yt's dependencies using the `conda
<https://conda.io/en/latest/>`_ package manager. To avoid needing a
compilation environment to run the install script, yt itself will also be
installed using ``conda``.

If you would like to customize your yt installation, you can edit the values of
several variables that are defined at the top of the script.

If you would like to build yt from source, you will need to edit the install
script and set ``INST_YT_SOURCE=1`` near the top. This will clone a copy of the
yt git repository and build yt form source. The default is
``INST_YT_SOURCE=0``, which installs yt from a binary conda package.

In addition, you can tell the install script to download and install some
additional packages --- currently these include
`PyX <http://pyx.sourceforge.net/>`_, the `Rockstar halo
finder <https://arxiv.org/abs/1110.4372>`_, `SciPy <https://www.scipy.org/>`_,
`Astropy <https://www.astropy.org/>`_,
`Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_,
and the necessary dependencies for
:ref:`unstructured mesh rendering <unstructured_mesh_rendering>`. The script has
all of the options for installing optional packages near the top of the
file. You should be able to open it and edit it without any knowledge of bash
syntax. For example, to install scipy, change ``INST_SCIPY=0`` to
``INST_SCIPY=1``.

To execute the install script, run:

.. code-block:: bash

  $ bash install_script.sh

Because the installer is downloading and building a variety of packages from
source, this will likely take a few minutes, especially if you have a slow
internet connection. You will get updates of its status at the command prompt
throughout.

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
your shell environment to use yt.

In particular, you will need to ensure that the installation's ``yt-conda/bin``
directory is prepended to your ``PATH`` environment variable.

For Bash-style shells, you can use the following command in a terminal session
to temporarily activate the yt installation:

.. code-block:: bash

  $ export PATH=/path/to/yt-conda/bin:$PATH

and on csh-style shells:

.. code-block:: csh

  $ setenv PATH /path/to/yt-conda/bin:$PATH

If you would like to permanently activate yt, you can also update the init file
appropriate for your shell and OS (e.g. .bashrc, .bash_profile, .cshrc, .zshrc)
to include the same command.

.. _updating-yt:

Updating yt and Its Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With many active developers, code development sometimes occurs at a furious
pace in yt.  To make sure you're using the latest version of the code, run
this command at a command-line:

.. code-block:: bash

  $ conda update yt

If you want to update your dependencies, run:

.. code-block:: bash

   $ conda update --all

If you have installed yt from source, you can use the following command to get
the latest development version of yt:

.. code-block:: bash

   $ yt update

.. _removing-yt:

Removing yt and Its Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because yt and its dependencies are installed in an isolated directory when
you use the script installer, you can easily remove yt and all of its
dependencies cleanly. Simply remove the install directory and its
subdirectories::

  $ rm -rf yt-conda

If you *really* had problems with the installation process, this is a last
defense for solving: remove and then fully :ref:`re-install <installing-yt>`
from the install script again.

.. _anaconda-installation:

Installing yt Using Anaconda
++++++++++++++++++++++++++++

For both the Anaconda and Miniconda installations, make sure that the Anaconda
``bin`` directory is in your path, and then issue:

.. code-block:: bash

  $ conda install -c conda-forge yt

which will install stable branch of yt along with all of its dependencies.

.. _conda-intel-python:

Using the Intel Distribution for Python from conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use conda, you can install yt with the
`Intel Distribution for Python <https://software.intel.com/en-us/distribution-for-python>`_
(recommended for performance in parallel computations on Intel architectures) instead of
the standard Anaconda distribution. First you need to add the intel channel:

.. code-block:: bash

   $ conda config --add channels intel

If you want, at this point you can create a separate environment and switch to it:

.. code-block:: bash

   $ conda create -c intel -n yt_intel
   $ conda activate yt_intel

Now you need to install the remaining yt dependencies in your current environment.
The following provides the Intel-optimized versions of these underlying packages:

.. code-block:: bash

   $ conda config --add channels intel
   $ conda install -c intel numpy scipy mpi4py cython git sympy ipython matplotlib netCDF4

Then you can install yt normally, either from the conda-forge channel as above, or from source (see below).

.. _distro-packages:

yt Distribution Packages
^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

  Since the third-party packages listed below are not officially supported by
  yt developers, support should not be sought out on the project mailing lists
  or Slack channels.  All support requests related to these packages should be
  directed to their official maintainers.

While we recommended installing ``yt`` with either the ``conda`` or ``pip``
package managers, a number of third-party packages exist for the distributions
listed below.  If you can't find your distro here, you can always install
``yt``'s current stable version using ``conda`` or ``pip``, or build the latest
development version from source.

.. image:: https://repology.org/badge/vertical-allrepos/python:yt.svg?header=yt%20packaging%20status
    :target: https://repology.org/project/python:yt/versions

.. note::

  Please be aware that the packages in some of these distributions may be out-of-date!

.. _conda-source-build:

Building yt from Source For Conda-based Installs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, ensure that you have all build dependencies installed in your current
conda environment:

.. code-block:: bash

  $ conda install -c conda-forge cython git sympy ipython matplotlib netCDF4

In addition, you will need a C compiler installed.

Clone the yt repository with:

.. code-block:: bash

  $ git clone https://github.com/yt-project/yt

Once inside the yt directory, update to the appropriate branch and run
``pip install -e .``. For example, the following commands will allow
you to see the tip of the development branch.

.. code-block:: bash

  $ git checkout main
  $ pip install -e .

This will make sure you are running a version of yt corresponding to the
most up-to-date source code.

.. note::

  Alternatively, you can replace ``pip install -e .`` with ``conda develop -b .``.

.. _windows-installation:

Installing yt on Windows
^^^^^^^^^^^^^^^^^^^^^^^^

Installation on 64-bit Microsoft Windows platforms is supported using Anaconda
(see :ref:`anaconda-installation`) and via ``pip``.

.. _source-installation:

Installing yt Using ``pip``
+++++++++++++++++++++++++++

If you already have a python installation that you manage using ``pip`` you can
install the latest release of yt by doing::

  $ pip install yt

If you do not have root access you may need to append ``--user`` to install to a
location in your home folder.

Installing yt from source
+++++++++++++++++++++++++

.. note::

  If you wish to install yt from source in a conda-based installation of yt,
  see :ref:`conda-source-build`.

To install yt from source, you must make sure you have yt's dependencies
installed on your system. Right now, the dependencies to build yt from
source include:

- ``git``
- A C compiler such as ``gcc`` or ``clang``
- ``Python >= 3.6``

In addition, building yt from source requires ``numpy`` and ``cython``
which can be installed with ``pip``:

.. code-block:: bash

  $ pip install --upgrade numpy cython

You may also want to install some of yt's optional dependencies, including
``jupyter``, ``h5py`` (which in turn depends on the HDF5 library), ``scipy``,
``astropy``, or ``cartopy``.

The source code for yt may be found on GitHub. If you prefer to install the
development version of yt instead of the latest stable release, you will need
``git`` to clone the official repo:

.. code-block:: bash

  $ git clone https://github.com/yt-project/yt
  $ cd yt
  $ git checkout main
  $ pip install . --user --install-option="--prefix="

.. note::

  If you maintain your own user-level python installation separate from the
  OS-level python installation, you can leave off ``--user --install-option="--prefix="``, although
  you might need ``sudo`` depending on where python is installed. See `This
  StackOverflow discussion
  <https://stackoverflow.com/questions/4495120/combine-user-with-prefix-error-with-setup-py-install>`_
  if you are curious why ``--install-option="--prefix="`` is necessary on some systems.

This will install yt into a folder in your home directory
(e.g. ``$HOME/.local/lib64/python<py_version>/site-packages`` on Linux,
``$HOME/Library/Python/<py_version>/lib/python/site-packages/`` on OSX) Please refer to
the ``setuptools`` documentation for the additional options.

If you are unable to locate the ``yt`` executable (i.e. executing ``yt version``
at the bash command line fails), then you likely need to add the
``$HOME/.local/bin`` (or the equivalent on your OS) to your PATH. Some Linux
distributions do not include this directory in the default search path.

If you choose this installation method, you do not need to run any activation
script since this will install yt into your global python environment.

If you will be modifying yt, you can also make the clone of the yt git
repository the "active" installed copy:

.. code-block:: bash

  $ git clone https://github.com/yt-project/yt
  $ cd yt
  $ git checkout main
  $ pip install -e . --user --install-option="--prefix="

As above, you can leave off ``--user --install-option="--prefix="`` if you want to install yt into
the default package install path.  If you do not have write access for this
location, you might need to use ``sudo``.

Build errors with ``setuptools`` or ``distribute``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building yt requires version 18.0 or higher of ``setuptools``. If you see error
messages about this package, you may need to update it. For example, with pip
via

.. code-block:: bash

  $ pip install --upgrade setuptools

or your preferred method. If you have ``distribute`` installed, you may also see
error messages for it if it's out of date. You can update with pip via

.. code-block:: bash

  $ pip install --upgrade distribute

or via your preferred method.

Keeping yt Updated via Git
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to maintain your yt installation via updates straight from the
GitHub repository or if you want to do some development on your own, we
suggest you check out some of the :ref:`development docs <contributing-code>`,
especially the sections on :ref:`Git <git-with-yt>` and
:ref:`building yt from source <building-yt>`.

You can also make use of the following command to keep yt up to date from the
command line:

.. code-block:: bash

  $ yt update

This will detect that you have installed yt from the git repository, pull any
changes from GitHub, and then recompile yt if necessary.

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
figure it out.  There is also information at :ref:`update-errors`.

If you like, this might be a good time to run the test suite, see :ref:`testing`
for more details.

.. _switching-between-yt-versions:

Switching versions of yt: ``yt-2.x``, ``stable``, and ``main`` branches
-------------------------------------------------------------------------

Here we explain how to switch between different development branches of yt.

If You Installed yt Using the Bash Install Script
+++++++++++++++++++++++++++++++++++++++++++++++++

The instructions for how to switch between branches depend on whether you ran
the install script with ``INST_YT_SOURCE=0`` (the default) or
``INST_YT_SOURCE=1``. You can determine which option you used by inspecting the
output:

.. code-block:: bash

  $ yt version

If the output from this command looks like:

.. code-block:: none

  The current version and changeset for the code is:

  ---
  Version = 3.2.3
  ---

i.e. it does not refer to a specific changeset hash, then you originally chose
``INST_YT_SOURCE=0``.

On the other hand, if the output from ``yt version`` looks like:

.. code-block:: none

  The current version and changeset for the code is:

  ---
  Version = 3.3-dev
  Changeset = d8eec89b2c86
  ---

i.e. it refers to a specific changeset in the yt git repository, then
you installed using ``INST_YT_SOURCE=1``.

Conda-based installs (``INST_YT_SOURCE=0``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case you can follow the instructions above to build yt from source under conda (see
:ref:`conda-source-build`).

Source-based installs (``INST_YT_SOURCE=1``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You already have the git repository, so you simply need to switch
which version you're using.  Navigate to the root of the yt git
repository, check out the desired version, and rebuild the source (some of the
C code requires a compilation step for big changes like this):

.. code-block:: bash

  $ cd yt-<machine>/src/yt-git
  $ git checkout <desired version>
  $ pip install -e .

Valid versions to jump to are described in :ref:`branches-of-yt`.

You can check which version of yt you have installed by invoking ``yt version``
at the command line.  If you encounter problems, see :ref:`update-errors`.

If You Installed yt Using from Source or Using pip
++++++++++++++++++++++++++++++++++++++++++++++++++

If you have installed python via ``pip``, remove
any extant installations of yt on your system and clone the git
repository of yt as described in :ref:`source-installation`.

.. code-block:: bash

  $ pip uninstall yt
  $ git clone https://github.com/yt-project/yt

Now, to switch between versions, you need to navigate to the root of the git yt
repository. Use git to update to the appropriate version and recompile.

.. code-block:: bash

  $ cd yt
  $ git checkout <desired-version>
  $ pip install . --user --install-option="--prefix="

Valid versions to jump to are described in :ref:`branches-of-yt`).

You can check which version of yt you have installed by invoking ``yt version``
at the command line.  If you encounter problems, see :ref:`update-errors`.
