.. _install-script:

Disclaimer
----------

There is no guarantee that the installation script will be maintained or
updated, and it may also be deprecated or removed in the future. Please see the
:ref:`prefered installation methods<installing-yt>`.



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

By default, the bash install script will create a Python environment based on
the `miniconda Python distribution <https://docs.conda.io/en/latest/miniconda.html>`_,
and will install yt's dependencies using the `conda
<https://conda.io/en/latest/>`_ package manager. To avoid needing a
compilation environment to run the install script, yt itself will also be
installed using conda.

If you would like to customize your yt installation, you can edit the values of
several variables that are defined at the top of the script.

If you would like to build yt from source, you will need to edit the install
script and set ``INST_YT_SOURCE=1`` near the top. This will clone a copy of the
yt git repository and build yt form source. The default is
``INST_YT_SOURCE=0``, which installs yt from a binary conda package.

In addition, you can tell the install script to download and install some
additional packages --- including
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

Installations of ``conda`` through miniconda or the anaconda distribution often
offer to make these changes to your path on your behalf.

.. _updating-yt:

Updating yt and Its Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With many active developers, code development sometimes occurs at a rapid
pace in yt.  To automatically update to the latest version of the code, run
this command at a command-line:

.. code-block:: bash

  $ conda update yt

If you want to update your dependencies, run:

.. code-block:: bash

   $ conda update --all


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


.. _switching-between-yt-versions:

Switching branches:
-------------------

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

Stable installs (``INST_YT_SOURCE=0``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case you can follow the :ref:`generic instructions <install-stable>`.

Source-based installs (``INST_YT_SOURCE=1``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You already have the git repository, so you simply need to switch which version
you're using.

In addition to its development branch (``main``), yt has a couple long-lived
branches. ``stable`` follows the stable releases as they are published,
``yt-2.x`` and ``yt-3.x`` are legacy branches for the 2.x and 3.x series
respectively.



Navigate to the root of the yt git repository, check out the
desired version with, for instance

.. code-block::bash

    $ git checkout stable

Finally, recompile: follow :ref:`the instructions to install from source
<install-from-source>`, skipping the git-cloning part.

You can check which version of yt is installed by invoking ``yt version`` at the
command line.
