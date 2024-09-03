.. _installing-yt:

Getting and Installing yt
=========================

.. contents::
   :depth: 2
   :local:
   :backlinks: none

.. _getting-yt:

Disclaimer
----------

The Python ecosystem offers many viable tools to setup isolated
Python environments, including but not restricted to

- `venv <https://docs.python.org/3/library/venv.html>`_ (part of the Python standard library)
- `Anaconda/conda <https://docs.conda.io/en/latest/index.html>`_
- `virtualenv <https://virtualenv.pypa.io/en/latest/>`_

We strongly recommend you choose and learn one. However, it is beyond the
scope of this page to cover every situation.

We will show you how to install a stable release or from source, using conda
or pip, and we will *assume* that you do so in an isolated environment.

Also note that each yt release supports a limited range of Python versions.
Here's a summary for most recent releases

+------------+------------+----------------+-----------------+
| yt release | Python 2.7 | Python3 min    | Python3 max     |
+============+============+================+=================+
| 4.4.x      | no         | 3.10.3         | 3.13 (expected) |
+------------+------------+----------------|-----------------|
| 4.3.x      | no         | 3.9.2          | 3.12            |
+------------+------------+----------------+-----------------+
| 4.2.x      | no         | 3.8            | 3.11            |
+------------+------------+----------------+-----------------+
| 4.1.x      | no         | 3.7            | 3.11            |
+------------+------------+----------------+-----------------+
| 4.0.x      | no         | 3.6            | 3.10            |
+------------+------------+----------------+-----------------+
| 3.6.x      | no         | 3.5            | 3.8             |
+------------+------------+----------------+-----------------+
| 3.5.x      | yes        | 3.4            | 3.5             |
+------------+------------+----------------+-----------------+

Minimum Python versions are strict requirements, while maximum
indicates the newest version for which the yt development team
provides pre-compiled binaries via PyPI and conda-forge.
It may be possible to compile existing yt versions under more
recent Python versions, though this is never guaranteed.

yt also adheres to `SPEC 0 <https://scientific-python.org/specs/spec-0000/>`_ as a soft
guideline for our support policy of core dependencies (Python, numpy, matplotlib ...).


Getting yt
----------

In this document we describe several methods for installing yt. The method that
will work best for you depends on your precise situation:

* If you need a stable build, see :ref:`install-stable`

* If you want to build the development version of yt see :ref:`install-from-source`.

.. _install-stable:

Installing a stable release
+++++++++++++++++++++++++++

The latest stable release can be obtained from PyPI with pip

.. code-block:: bash

  $ python -m pip install --upgrade pip
  $ python -m pip install --user yt


Or using the Anaconda/Miniconda Python distributions

.. code-block:: bash

  $ conda install --channel conda-forge yt



.. _install-additional:

Additional requirements (pip only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yt knows about several data file formats. In many cases (e.g. HDF5), additional
dependencies are needed to enable parsing, that are not installed with yt by default.
In order to install all required packages for a specific format alongside yt,
one can specify them as, for instance

.. code-block:: bash

  $ python -m pip install --upgrade pip
  $ python -m pip install --user "yt[ramses]"

Extra requirements can be combined, separated by commas (say ``yt[ramses,enzo_e]``).
Note that all format names are normalized to lower case.


.. _install-from-source:

Building from source
++++++++++++++++++++

To build yt from source, you need ``git``, and a C compiler (such as ``gcc``
or ``clang``).

Then run

.. code-block:: bash

  $ git clone https://github.com/yt-project/yt
  $ cd yt
  $ python -m pip install --upgrade pip
  $ python -m pip install --user -e .


.. _optional-runtime-deps:

Leveraging optional yt runtime dependencies
+++++++++++++++++++++++++++++++++++++++++++

Some relatively heavy runtime dependencies are not included in your build by
default as they may be irrelevant in your workflow. Common examples include
h5py, mpi4py, astropy or scipy. yt implements a on-demand import mechanism that
allows it to run even when they are not installed *until they're needed*, in
which case it will raise an ``ImportError``, pointing to the missing requirement.

If you wish to get everything from the start, you may specify it when building
yt as by appending ``[full]`` to the target name when calling pip, i.e.,

.. code-block:: bash

  $ # stable release
  $ python -m pip install --user yt[full]
  $ # from source
  $ python -m pip install --user -e .[full]


.. _testing-installation:

Testing Your Installation
+++++++++++++++++++++++++

To make sure everything is installed properly, try running yt at
the command line:

.. code-block:: bash

  $ python -c "import yt"

If this runs without raising errors, you have successfully installed yt. Congratulations!

Otherwise, read the error message carefully and follow any instructions it gives
you to resolve the issue. Do not hesitate to :ref:`contact us <asking-for-help>`
so we can help you figure it out.



.. _updating:

Updating yt
+++++++++++

For pip-based installations:

.. code-block:: bash

  $ python -m pip install --upgrade yt


For conda-based installations:

.. code-block:: bash

  $ conda update yt


For git-based installations (yt installed from source), we provide the following
one-liner facility

.. code-block:: bash

  $ yt update

This will pull any changes from GitHub, and recompile yt if necessary.


Uninstalling yt
+++++++++++++++

If you've installed via pip (either from Pypi or from source)

.. code-block:: bash

  $ python -m pip uninstall yt

Or with conda

.. code-block:: bash

  $ conda uninstall yt


TroubleShooting
---------------

If you are unable to locate the yt executable (i.e. executing ``yt version``
at the bash command line fails), then you likely need to add the
``$HOME/.local/bin`` (or the equivalent on your OS) to your PATH. Some Linux
distributions do not include this directory in the default search path.


Additional Resources
--------------------

.. _distro-packages:

yt Distribution Packages
++++++++++++++++++++++++

Some operating systems have yt pre-built packages that can be installed with the
system package manager. Note that the packages in some of these distributions
may be out of date.

.. note::

  Since the third-party packages listed below are not officially supported by
  yt developers, support should not be sought out on the project mailing lists
  or Slack channels.  All support requests related to these packages should be
  directed to their official maintainers.

While we recommended installing yt with either pip or conda, a number of
third-party packages exist for the distributions listed below.

.. image:: https://repology.org/badge/vertical-allrepos/python:yt.svg?header=yt%20packaging%20status
    :target: https://repology.org/project/python:yt/versions


Intel distribution for Python
+++++++++++++++++++++++++++++

A viable alternative to the installation based on Anaconda is the use of the
`Intel Distribution for Python
<https://software.intel.com/en-us/distribution-for-python>`_. For `Parallel
Computation
<http://yt-project.org/docs/dev/analyzing/parallel_computation.html>`_ on Intel
architectures, especially on supercomputers, a large `performance and
scalability improvement <https://arxiv.org/abs/1910.07855>`_ over several common
tasks has been demonstrated.   See `Parallel Computation
<http://yt-project.org/docs/dev/analyzing/parallel_computation.html>`_ for a
discussion on using yt in parallel. Leveraing this specialized distribution for
yt requires that you install some dependencies from the intel conda channel
before installing yt itself, like so

.. code-block:: bash

  $ conda install -c intel numpy scipy mpi4py cython git sympy ipython matplotlib netCDF4
  $ python -m install --user yt
