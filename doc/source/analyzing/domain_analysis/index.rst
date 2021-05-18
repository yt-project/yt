.. _domain-analysis:

Domain-Specific Analysis
========================

yt powers a number modules that provide specialized analysis tools
relevant to one or a few domains. Some of these are internal to yt,
but many exist as external packages, either maintained by the yt
project or independently.

Internal Analysis Modules
-------------------------

These modules exist within yt itself.

.. note::

   As of yt version 3.5, most of the astrophysical analysis tools
   have been moved to the :ref:`yt-astro` and :ref:`attic`
   packages. See below for more information.

.. toctree::
   :maxdepth: 2

   cosmology_calculator
   clump_finding
   xray_emission_fields

External Analysis Modules
-------------------------

These are external packages maintained by the yt project.

.. _yt-astro:

yt Astro Analysis
^^^^^^^^^^^^^^^^^

Source: https://github.com/yt-project/yt_astro_analysis

Documentation: https://yt-astro-analysis.readthedocs.io/

The ``yt_astro_analysis`` package houses most of the astrophysical
analysis tools that were formerly in the ``yt.analysis_modules``
import. These include halo finding, custom halo analysis, synthetic
observations, and exports to radiative transfer codes. See
:ref:`yt_astro_analysis:modules` for a list of available
functionality.

.. _attic:

yt Attic
^^^^^^^^

Source: https://github.com/yt-project/yt_attic

Documentation: https://yt-attic.readthedocs.io/

The ``yt_attic`` contains former yt analysis modules that have
fallen by the wayside. These may have small bugs or were simply
not kept up to date as yt evolved. Tools in here are looking for
a new owner and a new home. If you find something in here that
you'd like to bring back to life, either by adding it to
:ref:`yt-astro` or as part of your own package, you are welcome
to it! If you'd like any help, let us know! See
:ref:`yt_attic:attic_modules` for a list of inventory of the
attic.

Extensions
----------

There are a number of independent, yt-related packages for things
like visual effects, interactive widgets, synthetic absorption
spectra, X-ray observations, and merger-trees. See the
`yt Extensions <http://yt-project.org/extensions.html>`_ page for
a list of available extension packages.
