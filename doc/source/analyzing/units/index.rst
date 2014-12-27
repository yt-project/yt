.. _units:

Symbolic Units
==============

This section describes yt's symbolic unit capabilities. This is provided as
quick introduction for those who are already familiar with yt but want to learn
more about the unit system.  Please see :ref:`analyzing` and :ref:`visualizing`
for more detail about querying, analyzing, and visualizing data in yt.

Each subsection is a notebook.  To open these notebooks in a "live" IPython session
and execute the documentation interactively, you need to download the repository
and start the IPython notebook.

You will then need to navigate to :code:`$YT_HG/doc/source/units` (where $YT_HG
is the location of a clone of the yt mercurial repository), and then start an
IPython notebook server:

.. code:: bash
  
   $ ipython notebook

.. warning:: The pre-filled out notebooks are *far* less fun than running them
             yourself!

Here are the notebooks, which have been filled in for inspection:

.. toctree::
   :maxdepth: 1

   symbolic_units
   fields_and_unit_conversion
   comoving_units_and_code_units
   comparing_units_from_different_datasets
   units_and_plotting
   unit_equivalencies

.. note::

   The notebooks use sample datasets that are available for download at
   http://yt-project.org/data.  See :ref:`quickstart-introduction` for more
   details.

Let us know if you would like to contribute other example notebooks, or have
any suggestions for how these can be improved.
