.. _faq:


Frequently Asked Questions
==========================

.. contents::
   :depth: 2
   :local:
   :backlinks: none

Version & Installation
----------------------

.. _determining-version:

How can I tell what version of yt I'm using?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you run into problems with yt and you're writing to the mailing list
or contacting developers on Slack, they will likely want to know what version of
yt you're using.  Often times, you'll want to know both the yt version,
as well as the last changeset that was committed to the branch you're using.
To reveal this, go to a command line and type:

.. code-block:: bash

    $ yt version

The result will look something like this:

.. code-block:: bash

    yt module located at:
        /Users/mitchell/src/yt-conda/src/yt-git

    The current version of yt is:

    ---
    Version = 4.0.dev0
    Changeset = 9f947a930ab4
    ---
    This installation CAN be automatically updated.


For more information on this topic, see :ref:`updating`.

.. _yt-3.0-problems:

I upgraded to yt 4.0 but my code no longer works.  What do I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We've tried to keep the number of backward-incompatible changes to a minimum
with the release of yt-4.0, but because of the wide-reaching changes to how
yt manages data, there may be updates you have to make.
You can see many of the changes in :ref:`yt4differences`, and
in :ref:`transitioning-to-4.0` there are helpful tips on how to modify your scripts to update them.

Code Errors and Failures
------------------------

Python fails saying that it cannot import yt modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is commonly exhibited with an error about not being able to import code
that is part of yt. This is likely because the code that is failing to import
needs to be compiled or recompiled.

This error tends to occur when there are changes in the underlying Cython files
that need to be rebuilt, like after a major code update or when switching
between distant branches.

This is solved by running the install command again. See
:ref:`install-from-source`.


.. _faq-mpi4py:

yt complains that it needs the mpi4py module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For yt to be able to incorporate parallelism on any of its analysis (see
:ref:`parallel-computation`), it needs to be able to use MPI libraries.
This requires the ``mpi4py`` module to be installed in your version of python.
Unfortunately, installation of ``mpi4py`` is *just* tricky enough to elude the
yt batch installer.  So if you get an error in yt complaining about mpi4py
like:

.. code-block:: bash

    ImportError: No module named mpi4py

then you should install ``mpi4py``.  The easiest way to install it is through
the pip interface.  At the command line, type:

.. code-block:: bash

    $ python -m pip install mpi4py

What this does is it finds your default installation of Python (presumably
in the yt source directory), and it installs the mpi4py module.  If this
action is successful, you should never have to worry about your aforementioned
problems again.  If, on the other hand, this installation fails (as it does on
such machines as NICS Kraken, NASA Pleaides and more), then you will have to
take matters into your own hands.  Usually when it fails, it is due to pip
being unable to find your MPI C/C++ compilers (look at the error message).
If this is the case, you can specify them explicitly as per:

.. code-block:: bash

    $ env MPICC=/path/to/MPICC python -m pip install mpi4py

So for example, on Kraken, I switch to the gnu C compilers (because yt
doesn't work with the portland group C compilers), then I discover that
cc is the mpi-enabled C compiler (and it is in my path), so I run:

.. code-block:: bash

    $ module swap PrgEnv-pgi PrgEnv-gnu
    $ env MPICC=cc python -m pip install mpi4py

And voila!  It installs!  If this *still* fails for you, then you can
build and install from source and specify the mpi-enabled c and c++
compilers in the mpi.cfg file.  See the
`mpi4py installation page <https://mpi4py.readthedocs.io/en/stable/install.html>`_
for details.


Units
-----

.. _conversion-factors:

How do I convert between code units and physical units for my dataset?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting with yt-3.0, and continuing to yt-4.0, yt uses an internal symbolic
unit system.  In yt-3.0 this was bundled with the main yt codebase, and with
yt-4.0 it is now available as a separate package called `unyt
<https://unyt.readthedocs.org/>`_.  Conversion factors are tied up in the
``length_unit``, ``times_unit``, ``mass_unit``, and ``velocity_unit``
attributes, which can be converted to any arbitrary desired physical unit:

.. code-block:: python

    print("Length unit: ", ds.length_unit)
    print("Time unit: ", ds.time_unit)
    print("Mass unit: ", ds.mass_unit)
    print("Velocity unit: ", ds.velocity_unit)

    print("Length unit: ", ds.length_unit.in_units("code_length"))
    print("Time unit: ", ds.time_unit.in_units("code_time"))
    print("Mass unit: ", ds.mass_unit.in_units("kg"))
    print("Velocity unit: ", ds.velocity_unit.in_units("Mpc/year"))

So to accomplish the example task of converting a scalar variable ``x`` in
code units to kpc in yt-4.0, you can do one of two things.  If ``x`` is
already a YTQuantity with units in ``code_length``, you can run:

.. code-block:: python

    x.in_units("kpc")

However, if ``x`` is just a numpy array or native python variable without
units, you can convert it to a YTQuantity with units of ``kpc`` by running:

.. code-block:: python

    x = x * ds.length_unit.in_units("kpc")

For more information about unit conversion, see :ref:`units`.

How do I make a YTQuantity tied to a specific dataset's units?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to create a variable or array that is tied to a particular dataset
(and its specific conversion factor to code units), use the ``ds.quan`` (for
individual variables) and ``ds.arr`` (for arrays):

.. code-block:: python

    import yt

    ds = yt.load(filename)
    one_Mpc = ds.quan(1, "Mpc")
    x_vector = ds.arr([1, 0, 0], "code_length")

You can then naturally exploit the units system:

.. code-block:: python

    print("One Mpc in code_units:", one_Mpc.in_units("code_length"))
    print("One Mpc in AU:", one_Mpc.in_units("AU"))
    print("One Mpc in comoving kpc:", one_Mpc.in_units("kpccm"))

For more information about unit conversion, see :ref:`units`.

.. _accessing-unitless-data:

How do I access the unitless data in a YTQuantity or YTArray?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While there are numerous benefits to having units tied to individual
quantities in yt, they can also produce issues when simply trying to combine
YTQuantities with numpy arrays or native python floats that lack units.  A
simple example of this is::

    # Create a YTQuantity that is 1 kpc in length and tied to the units of
    # dataset ds
    >>> x = ds.quan(1, 'kpc')

    # Try to add this to some non-dimensional quantity
    >>> print(x + 1)

    YTUnitOperationError: The addition operator for YTArrays with units (kpc) and (1) is not well defined.

The solution to this means using the YTQuantity and YTArray objects for all
of one's computations, but this isn't always feasible.  A quick fix for this
is to just grab the unitless data out of a YTQuantity or YTArray object with
the ``value`` and ``v`` attributes, which return a copy, or with the ``d``
attribute, which returns the data itself:

.. code-block:: python

    x = ds.quan(1, "kpc")
    x_val = x.v
    print(x_val)

    array(1.0)

    # Try to add this to some non-dimensional quantity
    print(x + 1)

    2.0

For more information about this functionality with units, see :ref:`units`.

Fields
------

.. _faq-handling-log-vs-linear-space:

How do I modify whether or not yt takes the log of a particular field?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt sets up defaults for many fields for whether or not a field is presented
in log or linear space. To override this behavior, you can modify the
``field_info`` dictionary.  For example, if you prefer that ``density`` not be
logged, you could type:

.. code-block:: python

    ds = load("my_data")
    ds.index
    ds.field_info["gas", "density"].take_log = False

From that point forward, data products such as slices, projections, etc., would
be presented in linear space. Note that you have to instantiate ds.index before
you can access ds.field info.  For more information see the documentation on
:ref:`fields` and :ref:`creating-derived-fields`.

.. _faq-new-field:

I added a new field to my simulation data, can yt see it?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes! yt identifies all the fields in the simulation's output file
and will add them to its ``field_list`` even if they aren't listed in
:ref:`field-list`. These can then be accessed in the usual manner. For
example, if you have created a field for the potential called
``PotentialField``, you could type:

.. code-block:: python

   ds = load("my_data")
   ad = ds.all_data()
   potential_field = ad["PotentialField"]

The same applies to fields you might derive inside your yt script
via :ref:`creating-derived-fields`. To check what fields are
available, look at the properties ``field_list`` and ``derived_field_list``:

.. code-block:: python

   print(ds.field_list)
   print(ds.derived_field_list)

or for a more legible version, try:

.. code-block:: python

   for field in ds.derived_field_list:
       print(field)

.. _faq-add-field-diffs:

What is the difference between ``yt.add_field()`` and ``ds.add_field()``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The global ``yt.add_field()``
(:meth:`~yt.fields.field_info_container.FieldInfoContainer.add_field`)
function is for adding a field for every subsequent dataset that is loaded
in a particular python session, whereas ``ds.add_field()``
(:meth:`~yt.data_objects.static_output.Dataset.add_field`) will only add it
to dataset ``ds``.

Data Objects
------------

.. _ray-data-ordering:

Why are the values in my Ray object out of order?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the Ray objects
(:class:`~yt.data_objects.selection_data_containers.YTOrthoRay` and
:class:`~yt.data_objects.selection_data_containers.YTRay`) with AMR data
gives non-contiguous cell information in the Ray's data array. The
higher-resolution cells are appended to the end of the array.  Unfortunately,
due to how data is loaded by chunks for data containers, there is really no
easy way to fix this internally.  However, there is an easy workaround.

One can sort the ``Ray`` array data by the ``t`` field, which is the value of
the parametric variable that goes from 0 at the start of the ray to 1 at the
end. That way the data will always be ordered correctly. As an example you can:

.. code-block:: python

    my_ray = ds.ray(...)
    ray_sort = np.argsort(my_ray["t"])
    density = my_ray["gas", "density"][ray_sort]

There is also a full example in the :ref:`manual-line-plots` section of the
docs.

Developing
----------

.. _making-a-PR:

Someone asked me to make a Pull Request (PR) to yt.  How do I do that?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A pull request is the action by which you contribute code to yt.  You make
modifications in your local copy of the source code, then *request* that
other yt developers review and accept your changes to the main code base.
For a full description of the steps necessary to successfully contribute
code and issue a pull request (or manage multiple versions of the source code)
please see :ref:`sharing-changes`.

.. _making-an-issue:

Someone asked me to file an issue or a bug report for a bug I found.  How?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`reporting-a-bug` and :ref:`sharing-changes`.

Miscellaneous
-------------

.. _getting-sample-data:

How can I get some sample data for yt?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many different sample datasets can be found at https://yt-project.org/data/ .
These can be downloaded, unarchived, and they will each create their own
directory.  It is generally straight forward to load these datasets, but if
you have any questions about loading data from a code with which you are
unfamiliar, please visit :ref:`loading-data`.

To make things easier to load these sample datasets, you can add the parent
directory to your downloaded sample data to your *yt path*.
If you set the option ``test_data_dir``, in the section ``[yt]``,
in ``~/.config/yt/yt.toml``, yt will search this path for them.

This means you can download these datasets to ``/big_drive/data_for_yt`` , add
the appropriate item to ``~/.config/yt/yt.toml``, and no matter which directory you are
in when running yt, it will also check in *that* directory.

In many cases, these are also available using the ``load_sample`` command,
described in :ref:`loading-sample-data`.


.. _faq-scroll-up:

I can't scroll-up to previous commands inside python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the up-arrow key does not recall the most recent commands, there is
probably an issue with the readline library. To ensure the yt python
environment can use readline, run the following command:

.. code-block:: bash

   $ python -m pip install gnureadline

.. _faq-old-data:

.. _faq-log-level:

How can I change yt's log level?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt's default log level is ``INFO``. However, you may want less voluminous logging,
especially if you are in an IPython notebook or running a long or parallel script.
On the other hand, you may want it to output a lot more, since you can't figure out
exactly what's going wrong, and you want to output some debugging information.
The default yt log level can be changed using the :ref:`configuration-file`,
either by setting it in the ``$HOME/.config/yt/yt.toml`` file:

.. code-block:: bash

   $ yt config set yt log_level 10  # This sets the log level to "DEBUG"

which would produce debug (as well as info, warning, and error) messages, or at runtime:

.. code-block:: python

   yt.set_log_level("error")

This is the same as doing:

.. code-block:: python

   yt.set_log_level(40)

which in this case would suppress everything below error messages. For reference,
the numerical values corresponding to different log levels are:

.. csv-table::
   :header: Level, Numeric Value
   :widths: 10, 10

   ``CRITICAL``,50
   ``ERROR``,40
   ``WARNING``,30
   ``INFO``,20
   ``DEBUG``,10
   ``NOTSET``,0

Can I always load custom data objects, fields, quantities, and colormaps with every dataset?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`plugin-file` provides a means for always running custom code whenever
yt is loaded up.  This custom code can be new data objects, or fields, or
colormaps, which will then be accessible in any future session without having
modified the source code directly.  See the description in :ref:`plugin-file`
for more details.

How do I cite yt?
^^^^^^^^^^^^^^^^^

If you use yt in a publication, we'd very much appreciate a citation!  You
should feel free to cite the `ApJS paper
<https://ui.adsabs.harvard.edu/abs/2011ApJS..192....9T>`_ with the following BibTeX
entry: ::

   @ARTICLE{2011ApJS..192....9T,
      author = {{Turk}, M.~J. and {Smith}, B.~D. and {Oishi}, J.~S. and {Skory}, S. and
   	{Skillman}, S.~W. and {Abel}, T. and {Norman}, M.~L.},
       title = "{yt: A Multi-code Analysis Toolkit for Astrophysical Simulation Data}",
     journal = {The Astrophysical Journal Supplement Series},
   archivePrefix = "arXiv",
      eprint = {1011.3514},
    primaryClass = "astro-ph.IM",
    keywords = {cosmology: theory, methods: data analysis, methods: numerical },
        year = 2011,
       month = jan,
      volume = 192,
         eid = {9},
       pages = {9},
         doi = {10.1088/0067-0049/192/1/9},
      adsurl = {https://ui.adsabs.harvard.edu/abs/2011ApJS..192....9T},
     adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }
