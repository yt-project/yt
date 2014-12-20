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
or contacting developers on IRC, they will likely want to know what version of
yt you're using.  Oftentimes, you'll want to know both the yt version, 
as well as the last changeset that was committed to the branch you're using.  
To reveal this, go to a command line and type:

.. code-block:: bash
    
    $ yt version

    yt module located at:
        /Users/username/src/yt-x86_64/src/yt-hg
    The supplemental repositories are located at:
        /Users/username/src/yt-x86_64/src/yt-supplemental

    The current version and changeset for the code is:

    ---
    Version = 2.7-dev
    Changeset = 6bffc737a67a
    ---

    This installation CAN be automatically updated.
    yt dependencies were last updated on
    Wed Dec  4 15:47:40 MST 2013

    To update all dependencies, run "yt update --all".

If the changeset is displayed followed by a "+", it means you have made 
modifications to the code since the last changeset.

For more information on this topic, see :ref:`updating-yt`.

.. _yt-3.0-problems:

I upgraded to yt 3.0 but my code no longer works.  What do I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because there are a lot of backwards-incompatible changes in yt 3.0 (see 
:ref:`yt3differences`, it can
be a daunting effort in transitioning old scripts from yt 2.x to 3.0.
We have tried to describe the basic process of making that transition
in :ref:`transitioning-to-3.0`.  If you just want to change back to yt 2.x
for a while until you're ready to make the transition, you can follow
the instructions in :ref:`switching-between-yt-versions`.

Code Errors and Failures
------------------------

yt fails saying that it cannot import yt modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is commonly exhibited with this error: 
``ImportError: cannot import name obtain_rvec``.  This is likely because 
you need to rebuild the source.  You can do this automatically by running:

.. code-block:: bash

    cd $YT_HG
    python setup.py develop

where ``$YT_HG`` is the path to the yt mercurial repository.  

This error tends to occur when there are changes in the underlying cython
files that need to be rebuilt, like after a major code update or in switching
from 2.x to 3.x.  For more information on this, see 
:ref:`switching-between-yt-versions`.

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

    pip install mpi4py

What this does is it finds your default installation of python (presumably
in the yt source directory), and it installs the mpi4py module.  If this
action is successful, you should never have to worry about your aforementioned
problems again.  If, on the other hand, this installation fails (as it does on
such machines as NICS Kraken, NASA Pleaides and more), then you will have to
take matters into your own hands.  Usually when it fails, it is due to pip
being unable to find your MPI C/C++ compilers (look at the error message).
If this is the case, you can specify them explicitly as per:

.. code-block:: bash

    env MPICC=/path/to/MPICC pip install mpi4py

So for example, on Kraken, I switch to the gnu C compilers (because yt 
doesn't work with the portland group C compilers), then I discover that
cc is the mpi-enabled C compiler (and it is in my path), so I run:

.. code-block:: bash

    module swap PrgEnv-pgi PrgEnv-gnu
    env MPICC=cc pip install mpi4py

And voila!  It installs!  If this *still* fails for you, then you can 
build and install from source and specify the mpi-enabled c and c++ 
compilers in the mpi.cfg file.  See the 
`mpi4py installation page <http://mpi4py.scipy.org/docs/usrman/install.html>`_ 
for details.


Units
-----

.. _conversion-factors:

How do I get the convert between code units and physical units for my dataset?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converting between physical units and code units is a common task.  In yt-2.x,
the syntax for getting conversion factors was in the units dictionary 
(``pf.units['kpc']``). So in order to convert a variable ``x`` in code units to
kpc, you might run:

.. code-block:: python

    x = x*pf.units['kpc']

In yt-3.0, this no longer works.  Conversion factors are tied up in the 
``length_unit``, ``times_unit``, ``mass_unit``, and ``velocity_unit`` 
attributes, which can be converted to any arbitrary desired physical unit:

.. code-block:: python

    print "Length unit: ", ds.length_unit
    print "Time unit: ", ds.time_unit
    print "Mass unit: ", ds.mass_unit
    print "Velocity unit: ", ds.velocity_unit

    print "Length unit: ", ds.length_unit.in_units('code_length')
    print "Time unit: ", ds.time_unit.in_units('code_time')
    print "Mass unit: ", ds.mass_unit.in_units('kg')
    print "Velocity unit: ", ds.velocity_unit.in_units('Mpc/year')

So to accomplish the example task of converting a scalar variable ``x`` in 
code units to kpc in yt-3.0, you can do one of two things.  If ``x`` is 
already a YTQuantity with units in ``code_length``, you can run:

.. code-block:: python

    x.in_units('kpc')

However, if ``x`` is just a numpy array or native python variable without
units, you can convert it to a YTQuantity with units of ``kpc`` by running:

.. code-block:: python

    x = x*ds.length_unit.in_units('kpc')

For more information about unit conversion, see :ref:`data_selection_and_fields`.

How do I make a YTQuantity tied to a specific dataset's units?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to create a variable or array that is tied to a particular dataset
(and its specific conversion factor to code units), use the ``ds.quan`` (for 
individual variables) and ``ds.arr`` (for arrays):

.. code-block:: python

    import yt
    ds = yt.load(filename)
    one_Mpc = ds.quan(1, 'Mpc')
    x_vector = ds.arr([1,0,0], 'code_length')

You can then naturally exploit the units system:

.. code-block:: python

    print "One Mpc in code_units:", one_Mpc.in_units('code_length')
    print "One Mpc in AU:", one_Mpc.in_units('AU')
    print "One Mpc in comoving kpc:", one_Mpc.in_units('kpccm')

For more information about unit conversion, see :ref:`data_selection_and_fields`.

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
    >>> print x + 1
    
    YTUnitOperationError: The addition operator for YTArrays with units (kpc) and (1) is not well defined.

The solution to this means using the YTQuantity and YTArray objects for all 
of one's computations, but this isn't always feasible.  A quick fix for this 
is to just grab the unitless data out of a YTQuantity or YTArray object with
the ``value`` and ``v`` attributes, which return a copy, or with the ``d`` 
attribute, which returns the data itself:

.. code-block:: python

    x = ds.quan(1, 'kpc')
    x_val = x.v
    print x_val 

    array(1.0)

    # Try to add this to some non-dimensional quantity
    print x + 1

    2.0 

For more information about this functionality with units, see :ref:`data_selection_and_fields`.

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
    ds.field_info['density'].take_log = False

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

   print ds.field_list
   print ds.derived_field_list

or for a more legible version, try:

.. code-block:: python

   for field in ds.derived_field_list: 
       print field

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
(:class:`~yt.data_objects.selection_data_containers.YTOrthoRayBase` and 
:class:`~yt.data_objects.selection_data_containers.YTRayBase`) with AMR data 
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
    density = my_ray["density"][ray_sort]

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

Many different sample datasets can be found at http://yt-project.org/data/ .
These can be downloaded, unarchived, and they will each create their own
directory.  It is generally straight forward to load these datasets, but if
you have any questions about loading data from a code with which you are 
unfamiliar, please visit :ref:`loading-data`.

To make things easier to load these sample datasets, you can add the parent
directory to your downloaded sample data to your *yt path*.
If you set the option ``test_data_dir``, in the section ``[yt]``,
in ``~/.yt/config``, yt will search this path for them.

This means you can download these datasets to ``/big_drive/data_for_yt`` , add
the appropriate item to ``~/.yt/config``, and no matter which directory you are
in when running yt, it will also check in *that* directory.


.. _faq-scroll-up:

I can't scroll-up to previous commands inside python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the up-arrow key does not recall the most recent commands, there is
probably an issue with the readline library. To ensure the yt python
environment can use readline, run the following command:

.. code-block:: bash

   $ ~/yt/bin/pip install gnureadline

.. _faq-old-data:

yt seems to be plotting from old data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt does check the time stamp of the simulation so that if you
overwrite your data outputs, the new set will be read in fresh by
yt. However, if you have problems or the yt output seems to be
in someway corrupted, try deleting the ``.yt`` and
``.harray`` files from inside your data directory. If this proves to
be a persistent problem add the line:

.. code-block:: python

   from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"

to the very top of your yt script.  Turning off serialization is the default 
behavior in yt-3.0.

.. _faq-log-level:

How can I change yt's log level? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yt's default log level is ``INFO``. However, you may want less voluminous logging, especially
if you are in an IPython notebook or running a long or parallel script. On the other
hand, you may want it to output a lot more, since you can't figure out exactly what's going 
wrong, and you want to output some debugging information. The yt log level can be 
changed using the :ref:`configuration-file`, either by setting it in the 
``$HOME/.yt/config`` file:

.. code-block:: bash

   [yt]
   loglevel = 10 # This sets the log level to "DEBUG"
   
which would produce debug (as well as info, warning, and error) messages, or at runtime:

.. code-block:: python

   from yt.config import ytcfg
   ytcfg["yt","loglevel"] = "40" # This sets the log level to "ERROR"
   
which in this case would suppress everything below error messages. For reference, the numerical 
values corresponding to different log levels are:

.. csv-table:: 
   :header: Level, Numeric Value
   :widths: 10, 10
   
   ``CRITICAL``,50
   ``ERROR``,40
   ``WARNING``,30
   ``INFO``,20
   ``DEBUG``,10
   ``NOTSET``,0
   

.. _plugin-file:

Can I always load custom data objects, fields, and quantities with every dataset?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The plugin file is a means of modifying the available fields, quantities, data
objects and so on without modifying the source code of yt.  The plugin file
will be executed if it is detected.  It must be located in a ``.yt`` folder
in your home directory and be named ``my_plugins.py``:

.. code-block:: bash

   $HOME/.yt/my_plugins.py

The code in this file can add fields, define functions, define
datatypes, and on and on.  It is executed at the bottom of ``yt.mods``, and so
it is provided with the entire namespace available in the module ``yt.mods``.
For example, if I created a plugin file containing:

.. code-block:: python

   def _myfunc(field, data):
       return np.random.random(data["density"].shape)
   add_field("some_quantity", function=_myfunc, units='')

then all of my data objects would have access to the field "some_quantity".
Note that the units must be specified as a string, see
:ref:`data_selection_and_fields` for more details on units and derived fields.

.. note::

   Since the ``my_plugins.py`` is parsed inside of ``yt.mods``, you must import
   yt using ``yt.mods`` to use the plugins file.  If you import using
   ``import yt``, the plugins file will not be parsed.  You can tell that your
   plugins file is being parsed by watching for a logging message when you
   import yt.  Note that both the ``yt load`` and ``iyt`` command line entry
   points invoke ``from yt.mods import *``, so the ``my_plugins.py`` file
   will be parsed if you enter yt that way.

You can also define other convenience functions in your plugin file.  For
instance, you could define some variables or functions, and even import common
modules:

.. code-block:: python

   import os

   HOMEDIR="/home/username/"
   RUNDIR="/scratch/runs/"

   def load_run(fn):
       if not os.path.exists(RUNDIR + fn):
           return None
       return load(RUNDIR + fn)

In this case, we've written ``load_run`` to look in a specific directory to see
if it can find an output with the given name.  So now we can write scripts that
use this function:

.. code-block:: python

   from yt.mods import *

   my_run = load_run("hotgasflow/DD0040/DD0040")

And because we have imported from ``yt.mods`` we have access to the
``load_run`` function defined in our plugin file.

How do I cite yt?
^^^^^^^^^^^^^^^^^

If you use yt in a publication, we'd very much appreciate a citation!  You
should feel free to cite the `ApJS paper
<http://adsabs.harvard.edu/abs/2011ApJS..192....9T>`_ with the following BibTeX
entry: ::

   @ARTICLE{2011ApJS..192....9T,
      author = {{Turk}, M.~J. and {Smith}, B.~D. and {Oishi}, J.~S. and {Skory}, S. and 
   	{Skillman}, S.~W. and {Abel}, T. and {Norman}, M.~L.},
       title = "{yt: A Multi-code Analysis Toolkit for Astrophysical Simulation Data}",
     journal = {\apjs},
   archivePrefix = "arXiv",
      eprint = {1011.3514},
    primaryClass = "astro-ph.IM",
    keywords = {cosmology: theory, methods: data analysis, methods: numerical },
        year = 2011,
       month = jan,
      volume = 192,
       pages = {9-+},
         doi = {10.1088/0067-0049/192/1/9},
      adsurl = {http://adsabs.harvard.edu/abs/2011ApJS..192....9T},
     adsnote = {Provided by the SAO/NASA Astrophysics Data System}
   }
