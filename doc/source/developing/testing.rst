.. _testing:

Testing
=======

yt includes a testing suite which one can run on the codebase to assure that no
breaks in functionality have occurred.  This testing suite is based on the Nose_
testing framework.  The suite consists of two components, unit tests and answer
tests. Unit tests confirm that an isolated piece of functionality runs without
failure for inputs with known correct outputs.  Answer tests verify the
integration and compatibility of the individual code unit by generating output
from user-visible yt functions and comparing and matching the results against
outputs of the same function produced using older versions of the yt codebase.
This ensures consistency in results as development proceeds.

.. _nosetests:

The testing suite should be run locally by developers to make sure they aren't
checking in any code that breaks existing functionality.  To further this goal,
an automatic buildbot runs the test suite after each code commit to confirm
that yt hasn't broken recently.  To supplement this effort, we also maintain a
`continuous integration server <https://tests.yt-project.org>`_ that runs the
tests with each commit to the yt version control repository.

.. _unit_testing:

Unit Testing
------------

What do Unit Tests Do
^^^^^^^^^^^^^^^^^^^^^

Unit tests are tests that operate on some small set of machinery, and verify
that the machinery works.  yt uses the `Nose
<https://nose.readthedocs.io/en/latest/>`_ framework for running unit tests.  In
practice, what this means is that we write scripts that assert statements, and
Nose identifies those scripts, runs them, and verifies that the assertions are
true and the code runs without crashing.

How to Run the Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

One can run the unit tests very straightforwardly from any python interpreter
that can import the yt module:

.. code-block:: python

   import yt

   yt.run_nose()

If you are developing new functionality, it is sometimes more convenient to use
the Nose command line interface, ``nosetests``. You can run the unit tests
using ``nose`` by navigating to the base directory of the yt git
repository and invoking ``nosetests``:

.. code-block:: bash

   $ cd $YT_GIT
   $ nosetests

where ``$YT_GIT`` is the path to the root of the yt git repository.

If you want to specify a specific unit test to run (and not run the entire
suite), you can do so by specifying the path of the test relative to the
``$YT_GIT/yt`` directory -- note that you strip off one yt more than you
normally would!  For example, if you want to run the plot_window tests, you'd
run:

.. code-block:: bash

   $ nosetests yt/visualization/tests/test_plotwindow.py

Handling yt dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

Our dependencies are specified in ``setup.cfg``. Hard dependencies are found in
``options.install_requires``, while optional dependencies are specified in
``options.extras_require``. The ``full`` target contains the specs to run our
test suite, which are intended to be as modern as possible (we don't set upper
limits to versions unless we need to). The ``minimal`` target is used to check
that we don't break backward compatibility with old versions of upstream
projects by accident. It is intended to pin stricly our minimal supported
versions. The ``test`` target specifies the tools neeed to run the tests, but
not needed by yt itself.

**Python version support.**
When a new Python version is released, it takes about
a month or two for yt to support it, since we're dependent on bigger projects
like numpy and matplotlib. We vow to follow numpy's deprecation plan regarding
our supported versions for Python and numpy, defined formally in `NEP 29
<https://numpy.org/neps/nep-0029-deprecation_policy.html>`_. However, we try to
avoid bumping our minimal requirements shortly before a yt release.

**Third party dependencies.**
However, sometimes a specific version of a project that yt depends on
causes some breakage and must be blacklisted in the tests or a more
experimental project that yt depends on optionally might change sufficiently
that the yt community decides not to support an old version of that project.

**Note.**
Some of our optional dependencies are not trivial to install and their support
may vary across platforms. To manage such issue, we currently use requirement
files in additions to ``setup.cfg``. They are found in
``tests/*requirements.txt`` and used in ``tests/ci_install.sh``.

How to Write Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^

yt provides several pieces of testing assistance, all in the ``yt.testing``
module.  Describing them in detail is somewhat outside the scope of this
document, as in some cases they belong to other packages.  However, a few come
in handy:

* :func:`~yt.testing.fake_random_ds` provides the ability to create a random
  dataset, with several fields and divided into several different
  grids, that can be operated on.
* :func:`~yt.testing.assert_equal` can operate on arrays.
* :func:`~yt.testing.assert_almost_equal` can operate on arrays and accepts a
  relative allowable difference.
* :func:`~yt.testing.assert_allclose_units` raises an error if two arrays are
  not equal up to a desired absolute or relative tolerance. This wraps numpy's
  assert_allclose to correctly verify unit consistency as well.
* :func:`~yt.testing.amrspace` provides the ability to create AMR grid
  structures.
* :func:`~yt.testing.expand_keywords` provides the ability to iterate over
  many values for keywords.

To create new unit tests:

#. Create a new ``tests/`` directory next to the file containing the
   functionality you want to test and add an empty ``__init__.py`` file to
   it.
#. Inside that directory, create a new python file prefixed with ``test_`` and
   including the name of the functionality.
#. Inside that file, create one or more routines prefixed with ``test_`` that
   accept no arguments. The test function should do some work that tests some
   functionality and should also verify that the results are correct using
   assert statements or functions.
#. Tests can ``yield`` a tuple of the form ``function``, ``argument_one``,
   ``argument_two``, etc.  For example ``yield my_test, 'banana', 2.0`` would be
   captured by nose and the ``my_test`` function will be run with the provided
   arguments.
#. Use ``fake_random_ds`` to test on datasets, and be sure to test for
   several combinations of ``nproc``, so that domain decomposition can be
   tested as well.
#. Test multiple combinations of options by using the
   :func:`~yt.testing.expand_keywords` function, which will enable much
   easier iteration over options.

For an example of how to write unit tests, look at the file
``yt/data_objects/tests/test_covering_grid.py``, which covers a great deal of
functionality.

Debugging failing tests
^^^^^^^^^^^^^^^^^^^^^^^

When writing new tests, often one exposes bugs or writes a test incorrectly,
causing an exception to be raised or a failed test. To help debug issues like
this, ``nose`` can drop into a debugger whenever a test fails or raises an
exception. This can be accomplished by passing ``--pdb`` and ``--pdb-failures``
to the ``nosetests`` executable. These options will drop into the pdb debugger
whenever an error is raised or a failure happens, respectively. Inside the
debugger you can interactively print out variables and go up and down the call
stack to determine the context for your failure or error.

.. code-block:: bash

    nosetests --pdb --pdb-failures

In addition, one can debug more crudely using print statements. To do this,
you can add print statements to the code as normal. However, the test runner
will capture all print output by default. To ensure that output gets printed
to your terminal while the tests are running, pass ``-s`` to the ``nosetests``
executable.

Lastly, to quickly debug a specific failing test, it is best to only run that
one test during your testing session. This can be accomplished by explicitly
passing the name of the test function or class to ``nosetests``, as in the
following example:

.. code-block:: bash

    $ nosetests yt.visualization.tests.test_plotwindow:TestSetWidth

This nosetests invocation will only run the tests defined by the
``TestSetWidth`` class.

Finally, to determine which test is failing while the tests are running, it helps
to run the tests in "verbose" mode. This can be done by passing the ``-v`` option
to the ``nosetests`` executable.

All of the above ``nosetests`` options can be combined. So, for example to run
the ``TestSetWidth`` tests with verbose output, letting the output of print
statements come out on the terminal prompt, and enabling pdb debugging on errors
or test failures, one would do:

.. code-block:: bash

    $ nosetests --pdb --pdb-failures -v -s yt.visualization.tests.test_plotwindow:TestSetWidth

.. _answer_testing:

Answer Testing
--------------

What do Answer Tests Do
^^^^^^^^^^^^^^^^^^^^^^^

Answer tests test **actual data**, and many operations on that data, to make
sure that answers don't drift over time.  This is how we test frontends, as
opposed to operations, in yt.

.. _run_answer_testing:

How to Run the Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The very first step is to make a directory and copy over the data against which
you want to test. Next, add the config parameter ``test_data_dir`` pointing to
directory with the test data you want to test with, e.g.:

.. code-block:: bash

   $ yt config set yt test_data_dir /Users/tomservo/src/yt-data

We use a number of real-world datasets for the tests that must be downloaded and
unzipped in the ``test_data_dir`` path you have set. The test datasets, can be
downloaded from https://yt-project.org/data/. We do not explicitly list the
datasets we use in the tests here because the list of necessary datasets changes
regularly, instead you should take a look at the tests you would like to run and
make sure that the necessary data files are downloaded before running the tests.

To run the answer tests, you must first generate a set of test answers locally
on a "known good" revision, then update to the revision you want to test, and
run the tests again using the locally stored answers.

Let's focus on running the answer tests for a single frontend. It's possible to
run the answer tests for **all** the frontends, but due to the large number of
test datasets we currently use this is not normally done except on the yt
project's contiguous integration server.

.. code-block:: bash

   $ cd $YT_GIT
   $ nosetests --with-answer-testing --local --local-dir $HOME/Documents/test --answer-store --answer-name=local-tipsy yt.frontends.tipsy

This command will create a set of local answers from the tipsy frontend tests
and store them in ``$HOME/Documents/test`` (this can but does not have to be the
same directory as the ``test_data_dir`` configuration variable defined in your
``~/.config/yt/yt.toml`` file) in a file named ``local-tipsy``. To run the tipsy
frontend's answer tests using a different yt changeset, update to that
changeset, recompile if necessary, and run the tests using the following
command:

.. code-block:: bash

   $ nosetests --with-answer-testing --local --local-dir $HOME/Documents/test --answer-name=local-tipsy yt.frontends.tipsy

The results from a nose testing session are pretty straightforward to
understand, the results for each test are printed directly to STDOUT.  If a test
passes, nose prints a period, F if a test fails, and E if the test encounters an
exception or errors out for some reason.  Explicit descriptions for each test
are also printed if you pass ``-v`` to the ``nosetests`` executable.  If you
want to also run tests for the 'big' datasets, then you will need to pass
``--answer-big-data`` to ``nosetests``.  For example, to run the tests for the
OWLS frontend, do the following:

.. code-block:: bash

   $ nosetests --with-answer-testing --local --local-dir $HOME/Documents/test --answer-store --answer-big-data yt.frontends.owls


How to Write Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

Tests can be added in the file ``yt/utilities/answer_testing/framework.py`` .
You can find examples there of how to write a test.  Here is a trivial example:

.. code-block:: python

   #!python
   class MaximumValueTest(AnswerTestingTest):
       _type_name = "MaximumValue"
       _attrs = ("field",)

       def __init__(self, ds_fn, field):
           super().__init__(ds_fn)
           self.field = field

       def run(self):
           v, c = self.ds.find_max(self.field)
           result = np.empty(4, dtype="float64")
           result[0] = v
           result[1:] = c
           return result

       def compare(self, new_result, old_result):
           assert_equal(new_result, old_result)

What this does is calculate the location and value of the maximum of a
field.  It then puts that into the variable result, returns that from
``run`` and then in ``compare`` makes sure that all are exactly equal.

To write a new test:

* Subclass ``AnswerTestingTest``
* Add the attributes ``_type_name`` (a string) and ``_attrs``
  (a tuple of strings, one for each attribute that defines the test --
  see how this is done for projections, for instance)
* Implement the two routines ``run`` and ``compare``  The first
  should return a result and the second should compare a result to an old
  result.  Neither should yield, but instead actually return.  If you need
  additional arguments to the test, implement an ``__init__`` routine.
* Keep in mind that *everything* returned from ``run`` will be stored.  So if
  you are going to return a huge amount of data, please ensure that the test
  only gets run for small data.  If you want a fast way to measure something as
  being similar or different, either an md5 hash (see the grid values test) or
  a sum and std of an array act as good proxies.  If you must store a large
  amount of data for some reason, try serializing the data to a string
  (e.g. using ``numpy.ndarray.dumps``), and then compressing the data stream
  using ``zlib.compress``.
* Typically for derived values, we compare to 10 or 12 decimal places.
  For exact values, we compare exactly.

How To Write Answer Tests for a Frontend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add a new frontend answer test, first write a new set of tests for the data.
The Enzo example in ``yt/frontends/enzo/tests/test_outputs.py`` is
considered canonical.  Do these things:

* Create a new directory, ``tests`` inside the frontend's directory.

* Create a new file, ``test_outputs.py`` in the frontend's ``tests``
  directory.

* Create a new routine that operates similarly to the routines you can see
  in Enzo's output tests.

  * This routine should test a number of different fields and data objects.

  * The test routine itself should be decorated with
    ``@requires_ds(test_dataset_name)``. This decorator can accept the
    argument ``big_data=True`` if the test is expensive. The
    ``test_dataset_name`` should be a string containing the path you would pass
    to the ``yt.load`` function. It does not need to be the full path to the
    dataset, since the path will be automatically prepended with the location of
    the test data directory.  See :ref:`configuration-file` for more information
    about the ``test_data-dir`` configuration option.

  * There are ``small_patch_amr`` and ``big_patch_amr`` routines that you can
    yield from to execute a bunch of standard tests. In addition we have created
    ``sph_answer`` which is more suited for particle SPH datasets. This is where
    you should start, and then yield additional tests that stress the outputs in
    whatever ways are necessary to ensure functionality.

If you are adding to a frontend that has a few tests already, skip the first
two steps.

How to Write Image Comparison Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have a number of tests designed to compare images as part of yt. We make use
of some functionality from matplotlib to automatically compare images and detect
differences, if any. Image comparison tests are used in the plotting and volume
rendering machinery.

The easiest way to use the image comparison tests is to make use of the
``GenericImageTest`` class. This class takes three arguments:

* A dataset instance (e.g. something you load with ``yt.load`` or
  ``data_dir_load``)
* A function the test machinery can call which will save an image to disk. The
  test class will then find any images that get created and compare them with the
  stored "correct" answer.
* An integer specifying the number of decimal places to use when comparing
  images. A smaller number of decimal places will produce a less stringent test.
  Matplotlib uses an L2 norm on the full image to do the comparison tests, so
  this is not a pixel-by-pixel measure, and surprisingly large variations will
  still pass the test if the strictness of the comparison is not high enough.

You *must* decorate your test function with ``requires_ds``, otherwise the
answer testing machinery will not be properly set up.

Here is an example test function:

.. code-block:: python

   from matplotlib import pyplot as plt

   from yt.utilities.answer_testing.framework import (
       GenericImageTest,
       data_dir_load,
       requires_ds,
   )


   @requires_ds(my_ds)
   def test_my_ds():
       ds = data_dir_load(my_ds)

       def create_image(filename_prefix):
           plt.plot([1, 2], [1, 2])
           plt.savefig("%s_lineplot" % filename_prefix)

       test = GenericImageTest(ds, create_image, 12)

       # this ensures the test has a unique key in the
       # answer test storage file
       test.prefix = "my_unique_name"

       # this ensures a nice test name in nose's output
       test_my_ds.__name__ = test.description

       yield test

.. note:: The inner function ``create_image`` can create any number of images,
   as long as the corresponding filenames conform to the prefix.

Another good example of an image comparison test is the
``PlotWindowAttributeTest`` defined in the answer testing framework and used in
``yt/visualization/tests/test_plotwindow.py``. This test shows how a new answer
test subclass can be used to programmatically test a variety of different methods
of a complicated class using the same test class. This sort of image comparison
test is more useful if you are finding yourself writing a ton of boilerplate
code to get your image comparison test working.  The ``GenericImageTest`` is
more useful if you only need to do a one-off image comparison test.

Enabling Answer Tests on Jenkins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Before any code is added to or modified in the yt codebase, each incoming
changeset is run against all available unit and answer tests on our `continuous
integration server <https://tests.yt-project.org>`_. While unit tests are
autodiscovered by `nose <https://nose.readthedocs.io/en/latest/>`_ itself,
answer tests require definition of which set of tests constitute to a given
answer. Configuration for the integration server is stored in
*tests/tests.yaml* in the main yt repository:

.. code-block:: yaml

   answer_tests:
      local_artio_000:
         - yt/frontends/artio/tests/test_outputs.py
   # ...
   other_tests:
      unittests:
         - '-v'
         - '-s'

Each element under *answer_tests* defines answer name (*local_artio_000* in above
snippet) and specifies a list of files/classes/methods that will be validated
(*yt/frontends/artio/tests/test_outputs.py* in above snippet). On the testing
server it is translated to:

.. code-block:: bash

   $ nosetests --with-answer-testing --local --local-dir ... --answer-big-data \
      --answer-name=local_artio_000 \
      yt/frontends/artio/tests/test_outputs.py

If the answer doesn't exist on the server yet, ``nosetests`` is run twice and
during first pass ``--answer-store`` is added to the commandline.

Updating Answers
~~~~~~~~~~~~~~~~

In order to regenerate answers for a particular set of tests it is sufficient to
change the answer name in *tests/tests.yaml* e.g.:

.. code-block:: diff

   --- a/tests/tests.yaml
   +++ b/tests/tests.yaml
   @@ -25,7 +25,7 @@
        - yt/analysis_modules/halo_finding/tests/test_rockstar.py
        - yt/frontends/owls_subfind/tests/test_outputs.py

   -  local_owls_000:
   +  local_owls_001:
        - yt/frontends/owls/tests/test_outputs.py

      local_pw_000:

would regenerate answers for OWLS frontend.

When adding tests to an existing set of answers (like ``local_owls_000`` or ``local_varia_000``),
it is considered best practice to first submit a pull request adding the tests WITHOUT incrementing
the version number. Then, allow the tests to run (resulting in "no old answer" errors for the missing
answers). If no other failures are present, you can then increment the version number to regenerate
the answers. This way, we can avoid accidentally covering up test breakages.

Adding New Answer Tests
~~~~~~~~~~~~~~~~~~~~~~~

In order to add a new set of answer tests, it is sufficient to extend the
*answer_tests* list in *tests/tests.yaml* e.g.:

.. code-block:: diff

   --- a/tests/tests.yaml
   +++ b/tests/tests.yaml
   @@ -60,6 +60,10 @@
        - yt/analysis_modules/absorption_spectrum/tests/test_absorption_spectrum.py:test_absorption_spectrum_non_cosmo
        - yt/analysis_modules/absorption_spectrum/tests/test_absorption_spectrum.py:test_absorption_spectrum_cosmo

   +  local_gdf_000:
   +    - yt/frontends/gdf/tests/test_outputs.py
   +
   +
    other_tests:
      unittests:
