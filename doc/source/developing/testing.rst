.. _testing:

Testing
=======

yt includes a testing suite that one can run on the codebase to assure that no
breaks in functionality have occurred.  This testing suite is based on the
`Pytest <https://docs.pytest.org/en/stable/>`_
testing framework.  The suite consists of two components, unit tests and answer
tests. Unit tests confirm that an isolated piece of functionality runs without
failure for inputs with known correct outputs.  Answer tests verify the
integration and compatibility of the individual code unit by generating output
from user-visible yt functions and comparing and matching the results against
outputs of the same function produced using older versions of the yt codebase.
This ensures consistency in results as development proceeds.

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
that the machinery works. In practice, what this means is that we write scripts
that assert statements, and
Pytest identifies those scripts, runs them, and verifies that the assertions are
true and the code runs without crashing.

How to Run the Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

One can run the unit tests very straightforwardly from the root yt repository
using a python interpreter:

.. code-block:: python

   import pytest
   pytest.main()

If you are developing new functionality, it is sometimes more convenient to use
the Pytest command line interface, ``pytest``. You can run the unit tests
using ``pytest`` by navigating to the base directory of the yt git
repository and invoking ``pytest``:

.. code-block:: bash

   $ cd $YT_GIT
   $ pytest

where ``$YT_GIT`` is the path to the root of the yt git repository.

If you want to specify a specific unit test to run (and not run the entire
suite), you can do so by specifying the path of the test relative to the
``$YT_GIT/yt`` directory -- note that you strip off one yt more than you
normally would!  For example, if you want to run the plot_window tests, you'd
run:

.. code-block:: bash

   $ pytest yt/visualization/tests/test_plotwindow.py

Handling yt dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

We attempt to make yt compatible with a wide variety of upstream software
versions. However, sometimes a specific version of a project that yt depends on
causes some breakage and must be blacklisted in the tests or a more
experimental project that yt depends on optionally might change sufficiently
that the yt community decides not to support an old version of that project.

To handle cases like this, the versions of upstream software projects installed
on the machines running the yt test suite are pinned to specific version
numbers that must be updated manually. This prevents breaking the yt tests when
a new version of an upstream dependency is released and allows us to manage
updates in upstream projects at our pace.

If you would like to add a new dependency for yt (even an optional dependency)
or would like to update a version of a yt dependency, you must edit the
``tests/test_requirements.txt`` file, this path is relative to the root of the
repository. This file contains an enumerated list of direct dependencies and
pinned version numbers. For new dependencies, simply append the name of the new
dependency to the end of the file, along with a pin to the latest version
number of the package. To update a package's version, simply update the version
number in the entry for that package.

Finally, we also run a set of tests with "minimal" dependencies installed. Please make sure any new tests you add that depend on an optional dependency are properly set up so that the test is not run if the dependency is not installed. If for some reason you need to update the listing of packages that are installed for the "minimal" dependency tests, you will need to edit ``tests/test_minimal_requirements.txt``.

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
this, ``pytest`` can drop into a debugger whenever a test fails or raises an
exception. This can be accomplished by passing ``--pdb`` to the ``pytest``
executable. These options will drop into the pdb debugger
whenever an error is raised or a failure happens. Inside the
debugger you can interactively print out variables and go up and down the call
stack to determine the context for your failure or error.

.. code-block:: bash

    pytest --pdb

In addition, one can debug more crudely using print statements. To do this,
you can add print statements to the code as normal. However, the test runner
will capture all print output by default. To ensure that output gets printed
to your terminal while the tests are running, pass ``-s`` to the ``pytest``
executable.

Lastly, to quickly debug a specific failing test, it is best to only run that
one test during your testing session. This can be accomplished by explicitly
passing the name of the test function or class to ``pytest``, as in the
following example:

.. code-block:: bash

    $ pytest yt/visualization/tests/test_plotwindow.py::TestSetWidth

This pytest invocation will only run the tests defined by the
``TestSetWidth`` class.

Finally, to determine which test is failing while the tests are running, it helps
to run the tests in "verbose" mode. This can be done by passing the ``-v`` option
to the ``pytest`` executable.

All of the above ``pytest`` options can be combined. So, for example to run
the ``TestSetWidth`` tests with verbose output, letting the output of print
statements come out on the terminal prompt, and enabling pdb debugging on errors
or test failures, one would do:

.. code-block:: bash

    $ pytest --pdb -v -s yt/visualization/tests/test_plotwindow.py::TestSetWidth

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
   $ pytest --with-answer-testing --local-dir="~/Documents/test" --answer-store -k "TestTipsy"

This command will create a set of local answers from the tipsy frontend tests
and store them in ``~/Documents/test`` (this can but does not have to be the
same directory as the ``test_data_dir`` configuration variable defined in your
``~/.config/yt/ytrc`` file). The name of the file the answers are stored in is
specified in ``tests/tests.yaml``. To run the tipsy
frontend's answer tests using a different yt changeset, update to that
changeset, recompile if necessary, and run the tests using the following
command:

.. code-block:: bash

   $ pytest --with-answer-testing --local-dir="~/Documents/test" -k "TestTipsy"

The results from a pytest testing session are pretty straightforward to
understand, the results for each test are printed directly to STDOUT.  If a test
passes, pytest prints a period, F if a test fails, and E if the test encounters an
exception or errors out for some reason.  Explicit descriptions for each test
are also printed if you pass ``-v`` to the ``pytest`` executable.  If you
want to also run tests for the 'big' datasets, then you will need to pass
``--answer-big-data`` to ``pytest``.  For example, to run the tests for the
OWLS frontend, do the following:

.. code-block:: bash

   $ pytest --with-answer-testing --local-dir="~/Documents/test" --answer-store --answer-big-data -k "TestOwls"

Answer Formats
^^^^^^^^^^^^^^

By default, yt produces a hash of whatever result is returned from a test. This
is because, oftentimes, the test results are large arrays, so using a hash
facilitates easier comparison and much smaller, more human-readable answer
files using the ``yaml`` format.

However, there are times when you might want to actually see the result that
a test is returning, particularly during debugging. As such, yt provides the
``--answer-raw-arrays`` and ``--raw-answer-store`` command line options for
``pytest``. If both of these options are passed to ``pytest``, then, in
addition to saving a hashed version of the result, yt will also save the raw,
un-hashed version of the result, as well. These raw results are saved in an
hdf5 file in the directory specified by `--local-dir``. For example,

.. code-block:: bash

   $ pytest --with-answer-testing --answer-raw-arrays --raw-answer-store --local-dir="~/Documents/test" -k "TestTipsy"

will run the answer tests for the tipsy frontend and save the un-hased test
results in ``~/Documents/test`` with a file extension of ``.h5``.

To compare raw answers to a set of already saved raw answers, you simply omit
the ``--raw-answer-store`` option.

.. note::

   The raw answer files can get quite large!

How to Write Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

Tests can be added in the file ``yt/utilities/answer_testing/answer_tests.py`` .
You can find examples there of how to write a test.  Here is a trivial example:

.. code-block:: python

   #!python
   def maximum_value_test(ds, field):
      v, c = ds.find_max(field)
      result = np.empty(4, dtype="float64")
      result[0] = v
      result[1:] = c
      return result

What this does is calculate the location and value of the maximum of a
field.  It then puts that into the variable result and returns. If ``--answer-
store`` has been passed to ``pytest``, then a hash of ``result`` is saved.
If ``--answer-raw-arrays`` and ``--raw-answer-store`` have been passed to
``pytest`` then the un-hashed form of ``result`` is also saved.

If the answers are **not** being saved, then yt will look for an existing answer
file, read in the results, and then compare the new results to the old ones.
This is done via the ``_compare_raw_arrays`` and ``_compare_result`` functions
in ``yt/utilities/answer_testing/testing_utilities.py`` for raw answers and
hashed answers, respectively.

How To Write Answer Tests for a Frontend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add a new frontend answer test, first write a new set of tests for the data.
The Enzo example in ``yt/frontends/enzo/tests/test_outputs.py`` is
considered canonical.  Do these things:

* Create a new directory, ``tests`` inside the frontend's directory.

* Create a new file, ``test_outputs.py`` in the frontend's ``tests``
  directory.

* Create a new class with methods that operates similarly to the methods you
 can see in Enzo's ``TestEnzo`` class.

* The class should have the name ``TestFrontendName``, where ``FrontendName`` is
the name of the frontend you are writing tests for.

* The class should be decorated with ``@pytest.mark.answer_test``.

* These methods should test a number of different fields and data objects.

   * This is most easily accomplished by taking advantage of ``pytest``'s
   `parameterization features <https://docs.pytest.org/en/2.8.7/example/parametrize.html>`_

* If the methods are going to have their results saved, then the ``@pytest.mark.usefixtures("hashing")`` decorator should be used. This
fixture takes care of hashing the results, saving the results, and comparing
the results.

* If the method requires a large data file, it should be further decorated with
``@pytest.mark.big_data`` (see, for example, ``yt/frontends/owls/tests/test_outputs.py``).

How to Write Image Comparison Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have a number of tests designed to compare images as part of yt. We make use
of some functionality from matplotlib to automatically compare images and detect
differences, if any. Image comparison tests are used in the plotting and volume
rendering machinery.

The easiest way to use the image comparison tests is to make use of the
``generic_image`` function. This function takes three arguments:

* A function the test machinery can call which will save an image to disk. yt will then find any images that get created and compare them with the
  stored "correct" answer.

Here is an example test function:

.. code-block:: python

   from yt.utilities.answer_testing.answer_tests import generic_image
   from yt.utilities.answer_testing.testing_utilities import data_dir_load

   from matplotlib import pyplot as plt

   @pytest.mark.answer_test
   class TestFrontendName:
      @pytest.mark.usefixtures("hashing")
      def test_my_ds(self):
          ds = data_dir_load(my_ds)
          my_plot = plt.plot([1, 2], [1, ])

          def create_image(image_name):
            return my_plot.save(image_name)
         gi = generic_image(create_image)[0]
         self.hashes.update({"generic_image" : gi})

Another good example of an image comparison test is
``plot_window_attribute`` defined in ``answer_tests.py`` and used in
``yt/visualization/tests/test_plotwindow.py``.

Enabling Answer Tests on Jenkins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Before any code is added to or modified in the yt codebase, each incoming
changeset is run against all available unit and answer tests on our `continuous
integration server <https://tests.yt-project.org>`_. While unit tests are
autodiscovered by `pytest <https://docs.pytest.org/en/stable/>`_ itself,
answer tests require definition of which set of tests constitute to a given
answer. Configuration for the integration server is stored in
*tests/tests.yaml* in the main yt repository:

.. code-block:: yaml

   TestAdaptahop:
      - adaptahop_answers_000.yaml
      - adaptahop_answers_raw_000.h5
   # ...

where each answer test class is listed along with the file names to use for
both the hashed and raw answers.

The number in each file name refers to the answer file revision to use.


Updating Answers
~~~~~~~~~~~~~~~~

If you are trying to save answers for a test class that already has an answer
file on disk, you will get an error unless you specify the ``--force-overwrite``
option on the command line. This is **strongly** discouraged. You should almost
always instead bump the file revision number and re-run, e.g.,
``enzo_answers_000.yaml -> enzo_answers_001.yaml``.

When adding tests to an existing set of answers (like ``owls_answers_000.yaml``),
it is considered best practice to first submit a pull request adding the tests WITHOUT incrementing
the version number. Then, allow the tests to run (resulting in errors for the missing
answers). If no other failures are present, you can then increment the version number to regenerate
the answers. This way, we can avoid accidentally covering up test breakages.

Adding New Answer Tests
~~~~~~~~~~~~~~~~~~~~~~~

In order to add a new set of answer tests, it is sufficient to extend the
*answer_tests* list in *tests/tests.yaml* e.g.:

.. code-block:: yaml

   TestNewFrontend:
      - newfrontend_answers_000.yaml
      - newfrontend_answers_raw_000.h5
