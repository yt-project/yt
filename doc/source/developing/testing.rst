.. _testing:

Testing
=======

yt includes a testing suite that one can run on the code base to assure that no
breaks in functionality have occurred. This suite is based on the `pytest <https://docs.pytest.org/en/stable/>`_ testing framework and consists of two types of tests:

* Unit tests. These confirm that an isolated piece of functionality runs without failure for inputs with known correct outputs.

* Answer tests. These generate output from user-visible yt functions and compare this output against the output produced using an older, "known good", version of yt.

These tests ensure consistency in results as development proceeds.

The testing suite should be run locally by developers to make sure they aren't
checking in any code that breaks existing functionality. Additionally, as a safeguard, the test suite is automatically run after each commit that is pushed to the yt repository on github.

.. _unit_testing:

Unit Testing
------------

What Do Unit Tests Do
^^^^^^^^^^^^^^^^^^^^^

Unit tests are tests that operate on some small piece of machinery and verify
that the machinery works. In
practice, this means that we make assertions about a piece of machinery and then pytest runs the machinery, verifies that the assertions are true, and ensures that the code runs without crashing. An example of a unit test is :func:`~yt.fields.tests.test_fields.test_all_fields`.

How to Run the Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

One can run the unit tests by navigating to the base directory of the yt git
repository and invoking ``pytest``:

.. code-block:: bash

   $ cd $YT_GIT
   $ pytest

where ``$YT_GIT`` is the path to the root of the yt git repository.

If you only want to run tests in a specific file, you can do so by specifying the path of the test relative to the
``$YT_GIT/`` directory. For example, if you want to run the ``plot_window`` tests, you'd
run:

.. code-block:: bash

   $ pytest yt/visualization/tests/test_plotwindow.py

Additionally, if you only want to run, say, ``test_all_fields`` in ``plot_window.py``, you can do so by running:

.. code-block:: bash

    $ pytest yt/visualization/tests/test_plotwindow.py::test_all_fields

See this `link <https://docs.pytest.org/en/stable/usage.html?highlight=invocation>`_ for more on how to invoke pytest.


Unit Test Tools
^^^^^^^^^^^^^^^

yt provides several tools to assist with writing unit tests. These tools all reside in the ``yt.testing``
module.  Describing them all in detail is somewhat outside the scope of this
document. However, a few come
in handy:

* :func:`~yt.testing.fake_random_ds` provides the ability to create a random
  dataset, with several fields and divided into several different
  grids.
* :func:`~yt.testing.assert_equal` can operate on arrays.
* :func:`~yt.testing.assert_almost_equal` can operate on arrays and accepts a
  relative allowable difference.
* :func:`~yt.testing.assert_allclose_units` raises an error if two arrays are
  not equal up to a desired absolute or relative tolerance. This wraps numpy's
  :py:func:`numpy.testing.assert_allclose` to correctly verify unit consistency as well.
* :func:`~yt.testing.amrspace` provides the ability to create AMR grid
  structures.
* :func:`~yt.testing.expand_keywords` provides the ability to iterate over
  many values for keywords.

How To Write New Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create new unit tests:

#. Create a new ``tests/`` directory next to the file containing the
   machinery you want to test and add an empty ``__init__.py`` file to
   it.
#. Inside this new ``tests/`` directory, create a new python file prefixed with ``test_`` and
   including the name of the functionality.
#. Inside this new ``test_`` file, create one or more routines prefixed with ``test_`` that
   accept no arguments.
#. Each test function should do some work that tests some
   functionality and should also verify that the results are correct using
   assert statements or functions.
#. If a dataset is needed, use ``fake_random_ds`` and be sure to test for
   several combinations of ``nproc`` so that domain decomposition can be
   tested as well.
#. To iterate over multiple options, or combinations of options,
   use the :ref:`@pytest.mark.parametrize <@pytest.mark.parametrize>` decorator.

For an example of how to write unit tests, look at the file
``yt/data_objects/tests/test_covering_grid.py``, which covers a great deal of
functionality.

Debugging Failing Tests
^^^^^^^^^^^^^^^^^^^^^^^

When writing new tests, one often exposes bugs or writes a test incorrectly,
causing an exception to be raised or a failed test. To help debug issues like
this, ``pytest`` can drop into a debugger whenever a test fails or raises an
exception. This can be accomplished by passing ``--pdb``
to the ``pytest`` executable. Inside the
debugger you can interactively print out variables and go up and down the call
stack to determine the context for your failure or error.

.. code-block:: bash

    pytest --pdb

In addition, one can debug more crudely using print statements. To do this,
you can add print statements to the code as normal. However, the test runner
will capture all print output by default. To ensure that output gets printed
to your terminal while the tests are running, pass ``-s`` to the ``pytest``
executable.

.. code-block:: bash

    pytest -s

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

.. code-block:: bash

    $ pytest -v

All of the above ``pytest`` options can be combined. So, for example, to run
the ``TestSetWidth`` tests with verbose output, letting the output of print
statements come out on the terminal prompt, and enabling pdb debugging on errors
or test failures, one would do:

.. code-block:: bash

    $ pytest --pdb -v -s yt/visualization/tests/test_plotwindow.py::TestSetWidth

More pytest options can be found by using the ``-h`` option

.. code-block:: bash

    $ pytest -h

.. _answer_testing:

Answer Testing
--------------

What Do Answer Tests Do
^^^^^^^^^^^^^^^^^^^^^^^

Answer tests use `actual data <https://yt-project.org/data/>`_ to test reading, writing, and various manipulations of that data.

In order to ensure that each of these operations are performed correctly, yt has a set of yaml files containing the known correct results of having performed these operations. These files are called gold standard answer files. More generally, an answer file is a yaml file containing the results of having run the answer tests.

When the answer tests are run, their output is compared to the gold standard answers to ensure that the results of operating on data do not change over time.

.. _run_answer_testing:

How to Run the Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to run the answer tests locally:

* Create a directory to hold the data required by the answer tests one wishes to run

* Fill this directory with the required data

Next, yt needs to be made aware of where it can find the data. This is done by setting the config parameter ``test_data_dir`` to the
directory with the test data downloaded from https://yt-project.org/data/. For example,

.. code-block:: bash

   $ yt config set yt test_data_dir /Users/tomservo/src/yt-data

To run the answer tests locally, you must first generate a set of gold standard answer files using a "known good" version of yt. Once done, then update to the version of yt you want to test and
run the tests again so that the results can be compared to those contained in the gold standard answer files.

As an example, let's focus on running the answer tests for the tipsy frontend.

.. note::
    It's possible to run the answer tests for **all** the frontends, but due to the large number of test datasets we currently use this is not normally done except on the yt project's contiguous integration server.

.. code-block:: bash

   $ cd $YT_GIT
   $ pytest --with-answer-testing --answer-store --local-dir="$HOME/Documents/test" -k "TestTipsy"

This command will run the tipsy answer tests and, because of the presence of the ``--answer-store`` option, save the results in a local gold standard answer file in ``$HOME/Documents/test``.

.. note::
    The path specified by ``--local-dir`` can, but does not have to be, the same directory as the ``test_data_dir`` configuration variable.

The gold standard answer file will be named ``tipsy_answers_xyz.yaml``, where ``xyz`` denotes the version number of the gold standard answers.

.. note::
    The answer version number is determined by the ``answer_version`` attribute of the class being tested (e.g., ``TestTipsy.answer_version``).

.. note::
    Changes made to yt sometimes result in known, expected changes to the way certain operations behave. This necessitates updating the gold standard answer files. This process is accomplished by changing the version number specified in each answer test class (e.g., ``TestTipsy.answer_version``).

Once the gold standard answer file has been generated, update to the version of yt you wish to test, recompile if necessary, and run the tests using the following command:

.. code-block:: bash

   $ pytest --with-answer-testing --local-dir="$HOME/Documents/test" -k "TestTipsy"

The result of each test is printed to STDOUT. If a test passes, pytest prints a period. If a test fails or encounters an
exception or errors out for some reason, then an F is printed.  Explicit descriptions for each test
are also printed if you pass ``-v`` to the ``pytest`` executable. Similar to the unit tests, the ``-s`` and ``--pdb`` options can be passed, as well.


How to Write Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

To add a new answer test:

#. Create a new directory called ``tests`` inside the directory where the component you want to test resides and add an empty ``__init__.py`` file to it.

#. Create a new file in the ``tests`` directory that will hold the new answer tests. The name of the file should begin with ``test_``.

#. Create a new class whose name begins with ``Test`` (e.g., ``TestTipsy``).

#. Decorate the class with ``pytest.mark.answer_test``. This decorator is used to tell pytest which tests are answer tests.

   .. note::

      Tests that do not have this decorator are considered to be unit tests.

#. Add the following three attributes to the class: ``answer_file=None``, ``saved_hashes=None``, and ``answer_version=000``. These attributes are used by the ``hashing`` fixture (discussed below) to automate the creation of new answer files as well as facilitate the comparison to existing answers.

#. Add methods to the class that test a number of different fields and data objects.

#. If these methods are performing calculations or data manipulation, they should store the result in a ``ndarray``, if possible. This array should be be added to the ``hashes`` (see below) dictionary like so: ``self.hashes.update(<test_name>:<array>)``, where ``<test_name>`` is the name of the function from ``yt/utilities/answer_testing/answer_tests.py`` that is being used and ``<array>`` is the ``ndarray`` holding the result

If you are adding to a frontend that has tests already, simply add methods to the existing test class.

There are several things that can make the test writing process easier:

* ``yt/utilities/answer_testing/testing_utilities.py`` contains a large number of helper functions.
* Most frontends end up needing to test much of the same functionality as other frontends. As such, a list of functions that perform such work can be found in ``yt/utilities/answer_testing/answer_tests.py``.
* `Fixtures <https://docs.pytest.org/en/stable/fixture.html>`_! You can find the set of fixtures that have already been built for yt in ``$YT_GIT/conftest.py``. If you need/want to add additional fixtures, please add them there.
* The `parametrize decorator <https://docs.pytest.org/en/stable/example/parametrize.html?highlight=parametrizing%20tests>`_ is extremely useful for performing iteration over various combinations of test parameters. It should be used whenever possible.
    * The use of this decorator allows pytest to write the names and values of the test parameters to the generated answer files, which can make debugging failing tests easier, since one can easily see exactly which combination of parameters were used for a given test.
    * It is also possible to employ the ``requires_ds`` decorator to ensure that a test does not run unless a specific dataset is found, but not necessary. If the dataset is parametrized over, then the ``ds`` fixture found in the root ``conftest.py`` file performs the same check and marks the test as failed if the dataset isn't found.

Here is what a minimal example might look like for a new frontend:

.. code-block:: python

    # Content of yt/frontends/new_frontend/tests/test_outputs.py
    import pytest

    from yt.utilities.answer_testing.answer_tests import field_values

    # Parameters to test with
    ds1 = "my_first_dataset"
    ds2 = "my_second_dataset"
    field1 = ("Gas", "Density")
    field2 = ("Gas", "Temperature")
    obj1 = None
    obj2 = ("sphere", ("c", (0.1, "unitary")))

    @pytest.mark.answer_test
    class TestNewFrontend:
        answer_file = None
        saved_hashes = None
        answer_version = "000"

        @pytest.mark.usefixtures("hashing")
        @pytest.mark.parametrize("ds", [ds1, ds2], indirect=True)
        @pytest.mark.parametrize("field", [field1, field2], indirect=True)
        @pytest.mark.parametrize("dobj", [obj1, obj2], indirect=True)
        def test_fields(self, ds, field, dobj):
            self.hashes.update({"field_values": field_values(ds, field, dobj)})

Answer test examples can be found in ``yt/frontends/enzo/tests/test_outputs.py``.

How to Write Image Comparison Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many of yt's operations involve creating and manipulating images. As such, we have a number of tests designed to compare images. These tests employ functionality from matplotlib to automatically compare images and detect
differences, if any. Image comparison tests are used in the plotting and volume
rendering machinery.

The easiest way to use the image comparison tests is to make use of the
``generic_image`` function. As an argument, this function takes a function the test machinery can call which will save an image to disk. The test will then find any images that get created and compare them with the  stored "correct" answer.

Here is an example test function (from ``yt/visualization/tests/test_raw_field_slices.py``):

.. code-block:: python

    import pytest

    import yt
    from yt.utilities.answer_testing.answer_tests import generic_image
    from yt.utilities.answer_testing.testing_utilities import data_dir_load, requires_ds

    # Test data
    raw_fields = "Laser/plt00015"

    def compare(ds, field):
        def slice_image(im_name):
            sl = yt.SlicePlot(ds, "z", field)
            sl.set_log("all", False)
            image_file = sl.save(im_name)
            return image_file

        gi = generic_image(slice_image)
        # generic_image returns a list. In this case, there's only one entry,
        # which is a np array with the data we want
        return gi[0]

    @pytest.mark.answer_test
    @pytest.mark.usefixtures("temp_dir")
    class TestRawFieldSlices:
        answer_file = None
        saved_hashes = None
        answer_version = "000"

        @pytest.mark.usefixtures("hashing")
        @requires_ds(raw_fields)
        def test_raw_field_slices(self, field):
            ds = data_dir_load(raw_fields)
            gi = compare(ds, field)
            self.hashes.update({"generic_image": gi})

.. note::
    The inner function ``slice_image`` can create any number of images, as long as the corresponding filenames conform to the prefix.

Another good example of an image comparison test is the
``plot_window_attribute`` defined in the ``yt/utilities/answer_testing/answer_tests.py`` and used in
``yt/visualization/tests/test_plotwindow.py``. This sort of image comparison
test is more useful if you are finding yourself writing a ton of boilerplate
code to get your image comparison test working.  The ``generic_image`` function is
more useful if you only need to do a one-off image comparison test.

Updating Answers
~~~~~~~~~~~~~~~~

In order to regenerate answers for a particular set of tests it is sufficient to
change the ``answer_version`` attribute in the desired test class.

When adding tests to an existing set of answers (like ``local_owls_000.yaml`` or ``local_varia_000.yaml``),
it is considered best practice to first submit a pull request adding the tests WITHOUT incrementing
the version number. Then, allow the tests to run (resulting in "no old answer" errors for the missing
answers). If no other failures are present, you can then increment the version number to regenerate
the answers. This way, we can avoid accidentally covering up test breakages.

.. _handling_dependencies:

Handling yt Dependencies
------------------------

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

Finally, we also run a set of tests with "minimal" dependencies installed. When adding tests that depend on an optional dependency, you can wrap the test with the ``yt.testing.requires_module decorator`` to ensure it does not run during the minimal dependency tests (see yt/frontends/amrvac/tests/test_read_amrvac_namelist.py for a good example). If for some reason you need to update the listing of packages that are installed for the "minimal" dependency tests, you will need to edit ``tests/test_minimal_requirements.txt``.
