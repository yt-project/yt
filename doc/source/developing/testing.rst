.. _testing:

Testing
=======

yt includes a testing suite that one can run on the code base to ensure that no
breaks in functionality have occurred. This suite is based on the `pytest <https://docs.pytest.org/en/stable/>`_ testing framework and consists of two types of tests:

* Unit tests. These make sure that small pieces of code run (or fail) as intended in predictable contexts. See :ref:`unit_testing`.

* Answer tests. These generate outputs from the user-facing yt API and compare them against the outputs produced using an older, "known good", version of yt. See :ref:`answer_testing`.

These tests ensure consistency in results as development proceeds.

We recommend that developers run tests locally on changed features when developing to help ensure that the new code does not break any existing functionality. To further this goal and ensure that changes do not propagate errors or have unintentional consequences on the rest of the codebase, the full test suite is run through our continuous integration (CI) servers. CI is run on push on open pull requests on a variety of computational platforms using Github Actions and a `continuous integration server <https://tests.yt-project.org>`_ at the University of Illinois. The full test suite may take several hours to run, so we do not recommend running it locally.

.. _unit_testing:

Unit Testing
------------

What Do Unit Tests Do
^^^^^^^^^^^^^^^^^^^^^

Unit tests are tests that operate on some small piece of code and verify
that it behaves as intended. In
practice, this means that we write simple assertions (or ``assert`` statements) about results and then let pytest go through them. A test is considered a success when no error (in particular ``AssertionError``) occurs.

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

from the yt source code root directory.

Additionally, if you only want to run a specific test in a test file (rather than all of the tests contained in the file), such as, ``test_all_fields`` in ``plot_window.py``, you can do so by running:


.. code-block:: bash

    $ pytest yt/visualization/tests/test_plotwindow.py::test_all_fields

from the yt source code rood directory

See the pytest documentation for more on how to `invoke pytest <https://docs.pytest.org/en/stable/usage.html?highlight=invocation>`_ and `select tests <https://docs.pytest.org/en/stable/usage.html#specifying-tests-selecting-tests>`_.


Unit Test Tools
^^^^^^^^^^^^^^^

yt provides several helper functions and decorators to write unit tests. These tools all reside in the ``yt.testing``
module.  Describing them all in detail is outside the scope of this
document, as in some cases they belong to other packages.

How To Write New Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create new unit tests:

#. Create a new ``tests/`` directory next to the file containing the
   functionality you want to test and add an empty ``__init__.py`` file to
   it.
   #. If a ``tests/`` directory already exists, there is no need to create a new one.
#. Inside this new ``tests/`` directory, create a new python file prefixed with ``test_`` and
   including the name of the functionality or source file being tested.
   #. If a file testing the functionality you're interested in already exists, please add your tests to the existing there.
#. Inside this new ``test_`` file, create functions  prefixed with ``test_`` that
   accept no arguments.
#. Each test function should do some work that tests some
   functionality and should also verify that the results are correct using
   assert statements or functions.
#. If a dataset is needed, use ``fake_random_ds``, ``fake_amr_ds``, or ``fake_particle_ds`` (the former two of which have support for particles that may be utilized) and be sure to test for
   several combinations of ``nproc`` so that domain decomposition can be
   tested as well.
#. To iterate over multiple options, or combinations of options,
   use the `@pytest.mark.parametrize decorator <https://docs.pytest.org/en/6.2.x/parametrize.html#parametrizemark>`_.

For an example of how to write unit tests, look at the file
``yt/data_objects/tests/test_covering_grid.py``, which covers a great deal of
functionality.

Debugging Failing Tests
^^^^^^^^^^^^^^^^^^^^^^^
When writing new tests, one often exposes bugs or writes a test incorrectly,

causing an exception to be raised or a failed test. To help debug issues like
this, ``pytest`` can `drop into a debugger <https://docs.pytest.org/en/6.2.x/usage.html#dropping-to-pdb-python-debugger-on-failures>`_ whenever a test fails or raises an
exception.

In addition, one can debug more crudely using print statements. To do this,
you can add print statements to the code as normal. However, the test runner
will capture all print output by default. To ensure that output gets printed
to your terminal while the tests are running, pass ``-s`` (which will disable stdout and stderr capturing) to the ``pytest``
executable.

.. code-block:: bash

    $ pytest -s

Lastly, to quickly debug a specific failing test, it is best to only run that
one test during your testing session. This can be accomplished by explicitly
passing the name of the test function or class to ``pytest``, as in the
following example:

.. code-block:: bash

    $ pytest yt/visualization/tests/test_plotwindow.py::TestSetWidth

This pytest invocation will only run the tests defined by the
``TestSetWidth`` class. See the `pytest documentation <https://docs.pytest.org/en/6.2.x/usage.html>`_ for more on the various ways to invoke pytest.

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

    $ pytest yt/visualization/tests/test_plotwindow.py::TestSetWidth -v -s --pdb

More pytest options can be found by using the ``--help`` flag

.. code-block:: bash

    $ pytest --help

.. _answer_testing:

Answer Testing
--------------
.. note::
    This section documents answer tests run with ``pytest``. The plan is to
    switch to using ``pytest`` for answer tests at some point in the future,
    but currently (July 2024), answer tests are still implemented and run with
    ``nose``. We generally encourage developers to use ``pytest`` for any new
    tests, but if you need to change or update one of the older ``nose``
    tests, or are, e.g., writing a new frontend,
    an `older version of this documentation <https://yt-project.org/docs/4.0.0/developing/testing.html#answer-testing>`_
    decribes how the ``nose`` tests work.

.. note::
    Given that nose never had support for Python 3.10 (which as of yt 4.4 is our
    oldest supported version), it is necessary to patch it to get tests running.
    This is the command we run on CI to this end
    ``find .venv/lib/python3.10/site-packages/nose -name '*.py' -exec sed -i -e s/collections.Callable/collections.abc.Callable/g '{}' ';'``

What Do Answer Tests Do
^^^^^^^^^^^^^^^^^^^^^^^

Answer tests use `actual data <https://yt-project.org/data/>`_ to test reading, writing, and various manipulations of that data. Answer tests are how we test frontends, as opposed to operations, in yt.

In order to ensure that each of these operations are performed correctly, we store gold standard versions of yaml files called answer files. More generally, an answer file is a yaml file containing the results of having run the answer tests, which can be compared to a reference, enabling us to control that results do not drift over time.

.. _run_answer_testing:

How to Run the Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to run the answer tests locally:

* Create a directory to hold the data you'll be using for the answer tests you'll be writing or the answer tests you'll be running. This directory should be outside the yt git repository in a place that is logical to where you would normally store data.

* Add folders of the required data to this directory. Other yt data, such as ``IsolatedGalaxy``, can be downloaded to this directory as well.

* Tell yt where it can find the data. This is done by setting the config parameter ``test_data_dir`` to the path of the directory with the test data downloaded from https://yt-project.org/data/. For example,

.. code-block:: bash

   $ yt config set yt test_data_dir /Users/tomservo/src/yt-data

this should only need to be done once (unless you change where you're storing the data, in which case you'll need to repeat this step so yt looks in the right place).

* Generate or obtain a set of gold standard answer files. In order to generate gold standard answer files, wwitch to a "known good" version of yt and then run the answer tests as described below. Once done, switch back to the version of yt you wish to test.
* Now you're ready to run the answer tests!

As an example, let's focus on running the answer tests for the tipsy frontend. Let's also assume that we need to generate a gold standard answer file. To do this, we first switch to a "known good" version of yt and run the following command from the top of the yt git directory (i.e., ``$YT_GIT``) in order to generate the gold standard answer file:

.. note::
    It's possible to run the answer tests for **all** the frontends, but due to the large number of test datasets we currently use this is not normally done except on the yt project's contiguous integration server.

.. code-block:: bash

   $ cd $YT_GIT
   $ pytest --with-answer-testing --answer-store --local-dir="$HOME/Documents/test" -k "TestTipsy"

The ``--with-answer-testing`` tells pytest that we want to run answer tests. Without this option, the unit tests will be run instead of the answer tests. The ``--answer-store`` option tells pytest to save the results produced by each test to a local gold standard answer file. Omitting this option is how we tell pytest to compare the results to a gold standard. The ``--local-dir`` option specifies where the gold standard answer file will be saved (or is already located, in the case that ``--answer-store`` is omitted). The ``-k`` option tells pytest that we only want to run tests whose name matches the given pattern.

.. note::
    The path specified by ``--local-dir`` can, but does not have to be, the same directory as the ``test_data_dir`` configuration variable. It is best practice to keep the data that serves as input to yt separate from the answers produced by yt's tests, however.

.. note::
    The value given to the `-k` option (e.g., `"TestTipsy"`) is the name of the class containing the answer tests. You do not need to specify the path.

The newly generated gold standard answer file will be named ``tipsy_answers_xyz.yaml``, where ``xyz`` denotes the version number of the gold standard answers. The answer version number is determined by the ``answer_version`` attribute of the class being tested (e.g., ``TestTipsy.answer_version``).

.. note::
    Changes made to yt sometimes result in known, expected changes to the way certain operations behave. This necessitates updating the gold standard answer files. This process is accomplished by changing the version number specified in each answer test class (e.g., ``TestTipsy.answer_version``). The answer version for each test class can be found as the attribute `answer_version` of that class.

Once the gold standard answer file has been generated we switch back to the version of yt we want to test, recompile if necessary, and run the tests using the following command:

.. code-block:: bash

   $ pytest --with-answer-testing --local-dir="$HOME/Documents/test" -k "TestTipsy"

The result of each test is printed to STDOUT. If a test passes, pytest prints a period. If a test fails, encounters an
exception, or errors out for some reason, then an F is printed.  Explicit descriptions for each test
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


.. _update_image_tests:

Creating and Updating Image Baselines for pytest-mpl Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use `pytest-mpl <https://pytest-mpl.readthedocs.io/en/latest/>`_ for image comparison
tests. These tests take the form of functions, which must be decorated with
``@pytest.mark.mpl_image_compare`` and return a ``matplotlib.figure.Figure`` object.

The collection of reference images is kept as git submodule in ``tests/pytest_mpl_baseline/``.

There are 4 situations where updating reference images may be necessary

- adding new tests
- bugfixes
- intentional change of style in yt
- old baseline fails with a new version of matplotlib, but changes are not noticeable to the human eye

The process of updating images is the same in all cases. It involves opening two Pull Requests (PR)
that we'll number PR1 and PR2.

1. open a Pull Request (PR1) to yt's main repo with the code changes
2. wait for tests jobs to complete
3. go to the "Checks" tab on the PR page (``https://github.com/yt-project/yt/pull/<PR number>/checks``)
4. if all tests passed, you're done !
5. if tests other than image tests failed, fix them, and go back to step 2.
  Otherwise, if only image tests failed, navigate to the "Build and Tests" job summary page.
6. at the bottom of the page, you'll find "Artifacts".
   Download ``yt_pytest_mpl_results.zip``, unzip it and open ``fig_comparison.html`` therein;
   This document is an interactive report of the test job.
   Inspect failed tests results and verify that any differences are either intended or insignificant.
   If they are not, fix the code and go back to step 2
7. clone ``https://github.com/yt-project/yt_pytest_mpl_baseline.git`` and unzip the new baseline
8. Download the other artifact (``yt_pytest_mpl_new_baseline.zip``),
   unzip it within your clone of ``yt_pytest_mpl_baseline``.
9. create a branch, commit all changes, and open a Pull Request (PR2) to ``https://github.com/yt-project/yt_pytest_mpl_baseline``
   (PR2 should link to PR1)
10. wait for this second PR to be merged
11. Now it's time to update PR1: navigate back to your local copy of ``yt``'s main repository.
12. run the following commands

.. code-block:: bash

  $ git submodule update --init
  $ cd tests/pytest_mpl_baseline
  $ git checkout main
  $ git pull
  $ cd ../
  $ git add pytest_mpl_baseline
  $ git commit -m "update image test baseline"
  $ git push

13. go back to step 2. This time everything should pass. If not, ask for help !

.. note::
    Though it is technically possible to (re)generate reference images locally, it is
    best not to, because at a pixel level, matplotlib's behaviour is platform-dependent.
    By letting CI runners generate images, we ensure pixel-perfect comparison is possible
    in CI, which is where image comparison tests are most often run.


.. _deprecated_generic_image:

How to Write Image Comparison Tests (deprecated API)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
    this section describes deprecated API. New test code should follow :ref:`_update_image_tests`


Many of yt's operations involve creating and manipulating images. As such, we have a number of tests designed to compare images. These tests employ functionality from matplotlib to automatically compare images and detect
differences, if any. Image comparison tests are used in the plotting and volume
rendering machinery.

The easiest way to use the image comparison tests is to make use of the
``generic_image`` function. As an argument, this function takes a function the test machinery can call which will save an image to disk. The test will then find any images that get created and compare them with the stored "correct" answer.

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
        assert len(gi) == 1
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

Our dependencies are specified in ``pyproject.toml``. Hard dependencies are found in
``project.dependencies``, while optional dependencies are specified in
``project.optional-dependencies``. The ``full`` target contains the specs to run our
test suite, which are intended to be as modern as possible (we don't set upper
limits to versions unless we need to).

The ``test`` target specifies the tools needed to run the tests, but
not needed by yt itself.

Documentation and typechecking requirements are found in ``requirements/``,
and used in ``tests/ci_install.sh``.

**Python version support.**
We vow to follow numpy's deprecation plan regarding our supported versions for Python
and numpy, defined formally in
`NEP 29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_, but generally
support larger version intervals than recommended in this document.

**Third party dependencies.**
We attempt to make yt compatible with a wide variety of upstream software
versions.
However, sometimes a specific version of a project that yt depends on
causes some breakage and must be blacklisted in the tests or a more
experimental project that yt depends on optionally might change sufficiently
that the yt community decides not to support an old version of that project.

**Note.**
Some of our optional dependencies are not trivial to install and their support
may vary across platforms.

If you would like to add a new dependency for yt (even an optional dependency)
or would like to update a version of a yt dependency, you must edit the
``pyproject.toml`` file. For new dependencies, simply append the name of the new
dependency to the end of the file, along with a pin to the latest version
number of the package. To update a package's version, simply update the version
number in the entry for that package.

Finally, we also run a set of tests with "minimal" dependencies installed.
When adding tests that depend on an optional dependency, you can wrap the test
with the ``yt.testing.requires_module decorator`` to ensure it does not run
during the minimal dependency tests (see
``yt/frontends/amrvac/tests/test_read_amrvac_namelist.py`` for a good example).
If for some reason you need to update the listing of packages that are installed
for the "minimal" dependency tests, you will need to update
``minimal_requirements.txt``.
