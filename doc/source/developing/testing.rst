.. _testing:

Testing
=======

yt includes a testing suite which one can run on the codebase to assure that no
breaks in functionality have occurred.  This testing suite is based on the Pytest_
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
that the machinery works.  yt uses the `Pytest
<https://docs.pytest.org/en/latest/>`_ framework for running unit tests.  In
practice, what this means is that we write scripts that assert statements, and
Pytest identifies those scripts, runs them, and verifies that the assertions are
true and the code runs without crashing.

How to Run the Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

One can run the unit tests very straightforwardly from any python interpreter
that can import the yt module:

.. code-block:: python

   import yt
   yt.run_pytest()

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

To run a specific test from the plot_window test suite (for example, the ``on_off_compare`` test), you'd run:

.. code-block:: bash

   $ pytest yt/visualization/tests/test_plotwindow.py::test_on_off_compare

In general, if the test you'd like to run is a method of a class, you'd run:

.. code-block:: bash

   $ pytest path/to/test_module.py::TestClass::method_name

where ``test_module.py``, ``TestClass``, and ``method_name`` would be replaced by the desired module, class, and method names you'd like to test, respectively.

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

To create new unit tests:

#. Create a new ``tests/`` directory next to the file containing the
   functionality you want to test and add an empty ``__init__.py`` file to
   it.
#. Inside that directory, create a new python file prefixed with ``test_`` and
   including the name of the functionality.
#. Inside that file, create one or more routines prefixed with ``test_`` that
   accept no arguments (other than optional [pytest fixtures](https://docs.pytest.org/en/latest/fixture.html)). The test function should do some work that tests some
   functionality and should also verify that the results are correct using
   assert statements or functions.
#. Use ``fake_random_ds`` to test on datasets, and be sure to test for
   several combinations of ``nproc``, so that domain decomposition can be
   tested as well.
#. Test multiple combinations of options by using the
   [`pytest.mark.parametrize`](https://docs.pytest.org/en/latest/parametrize.html)
   decorator, which can be stacked and enables much easier iteration over options.
#. If using pytest fixtures, be sure to add them to a ``conftest.py`` file next to your ``test_`` file.

For an example of how to write unit tests, look at the file
``yt/data_objects/tests/test_covering_grid.py``, which covers a great deal of
functionality.

Debugging failing tests
^^^^^^^^^^^^^^^^^^^^^^^

When writing new tests, often one exposes bugs or writes a test incorrectly,
causing an exception to be raised or a failed test. To help debug issues like
this, ``pytest`` can drop into a debugger whenever a test fails or raises an
exception. This can be accomplished in one of three ways: by passing ``--pdb`` to the ``pytest`` executable, passing ``--trace`` to the ``pytest`` executable, or by inserting ``import pdb; pdb.set_trace()`` wherever you'd like to enter the debugger. The ``--pdb`` option will drop into the pdb debugger
whenever an error is raised or a failure happens. The ``--trace`` option will enter the debugger at the start of each test being run. These two options can be combined. Inside the debugger you can interactively print out variables and go up and down the call stack to determine the context for your failure or error.

.. code-block:: bash

    pytest --pdb --trace

In addition, one can debug more crudely using print statements. To do this,
you can add print statements to the code as normal. However, the test runner
will capture all print output by default. To ensure that output gets printed
to your terminal while the tests are running, pass ``-s`` to the ``pytest``
executable.

Finally, to determine which test is failing while the tests are running, it helps to run the tests in "verbose" mode. This can be done by passing the ``-v`` option to the ``pytest`` executable.

All of the above ``pytest`` options can be combined. So, for example to run
the ``TestSetWidth`` tests with verbose output, letting the output of print
statements come out on the terminal prompt, and enabling pdb debugging on errors
or test failures, one would do:

.. code-block:: bash

    $ pytest --pdb --trace -v -s yt/visualization/tests/test_plotwindow::TestSetWidth

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
   $ pytest --with-answer-testing --answer-store yt/frontends/tipsy

This command will create a set of local answers from the tipsy frontend tests
and store them in ``test_data_dir/answers``. The answers are stored in a file named ``tipsy_answers.yaml``. The structure of the answer file is:

.. code-block:: bash

   function_name:
      value_of_test_parameter_1:
         value_of_test_parameter_2:
         ...
            value_of_test_parameter_n: hash

where ``function_name`` is the name of the test function without the ``test_`` prefix. The hash is generated from the data produced by the test. Using hashes allows for simple asserts on strings, as opposed to expensive comparisons of potentially very large arrays. Additionally, it greatly reduces the storage footprint of the stored answers, since potentially very large arrays no longer need to be saved to disk.

To run the tipsy frontend's answer tests using a different yt changeset, update to that changeset, recompile if necessary, and run the tests as before.

The results from a pytest testing session are pretty straightforward to
understand, as the results for each test are printed directly to STDOUT. If you
want to also run tests for the 'big' datasets, then you will need to pass
``--answer-big-data`` to ``pytest``.  For example, to run the tests for the
OWLS frontend, do the following:

.. code-block:: bash

   $ pytest --with-answer-testing --answer-store --answer-big-data yt/frontends/owls


How to Write Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

Tests can be added in the file ``yt/utilities/answer_testing/framework.py`` .
You can find examples there of how to write a test.  The tests should be added as a method to the ``AnswerTest`` class. Here is a trivial example:

.. code-block:: python

   #!python
   def maximum_value_test(self, ds, field):
        v, c = ds.find_max(field)
        result = np.empty(4, dtype="float64")
        result[0] = v
        result[1:] = c
        return result.tostring()

What this does is calculate the location and value of the maximum of a
field.  It then puts that into the variable result and returns result's
binary data (which is used in the generation of the hash by the caller of
the test).

If we wanted to use this test on a specific frontend (Enzo, for example), we would add a method prefixed with ``test_`` to the ``TestEnzo`` class in ``yt/frontends/enzo/tests/test_outputs.py``. This method would look like:

.. code-block:: python

   @requires_file(enzo_tiny)
   def test_max_value(self):
      hashes = OrderedDict()
      hashes['max_value'] = OrderedDict()
      ds = yt.load(enzo_tiny)
      for f in ds.field_list:
         h = utils.generate_hash(self.maximum_value_test(ds, f))
         hashes['max_value_test'][f] = h
      hashes = {'test_max_value' : hashes}
      utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

What this does is load in the test data being used (in this case, the data stored in the ``enzo_tiny`` file) and then set up an ordered dictionary for storing the result of each test. The ordered dictionary should be keyed by the test name and then each test parameter, in turn. This makes it easy to identify which test and which combination of parameters causes a test to fail.

We then loop over every field in the dataset's field list and run the ``maximum_value_test`` on that field. The binary data returned by the test is then converted into a hash string by the function ``yt.utilities.answer_testing.utils.generate_hash``. This hash is then saved as the value in the nested dictionary keyed by the test name and the test parameter (the field name, in this case).

Lastly, once all of the tests have been run, we make a dictionary whose key is the  method name and whose value is the hash OrderedDict. We then call the function ``yt.utilities.answer_testing.utils.handle_hashes``, which uses the already set up answer testing directory, the name of the answer_file for the desired changeset, the generated hashes, and  the variable ``answer_store``.

``answer_store`` is set to ``False`` if the command line option ``--answer-store`` is not used and ``True`` if it is. If ``--answer-store`` is set, then ``handle_hashes`` will save the generated hashes to the ``answer_file`` file. If it is not set, it will compare the generated hashes with those assumed to be already saved in the ``answer_file`` file.

To write a new test:

* Add a method to the ``AnswerTest`` class
* Where applicable, the test should return hashable binary data

To use the new test:

* Define a new method in the desired frontend's ``test_outputs.py`` file prefixed with ``test_``
* If the test returns hashable data, set up an ``OrderedDict`` whose keys will be the test name and whose values will be another (or more than one) ``OrderedDict``. This allows for multiple tests to be called in the same method
* The keys of these nested ordered dictionaries should be the name or value of the test parameter
* Create a dictionary whose key is the method name and whose value is the hash OrderedDict. This is done because, for a given frontend, every test method writes its data to the same file. As such, this makes it easier to identify where the results from each method are in the answer file.
* Call ``utils.generate_hash`` on the binary data returned by the test
* Once all of the hashes have been generated, call ``handle_hashes`` on them to either save the generated answers or compare them to already saved answers

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
    ``@requires_ds(test_dataset_name)``. The ``test_dataset_name`` should be a string containing the path you would pass to the ``yt.load`` function. It does not need to be the full path to the dataset, since the path will be automatically prepended with the location of the test data directory.  See :ref:`configuration-file` for more information about the ``test_data-dir`` configuration option.

* There are ``small_patch_amr`` and ``big_patch_amr`` routines that you can
    use to execute a bunch of standard tests. In addition we have created
    ``sph_answer`` which is more suited for particle SPH datasets. This is where
    you should start, and then use additional tests that stress the outputs in
    whatever ways are necessary to ensure functionality.

* If the test uses a large dataset, the test routine should be decorated with ``@pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'), reason="--answer-big-data not set.")`` so that it is not run unless the ``--answer-big-data`` command line option is used

* If it's a new frontend, the test routines should all be methods of a class ``TestFrontendName`` that inherits from ``AnswerTest``.

* This class should be decorated with ``@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'), reason="--with-answer-testing not set.")`` so that it is not run unless the ``--with-answer-testing`` command line option is passed

If you are adding to a frontend that has a few tests already, skip the first
two steps.

How to Write Image Comparison Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have a number of tests designed to compare images as part of yt. We make use
of some functionality from matplotlib to automatically compare images and detect
differences, if any. Image comparison tests are used in the plotting and volume
rendering machinery.

The easiest way to use the image comparison tests is to make use of the
``generic_image_test`` method. This method takes one argument:

* The name of the saved image file

You *must* decorate your test function with ``requires_ds``, otherwise the
answer testing machinery will not be properly set up.

Here is an example test function:

.. code-block:: python

   from collections import OrderedDict
   import os
   import tempfile

   import yt
   import yt.utilities.answer_testing.framework as fw
   from yt.utilities.answer_testing import utils

   answer_file = 'myfrontend_answers.yaml'

   @pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'), reason="--with-answer-testing not set")
   class TestMyFrontend(fw.AnswerTest)
      @pytest.mark.usefixtures('temp_dir')
      @utils.requires_ds(my_ds)
      def test_my_ds():
          ds = utils.data_dir_load(my_ds)
          hashes = OrderedDict()
          hashes['generic_image'] = OrderedDict()
          for f in ds.field_list:
             tmpfd, tmpfname = tempfile.mkstemp(suffix='.png')
             os.close(tmpfd)
             plt = yt.ProjectionPlot(ds, 'z', [f])
             plt.savefig(tmpfname)
             h = utils.generate_hash(self.generic_image_test(tmpfname))
             hashes['generic_image'][f] = h
          hashes = {'test_my_ds': hashes}
          utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

Another good example of an image comparison test is the
``plot_window_attribute_test`` defined in the answer testing framework and used in
``yt/visualization/tests/test_plotwindow.py``. This test shows how a new answer
test subclass can be used to programmatically test a variety of different methods
of a complicated class using the same test class. This sort of image comparison
test is more useful if you are finding yourself writing a ton of boilerplate
code to get your image comparison test working.  The ``generic_image_test`` is
more useful if you only need to do a one-off image comparison test.

Enabling Answer Tests on Jenkins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Before any code is added to or modified in the yt codebase, each incoming
changeset is run against all available unit and answer tests on our `continuous
integration server <https://tests.yt-project.org>`_. While unit tests are
autodiscovered by `pytest`_ itself,
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

   $ pytest --with-answer-testing --answer-big-data \
      --answer-name=local_artio_000 \
      yt/frontends/artio/tests/test_outputs.py

If the answer doesn't exist on the server yet, ``pytest`` is run twice and
during first pass ``--answer-store`` is added to the command line.

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

Restricting Python Versions for Answer Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If for some reason a test can be run only for a specific version of python it is
possible to indicate this by adding a ``[py2]`` or ``[py3]`` tag. For example:

.. code-block:: yaml

   answer_tests:
      local_test_000:
         - yt/test_A.py  # [py2]
         - yt/test_B.py  # [py3]

would result in ``test_A.py`` being run only for *python2* and ``test_B.py``
being run only for *python3*.
