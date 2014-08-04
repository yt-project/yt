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
`continuous integration server <http://tests.yt-project.org>`_ that runs the
tests with each commit to the yt version control repository.

.. _unit_testing:

Unit Testing
------------

What do Unit Tests Do
^^^^^^^^^^^^^^^^^^^^^

Unit tests are tests that operate on some small set of machinery, and verify
that the machinery works.  yt uses the `Nose
<http://nose.readthedocs.org/en/latest/>`_ framework for running unit tests.
In practice, what this means is that we write scripts that ``yield``
assertions, and Nose identifies those scripts, runs them, and verifies that the
assertions are true.

How to Run the Unit Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

One can run the unit tests very straightforwardly from any python interpreter
that can import the yt module:

.. code-block:: python

   import yt
   yt.run_nose()

If you are developing new functionality, it is sometimes more convenient to use
the Nose command line interface, ``nosetests``. You can run the unit tests
using ``nose`` by navigating to the base directory of the yt mercurial
repository and invoking ``nosetests``:

.. code-block:: bash

   $ cd $YT_HG
   $ nosetests

where ``$YT_HG`` is the path to the root of the yt mercurial repository.

If you want to specify a specific unit test to run (and not run the entire
suite), you can do so by specifying the path of the test relative to the
``$YT_HG/yt`` directory -- note that you strip off one yt more than you
normally would!  For example, if you want to run the plot_window tests, you'd
run:

.. code-block:: bash

   $ nosetests visualization/tests/test_plotwindow.py

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
* :func:`~yt.testing.amrspace` provides the ability to create AMR grid
  structures.
* :func:`~yt.testing.expand_keywords` provides the ability to iterate over
  many values for keywords.

To create new unit tests:

#. Create a new ``tests/`` directory next to the file containing the
   functionality you want to test.  Be sure to add this new directory as a
   subpackage in the setup.py script located in the directory you're adding a
   new ``tests/`` folder to.  This ensures that the tests will be deployed in
   yt source and binary distributions.
#. Inside that directory, create a new python file prefixed with ``test_`` and
   including the name of the functionality.
#. Inside that file, create one or more routines prefixed with ``test_`` that
   accept no arguments.  These should ``yield`` a set of values of the form
   ``function``, ``arguments``.  For example ``yield assert_equal, 1.0, 1.0``
   would evaluate that 1.0 equaled 1.0.
#. Use ``fake_random_ds`` to test on datasets, and be sure to test for
   several combinations of ``nproc``, so that domain decomposition can be
   tested as well.
#. Test multiple combinations of options by using the
   :func:`~yt.testing.expand_keywords` function, which will enable much
   easier iteration over options.

For an example of how to write unit tests, look at the file
``yt/data_objects/tests/test_covering_grid.py``, which covers a great deal of
functionality.

.. _answer_testing:

Answer Testing
--------------

What do Answer Tests Do
^^^^^^^^^^^^^^^^^^^^^^^

Answer tests test **actual data**, and many operations on that data, to make
sure that answers don't drift over time.  This is how we will be testing
frontends, as opposed to operations, in yt.

.. _run_answer_testing:

How to Run the Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The very first step is to make a directory and copy over the data against which
you want to test.  Currently, we test:

* ``DD0010/moving7_0010`` (available in ``tests/`` in the yt distribution)
* ``IsolatedGalaxy/galaxy0030/galaxy0030``
* ``WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030``
* ``GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300``
* ``TurbBoxLowRes/data.0005.3d.hdf5``
* ``GaussianCloud/data.0077.3d.hdf5``
* ``RadAdvect/plt00000``
* ``RadTube/plt00500``

These datasets are available at http://yt-project.org/data/.

Next, modify the file ``~/.yt/config`` to include a section ``[yt]``
with the parameter ``test_data_dir``.  Set this to point to the
directory with the test data you want to compare.  Here is an example
config file:

.. code-block:: none

   [yt]
   test_data_dir = /Users/tomservo/src/yt-data

More data will be added over time.  To run the tests, you can import the yt
module and invoke ``yt.run_nose()`` with a new keyword argument:

__ run_answer_testing_

.. code-block:: python

   import yt
   yt.run_nose(run_answer_tests=True)

If you have installed yt using ``python setup.py develop`` you can also
optionally invoke nose using the ``nosetests`` command line interface:

.. code-block:: bash

   $ cd $YT_HG
   $ nosetests --with-answer-testing

In either case, the current gold standard results will be downloaded from the
rackspace cloud and compared to what is generated locally.  The results from a
nose testing session are pretty straightforward to understand, the results for
each test are printed directly to STDOUT. If a test passes, nose prints a
period, F if a test fails, and E if the test encounters an exception or errors
out for some reason.  If you want to also run tests for the 'big' datasets,
then you can use the ``answer_big_data`` keyword argument:

.. code-block:: python

   import yt
   yt.run_nose(run_answer_tests=True, answer_big_data=True)

or, in the base directory of the yt mercurial repository:

.. code-block:: bash

   $ nosetests --with-answer-testing --answer-big-data

It's also possible to only run the answer tests for one frontend.  For example,
to run only the enzo answers tests, one can do,

.. code-block:: bash

   $ nosetests --with-answer-testing yt.frontends.enzo

How to Write Answer Tests
^^^^^^^^^^^^^^^^^^^^^^^^^

Tests can be added in the file ``yt/utilities/answer_testing/framework.py`` .
You can find examples there of how to write a test.  Here is a trivial example:

.. code-block:: python

   #!python
   class MaximumValue(AnswerTestingTest):
       _type_name = "ParentageRelationships"
       _attrs = ("field",)
       def __init__(self, ds_fn, field):
           super(MaximumValue, self).__init__(ds_fn)
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

How to Add Data to the Testing Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add data to the testing suite, first write a new set of tests for the data.
The Enzo example in ``yt/frontends/enzo/tests/test_outputs.py`` is
considered canonical.  Do these things:

* Create a new directory, ``tests`` inside the frontend's directory.

* Create a new file, ``test_outputs.py`` in the frontend's ``tests``
  directory.

* Create a new routine that operates similarly to the routines you can see
  in Enzo's outputs.

  * This routine should test a number of different fields and data objects.

  * The test routine itself should be decorated with
    ``@requires_ds(file_name)``  This decorate can accept the argument
    ``big_data`` for if this data is too big to run all the time.

  * There are ``small_patch_amr`` and ``big_patch_amr`` routines that
    you can yield from to execute a bunch of standard tests.  This is where
    you should start, and then yield additional tests that stress the
    outputs in whatever ways are necessary to ensure functionality.

  * **All tests should be yielded!**

If you are adding to a frontend that has a few tests already, skip the first
two steps.

How to Upload Answers
^^^^^^^^^^^^^^^^^^^^^

To upload answers you can execute this command:

.. code-block:: bash

   $ nosetests --with-answer-testing frontends/enzo/ --answer-store --answer-name=whatever

The current version of the gold standard can be found in the variable
``_latest`` inside ``yt/utilities/answer_testing/framework.py``  As of
the time of this writing, it is ``gold007``  Note that the name of the
suite of results is now disconnected from the dataset's name, so you
can upload multiple outputs with the same name and not collide.

To upload answers, you **must** have the package boto installed, and you
**must** have an Amazon key provided by Matt.  Contact Matt for these keys.
