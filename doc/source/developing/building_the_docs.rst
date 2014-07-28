.. _docs_build:

==========================
Building the Documentation
==========================

The yt documentation makes heavy use of the sphinx documentation automation
suite.  Sphinx, written in python, was originally created for the documentation
of the python project and has many nice capabilities for managing the
documentation of python code.

While much of the yt documentation is static text, we make heavy use of
cross-referencing with API documentation that is automatically generated at
build time by sphinx.  We also use sphinx to run code snippets (e.g. the 
cookbook and the notebooks) and embed resulting images and example data.

Quick versus full documentation builds
--------------------------------------

Building the entire set of yt documentation is a laborious task, since you 
need to have a large number of packages in order to successfully execute
and render all of the notebooks and yt recipes drawing from every corner
of the yt source.  As an quick alternative, one can do a ``quick`` build
of the documentation, which eschews the need for downloading all of these
dependencies, but it only produces the static docs.  The static docs do 
not include the cookbook outputs and the notebooks, but this is good
enough for most cases of people testing out whether or not their documentation
contributions look OK before submitting them to the yt repository.

If you want to create the full documentation locally, then you'll need
to follow the instructions for building the ``full`` docs, so that you can
dynamically execute and render the cookbook recipes, the notebooks, etc.

Building the docs (quick)
-------------------------

You will need to have the yt repository available on your computer, which
is done by default if you have yt installed.  In addition, you need a 
current version of Sphinx_ (1.1.3) documentation software installed.

In order to tell sphinx not to do all of the dynamical building, you must
set the ``$READTHEDOCS`` environment variable to be True by typing this at 
the command line:

.. code-block:: bash

   export READTHEDOCS=True  # for bash
   setenv READTHEDOCS True  # for csh

This variable is set for automated builds on the free ReadTheDocs service but
can be used by anyone to force a quick, minimal build.

Now all you need to do is execute sphinx on the yt doc source.  Go to the 
documentation directory and build the docs:

.. code-block:: bash

   cd $YT_HG/doc
   make html

This will produce an html version of the documentation locally in the 
``$YT_HG/doc/build/html`` directory.  You can now go there and open
up ``index.html`` or whatever file you wish in your web browser.

Building the docs (full)
------------------------

As alluded to earlier, building the full documentation is a bit more involved
than simply building the static documentation.  

The full documentation makes heavy use of custom sphinx extensions to transform
recipes, notebooks, and inline code snippets into python scripts, IPython_
notebooks, or notebook cells that are executed when the docs are built.

To do this, we use IPython's nbconvert module to transform notebooks into
HTML. to simplify versioning of the notebook JSON format, we store notebooks in
an unevaluated state.  To generate evaluated notebooks, which could include
arbitrary output (text, images, HTML), we make use of runipy_, which provides
facilities to script notebook evaluation.

.. _runipy: https://github.com/paulgb/runipy
.. _IPython: http://ipython.org/

To build the full documentation, you will need yt, IPython, runipy, and all 
supplementary yt analysis modules installed. The following dependencies were 
used to generate the yt documentation during the release of yt 2.6 in late 2013.

- Sphinx_ 1.1.3
- IPython_ 1.1
- runipy_ (git hash f74458c2877)
- pandoc_ 1.11.1
- Rockstar halo finder 0.99.6
- SZpack_ 1.1.1
- ffmpeg_ 1.2.4 (compiled with libvpx support)
- JSAnimation_ (git hash 1b95cb3a3a)
- Astropy_ 0.2.5

.. _SZpack: http://www.cita.utoronto.ca/~jchluba/Science_Jens/SZpack/SZpack.html
.. _Astropy: http://astropy.org/
.. _Sphinx: http://sphinx-doc.org/
.. _pandoc: http://johnmacfarlane.net/pandoc/
.. _ffmpeg: http://www.ffmpeg.org/
.. _JSAnimation: https://github.com/jakevdp/JSAnimation

You will also need the full yt suite of `yt test data
<http://yt-project.org/data/>`_, including the larger datasets that are not used
in the answer tests.

You will need to ensure that your testing configuration is properly
configured and that all of the yt test data is in the testing directory.  See
:ref:`run_answer_testing` for more details on how to set up the testing
configuration.

Now that you have everything set up properly, go to the documentation directory
and build it using sphinx:

.. code-block:: bash

   cd $YT_HG/doc
   make html

If all of the dependencies are installed and all of the test data is in the
testing directory, this should churn away for a while (~ 1 hour) and 
eventually generate a docs build.  We suggest setting 
:code:`suppressStreamLogging = True` in your yt configuration (See 
:ref:`configuration-file`) to suppress large amounts of debug output from
yt.

To clean the docs build, use :code:`make clean`.  By default, :code:`make clean`
will not delete the autogenerated API docs, so use :code:`make fullclean` to
delete those as well.

Building the docs (hybrid)
--------------------------

It's also possible to create a custom sphinx build that builds a restricted set
of notebooks or scripts.  This can be accomplished by editing the Sphinx
:code:`conf.py` file included in the :code:`source` directory at the top level
of the docs.  The extensions included in the build are contained in the
:code:`extensions` list.  To disable an extension, simply remove it from the
list.  Doing so will raise a warning when sphinx encounters the directive in the
docs and will prevent sphinx from evaluating the directive.

As a concrete example, if one wanted to include the :code:`notebook`, and
:code:`notebook-cell` directives, but not the :code:`python-script` or
:code:`autosummary` directives, one would just need to comment out the lines
that append these extensions to the :code:`extensions` list. The resulting docs
build will be significantly quicker since it would avoid executing the lengthy
API autodocumentation as well as a large number of python script snippets in
the narrative docs.
