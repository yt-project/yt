.. _contributing-code:

How to Develop yt
=================

.. note:: If you already know how to use version control and are comfortable
   with handling it yourself, the quickest way to contribute to yt is to `fork
   us on BitBucket <http://hg.yt-project.org/yt/fork>`_, `make your changes
   <http://mercurial.selenic.com/>`_, and issue a `pull request
   <http://hg.yt-project.org/yt/pull>`_.  The rest of this document is just an
   explanation of how to do that.

yt is a community project!

We are very happy to accept patches, features, and bugfixes from any member of
the community!  yt is developed using mercurial, primarily because it enables
very easy and straightforward submission of changesets.  We're eager to hear
from you, and if you are developing yt, we encourage you to subscribe to the
`developer mailing list
<http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org>`_

Please feel free to hack around, commit changes, and send them upstream.  If
you're new to Mercurial, these three resources are pretty great for learning
the ins and outs:

* http://hginit.com/
* http://hgbook.red-bean.com/read/
* http://mercurial.selenic.com/

The commands that are essential for using mercurial include:

* ``hg commit`` which commits changes in the working directory to the
  repository, creating a new "changeset object."
* ``hg add`` which adds a new file to be tracked by mercurial.  This does
  not change the working directory.
* ``hg pull`` which pulls (from an optional path specifier) changeset
  objects from a remote source.  The working directory is not modified.
* ``hg push`` which sends (to an optional path specifier) changeset objects
  to a remote source.  The working directory is not modified.
* ``hg log`` which shows a log of all changeset objects in the current
  repository.  Use ``-g`` to show a graph of changeset objects and their
  relationship.
* ``hg update`` which (with an optional "revision" specifier) updates the
  state of the working directory to match a changeset object in the
  repository.
* ``hg merge`` which combines two changesets to make a union of their lines
  of development.  This updates the working directory.

Keep in touch, and happy hacking!  We also provide `doc/coding_styleguide.txt`
and an example of a fiducial docstring in `doc/docstring_example.txt`.  Please
read them before hacking on the codebase, and feel free to email any of the
mailing lists for help with the codebase.

.. _bootstrap-dev:

Submitting Changes
------------------

We provide a brief introduction to submitting changes here.  yt thrives on the
strength of its communities ( http://arxiv.org/abs/1301.7064 has further
discussion) and we encourage contributions from any user.  While we do not
discuss in detail version control, mercurial or the advanced usage of
BitBucket, we do provide an outline of how to submit changes and we are happy
to provide further assistance or guidance.

Licensing
+++++++++

yt has, with the 2.6 release, been `relicensed
<http://blog.yt-project.org/post/Relicensing.html>`_ under the BSD 3-clause
license.  Previously versions were released under the GPLv3.

All contributed code must be BSD-compatible.  If you'd rather not license in
this manner, but still want to contribute, please consider creating an external
package, which we'll happily link to.

.. _requirements-for-code-submission:

Requirements for Code Submission
++++++++++++++++++++++++++++++++

Modifications to the code typically fall into one of three categories, each of
which have different requirements for acceptance into the code base.  These
requirements are in place for a few reasons -- to make sure that the code is
maintainable, testable, and that we can easily include information about
changes in changelogs during the release procedure.  (See `YTEP-0008
<https://ytep.readthedocs.org/en/latest/YTEPs/YTEP-0008.html>`_ for more
detail.)

* New Features

  * New unit tests (possibly new answer tests) (See :ref:`testing`)
  * Docstrings in the source code for the public API
  * Addition of new feature to the narrative documentation (See :ref:`writing_documentation`)
  * Addition of cookbook recipe (See :ref:`writing_documentation`) 
  * Issue created on issue tracker, to ensure this is added to the changelog

* Extension or Breakage of API in Existing Features

  * Update existing narrative docs and docstrings (See :ref:`writing_documentation`) 
  * Update existing cookbook recipes (See :ref:`writing_documentation`) 
  * Modify of create new unit tests (See :ref:`testing`)
  * Issue created on issue tracker, to ensure this is added to the changelog

* Bug fixes

  * Unit test is encouraged, to ensure breakage does not happen again in the
    future. (See :ref:`testing`)
  * Issue created on issue tracker, to ensure this is added to the changelog

When submitting, you will be asked to make sure that your changes meet all of
these requirements.  They are pretty easy to meet, and we're also happy to help
out with them.  In :ref:`code-style-guide` there is a list of handy tips for
how to structure and write your code.

.. _mercurial-with-yt:

How to Use Mercurial with yt
++++++++++++++++++++++++++++

This document doesn't cover detailed mercurial use, but on IRC we are happy to
walk you through any troubles you might have.  Here are some suggestions
for using mercurial with yt:

* Named branches are to be avoided.  Try using bookmarks (``hg bookmark``) to
  track work.  (`More <http://mercurial.selenic.com/wiki/Bookmarks>`_)
* Make sure you set a username in your ``~/.hgrc`` before you commit any
  changes!  All of the tutorials above will describe how to do this as one of
  the very first steps.
* When contributing changes, you might be asked to make a handful of
  modifications to your source code.  We'll work through how to do this with
  you, and try to make it as painless as possible.
* Please avoid deleting your yt forks, as that eliminates the code review
  process from BitBucket's website.
* In all likelihood, you only need one fork.  To keep it in sync, you can
  sync from the website.  (See Bitbucket's `Blog Post
  <http://blog.bitbucket.org/2013/02/04/syncing-and-merging-come-to-bitbucket/>`_
  about this.)
* If you run into any troubles, stop by IRC (see :ref:`irc`) or the mailing
  list.

.. _building-yt:

Building yt
+++++++++++

If you have made changes to any C or Cython (``.pyx``) modules, you have to
rebuild yt.  If your changes have exclusively been to Python modules, you will
not need to re-build, but (see below) you may need to re-install.  

If you are running from a clone that is executable in-place (i.e., has been
installed via the installation script or you have run ``setup.py develop``) you
can rebuild these modules by executing:

.. code-block:: bash

  $ python2.7 setup.py develop

If you have previously "installed" via ``setup.py install`` you have to
re-install:

.. code-block:: bash

  $ python2.7 setup.py install

Only one of these two options is needed.

.. _windows-developing:

Developing yt on Windows
^^^^^^^^^^^^^^^^^^^^^^^^

If you plan to develop yt on Windows, it is necessary to use the `MinGW
<http://www.mingw.org/>`_ gcc compiler that can be installed using the `Anaconda
Python Distribution <https://store.continuum.io/cshop/anaconda/>`_. The libpython package must be
installed from Anaconda as well. These can both be installed with a single command:

.. code-block:: bash

  $ conda install libpython mingw

Additionally, the syntax for the setup command is slightly different; you must type:

.. code-block:: bash

  $ python2.7 setup.py build --compiler=mingw32 develop

or

.. code-block:: bash

  $ python2.7 setup.py build --compiler=mingw32 install

.. _sharing-changes:

Making and Sharing Changes
++++++++++++++++++++++++++

The simplest way to submit changes to yt is to do the following:

* Build yt from the mercurial repository
* Navigate to the root of the yt repository 
* Make some changes and commit them
* Fork the `yt repository on BitBucket <https://bitbucket.org/yt_analysis/yt>`_
* Push the changesets to your fork
* Issue a pull request.

Here's a more detailed flowchart of how to submit changes.

#. If you have used the installation script, the source code for yt can be
   found in ``$YT_DEST/src/yt-hg``.  Alternatively see
   :ref:`source-installation` for instructions on how to build yt from the
   mercurial repository. (Below, in :ref:`reading-source`, we describe how to
   find items of interest.)  
#. Edit the source file you are interested in and
   test your changes.  (See :ref:`testing` for more information.)
#. Fork yt on BitBucket.  (This step only has to be done once.)  You can do
   this at: https://bitbucket.org/yt_analysis/yt/fork.  Call this repository
   yt.
#. Commit these changes, using ``hg commit``.  This can take an argument
   which is a series of filenames, if you have some changes you do not want
   to commit.
#. If your changes include new functionality or cover an untested area of the
   code, add a test.  (See :ref:`testing` for more information.)  Commit
   these changes as well.
#. Push your changes to your new fork using the command::

      hg push -r . https://bitbucket.org/YourUsername/yt/
 
   If you end up doing considerable development, you can set an alias in the
   file ``.hg/hgrc`` to point to this path.

   .. note::
     Note that the above approach uses HTTPS as the transfer protocol
     between your machine and BitBucket.  If you prefer to use SSH - or
     perhaps you're behind a proxy that doesn't play well with SSL via
     HTTPS - you may want to set up an `SSH key`_ on BitBucket.  Then, you use
     the syntax ``ssh://hg@bitbucket.org/YourUsername/yt``, or equivalent, in
     place of ``https://bitbucket.org/YourUsername/yt`` in Mercurial commands.
     For consistency, all commands we list in this document will use the HTTPS
     protocol.

     .. _SSH key: https://confluence.atlassian.com/display/BITBUCKET/Set+up+SSH+for+Mercurial

#. Issue a pull request at
   https://bitbucket.org/YourUsername/yt/pull-request/new

During the course of your pull request you may be asked to make changes.  These
changes may be related to style issues, correctness issues, or even requesting
tests.  The process for responding to pull request code review is relatively
straightforward.

#. Make requested changes, or leave a comment indicating why you don't think
   they should be made.
#. Commit those changes to your local repository.
#. Push the changes to your fork::

      hg push https://bitbucket.org/YourUsername/yt/

#. Your pull request will be automatically updated.

.. _multiple-PRs:

Working with Multiple BitBucket Pull Requests
+++++++++++++++++++++++++++++++++++++++++++++

Once you become active developing for yt, you may be working on
various aspects of the code or bugfixes at the same time.  Currently,
BitBucket's *modus operandi* for pull requests automatically updates
your active pull request with every ``hg push`` of commits that are a
descendant of the head of your pull request.  In a normal workflow,
this means that if you have an active pull request, make some changes
locally for, say, an unrelated bugfix, then push those changes back to
your fork in the hopes of creating a *new* pull request, you'll
actually end up updating your current pull request!

There are a few ways around this feature of BitBucket that will allow
for multiple pull requests to coexist; we outline one such method
below.  We assume that you have a fork of yt at
``http://bitbucket.org/YourUsername/Your_yt`` (see
:ref:`sharing-changes` for instructions on creating a fork) and that
you have an active pull request to the main repository.

The main issue with starting another pull request is to make sure that
your push to BitBucket doesn't go to the same head as your
existing pull request and trigger BitBucket's auto-update feature.
Here's how to get your local repository away from your current pull
request head using `revsets <http://www.selenic.com/hg/help/revsets>`_
and your ``hgrc`` file:
   
#. Set up a Mercurial path for the main yt repository (note this is a convenience
   step and only needs to be done once).  Add the following to your
   ``Your_yt/.hg/hgrc``::

     [paths]
     upstream = https://bitbucket.org/yt_analysis/yt

   This will create a path called ``upstream`` that is aliased to the URL of the
   main yt repository.
#. Now we'll use revsets_ to update your local repository to the tip of the
   ``upstream`` path:

   .. code-block:: bash

      $ hg pull upstream
      $ hg update -r "remote(tip,'upstream')"

After the above steps, your local repository should be at the tip of
the main yt repository.  If you find yourself doing this a lot, it may
be worth aliasing this task in your ``hgrc`` file by adding something like::

  [alias]
  myupdate = update -r "remote(tip,'upstream')"

And then you can just issue ``hg myupdate`` to get at the tip of the main yt repository.

Make sure you are on the branch you want to be on, and then you can
make changes and ``hg commit`` them.  If you prefer working with
`bookmarks <http://mercurial.selenic.com/wiki/Bookmarks>`_, you may
want to make a bookmark before committing your changes, such as ``hg
bookmark mybookmark``.

To push to your fork on BitBucket if you didn't use a bookmark, you issue the following:

.. code-block:: bash

  $ hg push -r . -f https://bitbucket.org/YourUsername/Your_yt

The ``-r .`` means "push only the commit I'm standing on and any ancestors."  The
``-f`` is to force Mecurial to do the push since we are creating a new remote head.

Note that if you *did* use a bookmark, you don't have to force the push, but you do
need to push the bookmark; in other words do the following instead of the above:

.. code-block:: bash
		
   $ hg push -B mybookmark https://bitbucket.org/YourUsername/Your_yt

The ``-B`` means "publish my bookmark and any relevant changesets to the remote server."
		
You can then go to the BitBucket interface and issue a new pull request based on
your last changes, as usual.

How To Get The Source Code For Editing
--------------------------------------

yt is hosted on BitBucket, and you can see all of the yt repositories at
http://hg.yt-project.org/.  With the yt installation script you should have a
copy of Mercurial for checking out pieces of code.  Make sure you have followed
the steps above for bootstrapping your development (to assure you have a
bitbucket account, etc.)

In order to modify the source code for yt, we ask that you make a "fork" of the
main yt repository on bitbucket.  A fork is simply an exact copy of the main
repository (along with its history) that you will now own and can make
modifications as you please.  You can create a personal fork by visiting the yt
bitbucket webpage at https://bitbucket.org/yt_analysis/yt/ .  After logging in,
you should see an option near the top right labeled "fork".  Click this option,
and then click the fork repository button on the subsequent page.  You now have
a forked copy of the yt repository for your own personal modification.

This forked copy exists on the bitbucket repository, so in order to access
it locally, follow the instructions at the top of that webpage for that
forked repository, namely run at a local command line:

.. code-block:: bash

   $ hg clone http://bitbucket.org/<USER>/<REPOSITORY_NAME>

This downloads that new forked repository to your local machine, so that you
can access it, read it, make modifications, etc.  It will put the repository in
a local directory of the same name as the repository in the current working
directory.  You can see any past state of the code by using the hg log command.
For example, the following command would show you the last 5 changesets
(modifications to the code) that were submitted to that repository.

.. code-block:: bash

   $ cd <REPOSITORY_NAME>
   $ hg log -l 5

Using the revision specifier (the number or hash identifier next to each
changeset), you can update the local repository to any past state of the
code (a previous changeset or version) by executing the command:

.. code-block:: bash

   $ hg up revision_specifier

Lastly, if you want to use this new downloaded version of your yt repository as
the *active* version of yt on your computer (i.e. the one which is executed when
you run yt from the command line or the one that is loaded when you do ``import
yt``), then you must "activate" it using the following commands from within the
repository directory.

.. code-block:: bash

   $ cd <REPOSITORY_NAME>
   $ python2.7 setup.py develop

This will rebuild all C modules as well.

.. _reading-source:

How To Read The Source Code
---------------------------

If you just want to *look* at the source code, you may already have it on your
computer.  If you build yt using the install script, the source is available at
``$YT_DEST/src/yt-hg``.  See :ref:`source-installation` for more details about
to obtain the yt source code if you did not build yt using the install
script. 

The root directory of the yt mercurial repository contains a number of
subdirectories with different components of the code.  Most of the yt source
code is contained in the yt subdirectory.  This directory its self contains
the following subdirectories:

``frontends``
   This is where interfaces to codes are created.  Within each subdirectory of
   yt/frontends/ there must exist the following files, even if empty:

   * ``data_structures.py``, where subclasses of AMRGridPatch, Dataset
     and AMRHierarchy are defined.
   * ``io.py``, where a subclass of IOHandler is defined.
   * ``fields.py``, where fields we expect to find in datasets are defined
   * ``misc.py``, where any miscellaneous functions or classes are defined.
   * ``definitions.py``, where any definitions specific to the frontend are
     defined.  (i.e., header formats, etc.)

``fields``
   This is where all of the derived fields that ship with yt are defined.

``geometry`` 
   This is where geometric helpler routines are defined. Handlers
   for grid and oct data, as well as helpers for coordinate transformations
   can be found here.

``visualization``
   This is where all visualization modules are stored.  This includes plot
   collections, the volume rendering interface, and pixelization frontends.

``data_objects``
   All objects that handle data, processed or unprocessed, not explicitly
   defined as visualization are located in here.  This includes the base
   classes for data regions, covering grids, time series, and so on.  This
   also includes derived fields and derived quantities.

``analysis_modules``
   This is where all mechanisms for processing data live.  This includes
   things like clump finding, halo profiling, halo finding, and so on.  This
   is something of a catchall, but it serves as a level of greater
   abstraction that simply data selection and modification.

``gui``
   This is where all GUI components go.  Typically this will be some small
   tool used for one or two things, which contains a launching mechanism on
   the command line.

``utilities``
   All broadly useful code that doesn't clearly fit in one of the other
   categories goes here.

``extern`` 
   Bundled external modules (i.e. code that was not written by one of
   the yt authors but that yt depends on) lives here.


If you're looking for a specific file or function in the yt source code, use
the unix find command:

.. code-block:: bash

   $ find <DIRECTORY_TREE_TO_SEARCH> -name '<FILENAME>'

The above command will find the FILENAME in any subdirectory in the
DIRECTORY_TREE_TO_SEARCH.  Alternatively, if you're looking for a function
call or a keyword in an unknown file in a directory tree, try:

.. code-block:: bash

   $ grep -R <KEYWORD_TO_FIND> <DIRECTORY_TREE_TO_SEARCH>

This can be very useful for tracking down functions in the yt source.

.. _code-style-guide:

Code Style Guide
----------------

To keep things tidy, we try to stick with a couple simple guidelines.

General Guidelines
++++++++++++++++++

* In general, follow `PEP-8 <http://www.python.org/dev/peps/pep-0008/>`_ guidelines.
* Classes are ConjoinedCapitals, methods and functions are
  ``lowercase_with_underscores.``
* Use 4 spaces, not tabs, to represent indentation.
* Line widths should not be more than 80 characters.
* Do not use nested classes unless you have a very good reason to, such as
  requiring a namespace or class-definition modification.  Classes should live
  at the top level.  ``__metaclass__`` is exempt from this.
* Do not use unnecessary parentheses in conditionals.  ``if((something) and
  (something_else))`` should be rewritten as ``if something and
  something_else``.  Python is more forgiving than C.
* Avoid copying memory when possible. For example, don't do ``a =
  a.reshape(3,4)`` when ``a.shape = (3,4)`` will do, and ``a = a * 3`` should be
  ``np.multiply(a, 3, a)``.
* In general, avoid all double-underscore method names: ``__something`` is
  usually unnecessary.
* Doc strings should describe input, output, behavior, and any state changes
  that occur on an object.  See the file `doc/docstring_example.txt` for a
  fiducial example of a docstring.

API Guide
+++++++++

* Do not import "*" from anything other than ``yt.funcs``.
* Internally, only import from source files directly; instead of: ``from
  yt.visualization.api import SlicePlot`` do
  ``from yt.visualization.plot_window import SlicePlot``.
* Numpy is to be imported as ``np``.
* Do not use too many keyword arguments.  If you have a lot of keyword
  arguments, then you are doing too much in ``__init__`` and not enough via
  parameter setting.
* In function arguments, place spaces before commas.  ``def something(a,b,c)``
  should be ``def something(a, b, c)``.
* Don't create a new class to replicate the functionality of an old class --
  replace the old class.  Too many options makes for a confusing user
  experience.
* Parameter files external to yt are a last resort.
* The usage of the ``**kwargs`` construction should be avoided.  If they
  cannot be avoided, they must be explained, even if they are only to be
  passed on to a nested function.
* Constructor APIs should be kept as *simple* as possible.
* Variable names should be short but descriptive.
* No global variables!

Variable Names and Enzo-isms
++++++++++++++++++++++++++++

* Avoid Enzo-isms.  This includes but is not limited to:

  + Hard-coding parameter names that are the same as those in Enzo.  The
    following translation table should be of some help.  Note that the
    parameters are now properties on a Dataset subclass: you access them
    like ``ds.refine_by`` .

    - ``RefineBy `` => `` refine_by``
    - ``TopGridRank `` => `` dimensionality``
    - ``TopGridDimensions `` => `` domain_dimensions``
    - ``InitialTime `` => `` current_time``
    - ``DomainLeftEdge `` => `` domain_left_edge``
    - ``DomainRightEdge `` => `` domain_right_edge``
    - ``CurrentTimeIdentifier `` => `` unique_identifier``
    - ``CosmologyCurrentRedshift `` => `` current_redshift``
    - ``ComovingCoordinates `` => `` cosmological_simulation``
    - ``CosmologyOmegaMatterNow `` => `` omega_matter``
    - ``CosmologyOmegaLambdaNow `` => `` omega_lambda``
    - ``CosmologyHubbleConstantNow `` => `` hubble_constant``

  + Do not assume that the domain runs from 0 to 1.  This is not true
    for many codes and datasets.
