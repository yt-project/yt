How to Do a Release
-------------------

Periodically, the yt development community issues new releases. yt loosely follows
`semantic versioning <https://semver.org/>`_. The type of release can be read off
from the version number used. Version numbers should follow the scheme
``MAJOR.MINOR.PATCH``. There are three kinds of possible releases:

* Bugfix releases

  These releases should contain only fixes for bugs discovered in
  earlier releases and should not contain new features or API changes. Bugfix
  releases only increment the ``PATCH`` version number. Bugfix releases are
  generated from a dedicated backport branch that contains cherry-picked
  commits for bug fixes (handled semi-automatically as pull requests are
  merged), see :ref:`doing-a-bugfix-release` for more details. Version
  ``3.2.2`` is a bugfix release.

* Minor releases

  Minor releases include new features and fixes for bugs if the fix is
  determined to be too invasive for a bugfix release. Minor releases
  should *not* include backwards-incompatible changes and should not change APIs.  If an API change
  is deemed to be necessary, the old API should continue to function but might
  trigger deprecation warnings. Minor releases should happen by creating a new series
  branch off of ``main`` branch. Minor releases should increment the
  ``MINOR`` version number and reset the ``PATCH`` version number to zero.
  Version ``3.3.0`` is a minor release and has a corresponding series branch ``yt-3.3.x``.
  After release, the new series branch becomes the backport branch for any bug fix releases
  in the series.

* Major releases

  These releases happen when the development community decides to make major
  backwards-incompatible changes intentionally. In principle a major version release could
  include arbitrary changes to the library. Major version releases should only
  happen after extensive discussion and vetting among the developer and user
  community. Like minor releases, a major release should happen by creating
  a new series branch off of the ``main`` branch from which to release. Major releases should
  increment the ``MAJOR`` version number and reset the ``MINOR`` and ``PATCH``
  version numbers to zero. Version ``4.0.0`` is a major release.

The job of doing a release differs depending on the kind of release. Below, we
describe the necessary steps for each kind of release in detail. Several of the
following steps require that you have write privileges for the main yt GitHub
repository (if you're reading this, you probably already do).

.. note::

  This documentation assumes that your local copy of the yt repository has the
  main yt repository as the ``upstream`` remote (the default configuration when
  forking with ``gh fork yt-project/yt``). You can double check with:
  ``git remote -v``, which will display all the remote sources you have setup
  as well as their addresses. If you are missing an ``upstream`` remote (for
  example, after a standard ``git clone`` of your yt fork), you can add it
  with ``git remote add upstream git@github.com:yt-project/yt.git``.

.. _prepping-release-notes:

Prepping Release Notes
~~~~~~~~~~~~~~~~~~~~~~

Before starting the release process, it's useful to create a draft of a GitHub release. Doing this
at the start lets you share the upcoming release notes with other developers before the release
actually happens. To create a new draft release, go to
https://github.com/yt-project/yt/releases and click "Draft a new release". Use the version tag
as the release title (e.g., ``yt-4.4.1``). For the target tag, enter the tag that you will use
(see  :ref:`tagging-a-release-release`) and select the branch from which the release will be cut.

Now it's time to generate some release notes. To get a nice starting point, try using the ``uv`` script
`draft_yt_release_notes.py <https://gist.github.com/chrishavlin/248adea4296abb7bcdbaac952f304cf0>`_. You can
run the script directly from the gist using ``uv`` and the url to the raw script (check the readme at the gist).
This script will pull all the issues and pull requests that have been tagged to a specified GitHub milestone
and do some initial categorizing to produce a decent draft for release notes. You can then create a draft
release manually via the GitHub interface and copy in your draft notes, or if using ``gh``, you can do so
from the command line with ``gh release create --draft --notes-file <file> --target <branch>``. The initial
notes created by the ``draft_yt_release_notes.py`` will still need some manual attention: you should update
the release summary text, add any highlights if desired, add PRs that are missing. At present, you'll also
need to manually sort the frontend-specific changes by frontend to match previous release notes.

When updating the draft via the GitHub interface, click "save" at any point to save a draft.
Do NOT publish the draft yet.

.. _doing-a-bugfix-release:

Doing a Bugfix Release
~~~~~~~~~~~~~~~~~~~~~~

As described above, bugfix releases are regularly scheduled updates for minor
releases to ensure fixes for bugs make their way out to users in a timely
manner. Since bugfix releases should not include new features, we do not issue
bugfix releases by simply merging from the development ``main`` branch into
the designated backport branch for the series. Instead, commits are cherry-picked
from the ``main`` branch to a backport branch, and the backport branch is tagged
with the release tags.

Backport branches are named after the minor version they descend from, followed by
an ``x``. For instance, ``yt-4.0.x`` is the backport branch for all releases in the 4.0 series.
The backport branches are initially created during the release of a new minor or major
version (see :ref:`doing-a-minor-or-major-release`).

Backporting bugfixes can be done automatically using the `MeeseeksBox bot
<https://meeseeksbox.github.io>`_ with a GitHub milestone linked to a backport branch.
To set up a new milestone linked to a backport branch, click the New Milestone button
on the `GitHub interface <https://github.com/yt-project/yt/milestones>`_` then name it
with the version for the future release and a note within the description field:
``on-merge: backport to <name of backport branch>``, for example:
``on-merge: backport to yt-4.0.x``.

Then, every PR that was triaged into the milestone will be replicated as a
backport PR by the bot when it's merged into main. Some backports are non-trivial and
require human attention; if conflicts occur, the bot will provide detailed
instructions to perfom the task manually. If you forget to assign a backport branch label
before merging a PR, you can tag the bot in a comment on the merged PR to have it
create a new backport PR for the already merged PR (see `here <https://github.com/scientific-python/MeeseeksDev>`_
for a list of commands understood by the Meeseeks bot).

In short, a manual backport consist of 4 steps

- checking out the backport branch locally
- create a new branch from there
- cherry-picking the merge commit from the original PR with ``git cherry-pick -m1 <commit sha>``
- opening a PR to the backport branch

.. _doing-a-minor-or-major-release:

Doing a Minor or Major Release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is simpler than a bugfix release. First, make sure that all
deprecated features targeted for removal in the new release are removed from the
``main`` branch, ideally in a single PR. Such a PR can be issued at any point
between the previous minor or major release and the new one. Then, create a new
series branch off of the ``main`` branch (for example ``git switch -c yt-4.5.x``)
and push the new branch up to the yt repository.

.. code-block:: bash
  git fetch upstream
  git switch upstream/main
  git switch -c yt-4.5.x
  git push --set-upstream upstream yt-4.5.x

After the series branch is up, you will bump the version number and generate a git tag
as described below.

After the completion of the release, the new series branch becomes the
backport branch for subsequent bugfix releases.

Incrementing Version Numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before creating the tag for the release, you must increment the version numbers
that are hard-coded in a few files in the yt source so that version metadata
for the code is generated correctly. The paths relative to the root of the
repository for the three files that need to be edited are:

* ``doc/source/conf.py``

  The ``version`` and ``release`` variables need to be updated.

* ``pyproject.toml``

  The ``version`` variable needs to be updated

* ``yt/_version.py``

  The ``__version__`` variable must be updated.

To update these files, check out and update the branch that will be released (``main``
if doing a major or minor release, or the backport branch if doing a bugfix release).

Once these files have been updated, commit these updates and submit a pull request
targeting the branch that will be released. This is the commit we
will tag for the release.


.. _tagging-a-release:

Tagging a Release
~~~~~~~~~~~~~~~~~

After incrementing version numbers, checkout and update the branch that will be released
and actually create the tag by issuing the following command:

.. code-block:: bash

   git tag <tag-name>

Where ``<tag-name>`` follows the project's naming scheme for tags
(e.g. ``yt-3.2.1``). Once you are done, you will need to push the
tag to github::

  git push origin --tag

This assumes that you have configured the remote ``origin`` to point at the main
yt git repository. If you are doing a minor or major version number release, you
will also need to update back to the development branch and update the
development version numbers in the same files.


Publishing
~~~~~~~~~~

We distribute yt on two main channels: PyPI.org and conda-forge, in this order.

PyPI
++++

The publication process for PyPI is automated for the most part, via Github
actions, using ``.github/workflows/wheels.yaml``. Specifically, a release is
pushed to PyPI when a new git tag starting with ``yt-`` is pushed to the main
repo. Let's review the details here.

PyPI releases contain the source code (as a tarball), and wheels. Wheels are
compiled distributions of the source code. They are OS specific as well as
Python-version specific. Producing wheels for every supported combination of OS
and Python versions is done with `cibuildwheels
<https://cibuildwheel.readthedocs.org>`_

Upload to PyPI is automated via Github Actions `upload-artifact
<https://github.com/actions/upload-artifact>`_ and `download-artifact
<https://github.com/actions/upload-artifact>`_.

Note that automated uploads are currently perfomed using Matt Turk's
credentials.

If that worked, you can skip to the next section. Otherwise, upload can be
perfomed manually by first downloading the artifacts ``wheels`` and ``tarball``
from the workflow webpage, then at the command line (make sure that the
``dist`` directory doesn't exist or is empty)

.. code-block:: bash

   unzip tarball.zip -d dist
   unzip wheels.zip -d dist
   python -m pip install --upgrade twine
   twine upload dist/*

You will be prompted for your PyPI credentials and then the package should
upload. Note that for this to complete successfully, you will need an account on
PyPI and that account will need to be registered as an "owner" or "maintainer"
of the yt package.


``conda-forge``
+++++++++++++++

Conda-forge packages for yt are managed via the yt feedstock, located at
https://github.com/conda-forge/yt-feedstock. When a release is pushed to PyPI a
bot should detect a new version and issue a PR to the feedstock with the new
version automatically. When this feedstock is updated, make sure that the
SHA256 hash of the tarball matches the one you uploaded to PyPI and that
the version number matches the one that is being released.

In case the automated PR fails CI, feedstock maintainers are allowed to push to
the bot's branch with any fixes required.

Should you need to update the feedstock manually, you will
need to update the ``meta.yaml`` file located in the ``recipe`` folder in the
root of the feedstock repository. Most likely you will only need to update the
version number and the SHA256 hash of the tarball. If yt's dependencies change
you may also need to update the recipe. Once you have updated the recipe,
propose a pull request on github and merge it once all builds pass.


Announcing
~~~~~~~~~~

After the release is uploaded to `PyPI <https://pypi.org/project/yt/#files>`_ and
`conda-forge <https://anaconda.org/conda-forge/yt>`_, you should publish the
GitHub draft release (see :ref:`prepping-release-notes`) and then you should
send out an announcement e-mail to the yt mailing lists and the yt slack. Optionally,
for major and minor releases, you can send the announcement to other interested
mailing lists. Include a brief overview of the changes in the release and link to
the GitHub release.
