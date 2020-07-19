.. _extensions:

Extension Packages
==================

.. note:: For some additional discussion, see `YTEP-0029
          <https://ytep.readthedocs.io/en/latest/YTEPs/YTEP-0029.html>`_, where
          this plan was designed.

As of version 3.3 of yt, we have put into place new methods for easing the
process of developing "extensions" to yt.  Extensions might be analysis
packages, visualization tools, or other software projects that use yt as a base
engine but that are versioned, developed and distributed separately.  This
brings with it the advantage of retaining control over the versioning,
contribution guidelines, scope, etc, while also providing a mechanism for
disseminating information about it, and potentially a method of interacting
with other extensions.

We have created a few pieces of infrastructure for developing extensions,
making them discoverable, and distributing them to collaborators.

If you have a module you would like to retain some external control over, or
that you don't feel would fit into yt, we encourage you to build it as an
extension module and distribute and version it independently.

Hooks for Extensions
--------------------

Starting with version 3.3 of yt, any package named with the prefix ``yt_`` is
importable from the namespace ``yt.extensions``.  For instance, the
``yt_interaction`` package ( https://bitbucket.org/data-exp-lab/yt_interaction
) is importable as ``yt.extensions.interaction``.

In subsequent versions, we plan to include in yt a catalog of known extensions
and where to find them; this will put discoverability directly into the code
base.

Extension Template
------------------

A template for starting an extension module (or converting an existing set of
code to an extension module) can be found at
https://github.com/yt-project/yt_extension_template .

To get started, download a zipfile of the template (
https://codeload.github.com/yt-project/yt_extension_template/zip/master ) and
follow the directions in ``README.md`` to modify the metadata.

Distributing Extensions
-----------------------

We encourage you to version on your choice of hosting platform (Bitbucket,
GitHub, etc), and to distribute your extension widely.  We are presently
working on deploying a method for listing extension modules on the yt webpage.
