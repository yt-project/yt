.. _quickstart:

yt Quickstart
=============

The quickstart is a series of worked examples of how to use much of the
functionality of yt.  These are simple, short introductions to give you a taste
of what the code can do and are not meant to be detailed walkthroughs.

There are two ways in which you can go through the quickstart: interactively and
non-interactively.  We recommend the interactive method, but if you're pressed
on time, you can non-interactively go through the linked pages below and view the
worked examples.

To execute the quickstart interactively, you have a couple of options: 1) run
the notebook from your own system or 2) run it from the url
https://girder.hub.yt/#raft/5b5b4686323d12000122aa8a.
Option 1 requires an existing installation of yt (see
:ref:`installing-yt`), a copy of the yt source (which you may
already have depending on your installation choice), and a download of the
tutorial data-sets (total about 3 GB). If you know you are going to be a yt user
and have the time to download the data-sets, option 1 is a good choice. However,
if you're only interested in getting a feel for yt and its capabilities, or you
already have yt but don't want to spend time downloading the data, go ahead to
https://girder.hub.yt/#raft/5b5b4686323d12000122aa8a.

If you're running the tutorial from your own system and you do not already have
the yt repository, the easiest way to get the repository is to clone it using
git:

.. code-block:: bash

   git clone https://github.com/yt-project/yt

Now start the IPython notebook from within the repository (we presume you have
yt and [jupyter](https://jupyter.org/) installed):

.. code-block:: bash

   cd yt/doc/source/quickstart
   yt notebook

This command will give you information about the notebook server and how to
access it.  You will basically just pick a password (for security reasons) and then
redirect your web browser to point to the notebook server.
Once you have done so, choose "Introduction" from the list of
notebooks, which includes an introduction and information about how to download
the sample data.

.. warning:: The pre-filled out notebooks are *far* less fun than running them
             yourselves!  Check out the repo and give it a try.

Here are the notebooks, which have been filled in for inspection:

.. toctree::
   :maxdepth: 1

   introduction
   data_inspection
   simple_visualization
   data_objects_and_time_series
   derived_fields_and_profiles
   volume_rendering

.. note::

   The notebooks use sample datasets that are available for download at
   https://yt-project.org/data.  See :ref:`quickstart-introduction` for more
   details.

Let us know if you would like to contribute other example notebooks, or have
any suggestions for how these can be improved.
