.. _bootcamp:

yt Bootcamp
===========

yt Bootcamp is a series of worked examples of how to use much of the 
funtionality of yt.  These are not meant to be detailed walkthroughs but simple,
short introductions to give you a taste of what the code can do.

There are two ways in which you can go through the bootcamp: interactively and 
non-interactively.  We recommend the interactive method, but if you're pressed 
on time, you can non-interactively go through the following pages and view the 
worked examples.

To execute the bootcamp interactively, you need to download the repository 
and start the IPython notebook.  The easiest way to get the repository is to 
use your already-installed mercurial program to grab it:

.. code-block:: bash

   hg clone https://bitbucket.org/yt_analysis/yt-doc

Now start the IPython notebook from within the repository:

.. code-block:: bash

   cd yt-doc/source/bootcamp
   yt notebook

This command will give you information about the notebook server and how to
access it.  You will basically just redirect your web browser to point to it.
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
   http://yt-project.org/data.  See :ref:`bootcamp-introduction` for more
   details.

Let us know if you would like to contribute other example notebooks, or have
any suggestions for how these can be improved.
