.. _command-line:

Command-Line Usage
------------------

.. _interactive-prompt:

Interactive Prompt
~~~~~~~~~~~~~~~~~~

The interactive prompt offers a number of excellent opportunities for
exploration of data.  While there are challenges for repeatability, and some
operations will be more challenging to operate in parallel, interactive prompts
can be exceptionally useful for debugging, exploring, and tweaking plots.

You can start up an empty shell, with a handful of useful yt utilities (such as
tab-completion and pre-imported modules) by executing:

.. code-block:: bash

   iyt

The other option, which is shorthand for "iyt plus dataset loading" is to use
the command-line tool (see :ref:`command-line`) with the ``load`` subcommand
and to specify a dataset.  For instance:

.. code-block:: bash

   yt load cosmoSim_coolhdf5_chk_0026

or

.. code-block:: bash

   yt load DD0030/DD0030

This will spawn ``iyt``, but the dataset given on the command line will
already be in the namespace as ``ds``.  With interactive mode, you can use the
``pylab`` module to interactively plot.

Command-line Functions
~~~~~~~~~~~~~~~~~~~~~~

The :code:`yt` command-line tool allows you to access some of yt's basic
funcionality without opening a python interpreter.  The tools is a collection of
subcommands.  These can quickly making plots of slices and projections through a
dataset, updating yt's codebase, print basic statistics about a dataset, laucnh
an IPython notebook session, and more.  To get a quick list of what is
available, just type:

.. code-block:: bash

   yt -h

This will print the list of available subcommands,

.. code-block:: bash

    help                Print help message
    bootstrap_dev       Bootstrap a yt development environment
    bugreport           Report a bug in yt
    hub_register        Register a user on the Hub: http://hub.yt-project.org/
    hub_submit          Submit a mercurial repository to the yt Hub
                        (http://hub.yt-project.org/), creating a BitBucket
                        repo in the process if necessary.
    instinfo            Get some information about the yt installation
    version             Get some information about the yt installation (this
                        is an alias for instinfo).
    load                Load a single dataset into an IPython instance
    mapserver           Serve a plot in a GMaps-style interface
    pastebin            Post a script to an anonymous pastebin
    pastebin_grab       Print an online pastebin to STDOUT for local use.
    upload_notebook     Upload an IPython notebook to hub.yt-project.org.
    plot                Create a set of images
    rpdb                Connect to a currently running (on localhost) rpd
                        session. Commands run with --rpdb will trigger an rpdb
                        session with any uncaught exceptions.
    notebook            Run the IPython Notebook
    stats               Print stats and max/min value of a given field (if
                        requested), for one or more datasets (default field is
                        Density)
    update              Update the yt installation to the most recent version
    upload_image        Upload an image to imgur.com. Must be PNG.


To execute any such function, simply run:

.. code-block:: bash

   yt <subcommand>

Finally, to identify the options associated with any of these subcommand, run:

.. code-block:: bash

   yt <subcommand> -h

Plotting from the command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we'll discuss plotting from the command line, then we will give a brief
summary of the functionality provided by each command line subcommand. This
example uses the :code:`DD0010/moving7_0010` dataset distributed in the yt
mercurial repository.

First let's see what our options are for plotting:

.. code-block:: bash

  $ yt plot --help

There are many!  We can choose whether we want a slice (default) or a
projection (``-p``), the field, the colormap, the center of the image, the
width and unit of width of the image, the limits, the weighting field for
projections, and on and on.  By default the plotting command will execute the
same thing along all three axes, so keep that in mind if it takes three times
as long as you'd like!  The center of a slice defaults to the center of
the domain, so let's just give that a shot and see what it looks like:

.. code-block:: bash

  $ yt plot DD0010/moving7_0010

Well, that looks pretty bad!  What has happened here is that the center of the
domain only has some minor shifts in density, so the plot is essentially
incomprehensible.  Let's try it again, but instead of slicing, let's project.
This is a line integral through the domain, and for the density field this
becomes a column density.:

.. code-block:: bash

  $ yt plot -p DD0010/moving7_0010

Now that looks much better!  Note that all three axes' projections appear
nearly indistinguishable, because of how the two spheres are located in the
domain.  We could center our domain on one of the spheres and take a slice, as
well.  Now let's see what the domain looks like with grids overlaid, using the
``--show-grids`` option.:

.. code-block:: bash

  $ yt plot --show-grids -p DD0010/moving7_0010

We can now see all the grids in the field of view.

Command-line subcommand summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

help
++++

Help lists all of the various command-line options in yt.


bugreport         
+++++++++

Encountering a bug in your own code can be a big hassle, but it can be 
exponentially worse to find it in someone else's.  That's why we tried to 
make it as easy as possible for users to report bugs they find in yt.  
After you go through the necessary channels to make sure you're not just
making a mistake (see :ref:`asking-for-help`), you can submit bug 
reports using this nice utility.

instinfo and version
++++++++++++++++++++

This gives information about where your yt installation is, what version
and changeset you're using and more. 

load
++++

This will start the iyt interactive environment with your specified 
dataset already loaded.  See :ref:`interactive-prompt` for more details.

mapserver
+++++++++

Ever wanted to interact with your data using the 
`google maps <http://maps.google.com/>`_ interface?  Now you can by using the
yt mapserver.  See :ref:`mapserver` for more details.

pastebin and pastebin_grab
++++++++++++++++++++++++++

The `pastebin <http://paste.yt-project.org/>`_ is an online location where 
you can anonymously post code snippets and error messages to share with 
other users in a quick, informal way.  It is often useful for debugging 
code or co-developing.  By running the ``pastebin`` subcommand with a 
text file, you send the contents of that file to an anonymous pastebin; 

.. code-block:: bash

   yt pastebin my_script.py

By running the ``pastebin_grab`` subcommand with a pastebin number 
(e.g. 1768), it will grab the contents of that pastebin 
(e.g. the website http://paste.yt-project.org/show/1768 ) and send it to 
STDOUT for local use.  See :ref:`pastebin` for more information.

.. code-block:: bash

   yt pastebin_grab 1768

plot
++++

This command generates one or many simple plots for a single dataset.  
By specifying the axis, center, width, etc. (run ``yt help plot`` for 
details), you can create slices and projections easily at the 
command-line.

upload_notebook
+++++++++++++++

This command will accept the filename of a ``.ipynb`` file (generated from an
IPython notebook session) and upload it to the `yt hub
<http://hub.yt-project.org/>` where others will be able to view it, and
download it.  This is an easy method for recording a sequence of commands,
their output, narrative information, and then sharing that with others.  These
notebooks will be viewable online, and the appropriate URLs will be returned on
the command line.


rpdb
++++

Connect to a currently running (on localhost) rpd session.

notebook
++++++++

Launches an IPython notebook server and prints out instructions on how to open
an ssh tunnel to connect to the notebook server with a web browser.  This is
most useful when you want to run an IPython notebook using CPUs on a remote
host.

stats
+++++

This subcommand provides you with some basic statistics on a given dataset.
It provides you with the number of grids and cells in each level, the time
of the dataset, the resolution, and the maximum density in a variety of units.
It is tantamount to performing the ``print_stats()`` inside of yt.

update
++++++

This subcommand updates the yt installation to the most recent version for
your repository (e.g. stable, 2.0, development, etc.).  Adding the ``--all`` 
flag will update the dependencies as well.

.. _upload-image:

upload_image
++++++++++++

Images are often worth a thousand words, so when you're trying to 
share a piece of code that generates an image, or you're trying to 
debug image-generation scripts, it can be useful to send your
co-authors a link to the image.  This subcommand makes such sharing 
a breeze.  By specifying the image to share, ``upload_image`` automatically
uploads it anonymously to the website `imgur.com <http://imgur.com/>`_ and
provides you with a link to share with your collaborators.  Note that the
image *must* be in the PNG format in order to use this function.
