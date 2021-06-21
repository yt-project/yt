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
functionality without opening a python interpreter.  The tools is a collection of
subcommands.  These can quickly making plots of slices and projections through a
dataset, updating yt's codebase, print basic statistics about a dataset, launch
an IPython notebook session, and more. To get a quick list of what is
available, just type:

.. code-block:: bash

   yt -h

This will print the list of available subcommands,

.. config_help:: yt

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
git repository.

First let's see what our options are for plotting:

.. code-block:: bash

  $ yt plot --help

There are many! We can choose whether we want a slice (default) or a
projection (``-p``), the field, the colormap, the center of the image, the
width and unit of width of the image, the limits, the weighting field for
projections, and on and on. By default the plotting command will execute the
same thing along all three axes, so keep that in mind if it takes three times
as long as you'd like! The center of a slice defaults to the center of
the domain, so let's just give that a shot and see what it looks like:

.. code-block:: bash

  $ yt plot DD0010/moving7_0010

Well, that looks pretty bad! What has happened here is that the center of the
domain only has some minor shifts in density, so the plot is essentially
incomprehensible. Let's try it again, but instead of slicing, let's project.
This is a line integral through the domain, and for the density field this
becomes a column density:

.. code-block:: bash

  $ yt plot -p DD0010/moving7_0010

Now that looks much better! Note that all three axes' projections appear
nearly indistinguishable, because of how the two spheres are located in the
domain. We could center our domain on one of the spheres and take a slice, as
well. Now let's see what the domain looks like with grids overlaid, using the
``--show-grids`` option:

.. code-block:: bash

  $ yt plot --show-grids -p DD0010/moving7_0010

We can now see all the grids in the field of view. If you want to
annotate your plot with a scale bar, you can use the
``--show-scale-bar`` option:

.. code-block:: bash

  $ yt plot --show-scale-bar -p DD0010/moving7_0010


Command-line subcommand summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

help
++++

Help lists all of the various command-line options in yt.

instinfo and version
++++++++++++++++++++

This gives information about where your yt installation is, what version
and changeset you're using and more.

load
++++

This will start the iyt interactive environment with your specified
dataset already loaded. See :ref:`interactive-prompt` for more details.

mapserver
+++++++++

Ever wanted to interact with your data using a
`google maps <http://maps.google.com/>`_-style interface?  Now you can by using the
yt mapserver. See :ref:`mapserver` for more details.

pastebin and pastebin_grab
++++++++++++++++++++++++++

The `pastebin <http://paste.yt-project.org/>`_ is an online location where
you can anonymously post code snippets and error messages to share with
other users in a quick, informal way. It is often useful for debugging
code or co-developing. By running the ``pastebin`` subcommand with a
text file, you send the contents of that file to an anonymous pastebin;

.. code-block:: bash

   yt pastebin my_script.py

By running the ``pastebin_grab`` subcommand with a pastebin number
(e.g. 1768), it will grab the contents of that pastebin
(e.g. the website http://paste.yt-project.org/show/1768 ) and send it to
STDOUT for local use.  See :ref:`pastebin` for more information.

.. code-block:: bash

   yt pastebin_grab 1768

upload
++++++

Upload a file to a public curldrop instance. Curldrop is a simple web
application that allows you to upload and download files straight from your
Terminal with an http client like e.g. curl. It was initially developed by
`Kevin Kennell <https://github.com/kennell/curldrop>`_ and later forked and
adjusted for yt’s needs. After a successful upload you will receive a url that
can be used to share the data with other people.

.. code-block:: bash

   yt upload my_file.tar.gz

plot
++++

This command generates one or many simple plots for a single dataset.
By specifying the axis, center, width, etc. (run ``yt help plot`` for
details), you can create slices and projections easily at the
command-line.

rpdb
++++

Connect to a currently running (on localhost) rpdb session. See
:ref:`remote-debugging` for more info.

notebook
++++++++

Launches a Jupyter notebook server and prints out instructions on how to open
an ssh tunnel to connect to the notebook server with a web browser. This is
most useful when you want to run a Jupyter notebook using CPUs on a remote
host.

stats
+++++

This subcommand provides you with some basic statistics on a given dataset.
It provides you with the number of grids and cells in each level, the time
of the dataset, and the resolution. It is tantamount to calling the
``Dataset.print_stats`` method.

Additionally, there is the option to print the minimum, maximum, or both for
a given field. The field is assumed to be density by default:

.. code-block:: bash

   yt stats GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150 --max --min

or a different field can be specified using the ``-f`` flag:

.. code-block:: bash

   yt stats GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150 --max --min -f gas,temperature

The field-related stats output from this command can be directed to a file using
the ``-o`` flag:

.. code-block:: bash

   yt stats GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150 --max -o out_stats.dat

update
++++++

This subcommand updates the yt installation to the most recent version for
your repository (e.g. stable, 2.0, development, etc.). Adding the ``--all``
flag will update the dependencies as well.

.. _upload-image:

upload_image
++++++++++++

Images are often worth a thousand words, so when you're trying to
share a piece of code that generates an image, or you're trying to
debug image-generation scripts, it can be useful to send your
co-authors a link to the image.  This subcommand makes such sharing
a breeze.  By specifying the image to share, ``upload_image`` automatically
uploads it anonymously to the website `imgur.com <https://imgur.com/>`_ and
provides you with a link to share with your collaborators.  Note that the
image *must* be in the PNG format in order to use this function.

delete_image
++++++++++++

The image uploaded using ``upload_image`` is assigned with a unique hash that
can be used to remove it. This subcommand provides an easy way to send a delete
request directly to the `imgur.com <https://imgur.com/>`_.

download
~~~~~~~~

This subcommand downloads a file from https://yt-project.org/data. Using ``yt download``,
one can download a file to:

* ``"test_data_dir"``: Save the file to the location specified in
  the ``"test_data_dir"`` configuration entry for test data.
* ``"supp_data_dir"``: Save the file to the location specified in
  the ``"supp_data_dir"`` configuration entry for supplemental data.
* Any valid path to a location on disk, e.g. ``/home/jzuhone/data``.

Examples:

.. code-block:: bash

   $ yt download apec_emissivity_v2.h5 supp_data_dir

.. code-block:: bash

   $ yt download GasSloshing.tar.gz test_data_dir

.. code-block:: bash

   $ yt download ZeldovichPancake.tar.gz /Users/jzuhone/workspace

If the configuration values ``"test_data_dir"`` or ``"supp_data_dir"`` have not
been set by the user, an error will be thrown.

Config helper
~~~~~~~~~~~~~

The :code:`yt config` command-line tool allows you to modify and access yt's
configuration without manually locating and opening the config file in an editor.
To get a quick list of available commands, just type:

.. code-block:: bash

   yt config -h

This will print the list of available subcommands:

.. config_help:: yt config


Since yt version 4, the configuration file is located in ``$XDG_CONFIG_HOME/yt/yt.toml`` adhering to the
`XDG Base Directory Specification
<https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html>`_.
Unless customized, this defaults to ``$HOME/.config/`` on Unix-like systems (macOS, Linux, ...).
The old configuration file (``$XDG_CONFIG_HOME/yt/ytrc``) is deprecated.
In order to perform an automatic migration of the old config, you are
encouraged to run:

.. code-block:: bash

   yt config migrate

This will convert your old config file to the toml format. The original file
will be moved to ``$XDG_CONFIG_HOME/yt/ytrc.bak``.

Examples
++++++++

Listing current content of the config file:

.. code-block:: bash

   $ yt config list
   [yt]
   log_level = 50

Obtaining a single config value by name:

.. code-block:: bash

   $ yt config get yt log_level
   50

Changing a single config value:

.. code-block:: bash

   $ yt config set yt log_level 10

Removing a single config entry:

.. code-block:: bash

   $ yt config rm yt log_level
