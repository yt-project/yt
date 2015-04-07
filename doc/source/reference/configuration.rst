Customizing yt: The Configuration and Plugin Files
==================================================

yt features ways to customize it to your personal preferences in terms of 
how much output it displays, loading custom fields, loading custom colormaps, 
accessing test datasets regardless of where you are in the file system, etc.
This customization is done through :ref:`configuration-file` and 
:ref:`plugin-file` both of which exist in your ``$HOME/.yt`` directory.

.. _configuration-file:

The Configuration File
----------------------

The configuration is a simple text file setting internal yt variables to
custom default values to be used in future sessions.

Configuration File Format
^^^^^^^^^^^^^^^^^^^^^^^^^

yt will look for and recognize the file ``$HOME/.yt/config`` as a configuration
file, containing several options that can be modified and adjusted to control
runtime behavior.  For example, a sample ``$HOME/.yt/config`` file could look
like:

.. code-block:: none
    
   [yt]
   loglevel = 1
   maximumstoreddatasets = 10000

This configuration file would set the logging threshold much lower, enabling
much more voluminous output from yt.  Additionally, it increases the number of
datasets tracked between instantiations of yt.

Configuration Options At Runtime
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to setting parameters in the configuration file itself, you can set
them at runtime.  

.. warning:: Several parameters are only accessed when yt starts up: therefore,
   if you want to modify any configuration parameters at runtime, you should
   execute the appropriate commands at the *very top* of your script!

This involves importing the configuration object and then setting a given
parameter to be equal to a specific string.  Note that even for items that
accept integers, floating points and other non-string types, you *must* set
them to be a string or else the configuration object will consider them broken.

Here is an example script, where we adjust the logging at startup:

.. code-block:: python

   import yt
   yt.funcs.mylog.setLevel(1)

   ds = yt.load("my_data0001")
   ds.print_stats()

This has the same effect as setting ``loglevel = 1`` in the configuration
file. Note that a log level of 1 means that all log messages are printed to
stdout.  To disable logging, set the log level to 50.

Setting Configuration On the Command Line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Options can also be set directly on the command line by specifying a
command-line option.  For instance, if you are running the script
``my_script.py`` you can specify a configuration option with the ``--config``
argument.  As an example, to lower the log level (thus making it more verbose)
you can specify:

.. code-block:: bash

   $ python2.7 my_script.py --config loglevel=1

Any configuration option specific to yt can be specified in this manner.  One
common configuration option would be to disable serialization:

.. code-block:: bash

   $ python2.7 my_script.py --config serialize=False

This way projections are always re-created.

Available Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following external parameters are available.  A number of parameters are
used internally.

* ``coloredlogs`` (default: ``'False'``): Should logs be colored?
* ``loadfieldplugins`` (default: ``'True'``): Do we want to load the plugin file?
* ``pluginfilename``  (default ``'my_plugins.py'``) The name of our plugin file.
* ``logfile`` (default: ``'False'``): Should we output to a log file in the
  filesystem?
* ``loglevel`` (default: ``'20'``): What is the threshold (0 to 50) for
  outputting log files?
* ``test_data_dir`` (default: ``'/does/not/exist'``): The default path the
  ``load()`` function searches for datasets when it cannot find a dataset in the
  current directory.
* ``notebook_password`` (default: empty): If set, this will be fed to the
  IPython notebook created by ``yt notebook``.  Note that this should be an
  sha512 hash, not a plaintext password.  Starting ``yt notebook`` with no
  setting will provide instructions for setting this.
* ``serialize`` (default: ``'True'``): Are we allowed to write to the ``.yt`` file?
* ``sketchfab_api_key`` (default: empty): API key for http://sketchfab.com/ for
  uploading AMRSurface objects.
* ``suppressStreamLogging`` (default: ``'False'``): If true, execution mode will be
  quiet.
* ``stdoutStreamLogging`` (default: ``'False'``): If true, logging is directed
  to stdout rather than stderr
* ``skip_dataset_cache`` (default: ``'False'``): If true, automatic caching of datasets
  is turned off.

.. _plugin-file:

The Plugin File
---------------

The plugin file is a means of creating custom fields, quantities, data 
objects, colormaps, and other code classes and objects to be used in future
yt sessions without modifying the source code directly.  


.. note::

   The ``my_plugins.py`` is only parsed inside of ``yt.mods``, so in order
   to use it, you must load yt with either: ``import yt.mods as yt``
   or ``from yt.mods import *``.  You can tell that your
   plugins file is being parsed by watching for a logging message when you
   import yt.  Note that both the ``yt load`` and ``iyt`` command line entry
   points invoke ``from yt.mods import *``, so the ``my_plugins.py`` file
   will be parsed if you enter yt that way.

Plugin File Format
^^^^^^^^^^^^^^^^^^

yt will look for and recognize the file ``$HOME/.yt/my_plugins`` as a plugin
file, which should contain python code.  If accessing yt functions and classes
they will not require the ``yt.`` prefix, because of how they are loaded.
It is executed at the bottom of ``yt.mods``, and so
it is provided with the entire namespace available in the module ``yt.mods``.

For example, if I created a plugin file containing:

.. code-block:: python

   def _myfunc(field, data):
       return np.random.random(data["density"].shape)
   add_field("random", function=_myfunc, units='auto')

then all of my data objects would have access to the field ``some_quantity``.

You can also define other convenience functions in your plugin file.  For
instance, you could define some variables or functions, and even import common
modules:

.. code-block:: python

   import os

   HOMEDIR="/home/username/"
   RUNDIR="/scratch/runs/"

   def load_run(fn):
       if not os.path.exists(RUNDIR + fn):
           return None
       return load(RUNDIR + fn)

In this case, we've written ``load_run`` to look in a specific directory to see
if it can find an output with the given name.  So now we can write scripts that
use this function:

.. code-block:: python

   import yt.mods as yt

   my_run = yt.load_run("hotgasflow/DD0040/DD0040")

And because we have imported from ``yt.mods`` we have access to the
``load_run`` function defined in our plugin file.

Adding Custom Colormaps
^^^^^^^^^^^^^^^^^^^^^^^

To add custom :ref:`colormaps` to your plugin file, you must use the
:func:`~yt.visualization.color_maps.make_colormap` function to generate a 
colormap of your choice and then add it to the plugin file.  You can see
an example of this in :ref:`custom-colormaps`.  Remember that you don't need
to prefix commands in your plugin file with ``yt.``, but you'll only be
able to access the colormaps when you load the ``yt.mods`` module, not simply
``yt``.
