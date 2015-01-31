.. _configuration-file:

Configuration File
==================

Configuration File Format
-------------------------

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
--------------------------------

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
-----------------------------------------

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
-------------------------------

The following external parameters are available.  A number of parameters are
used internally.

* ``coloredlogs`` (default: ``'False'``): Should logs be colored?
* ``loadfieldplugins`` (default: ``'True'``): Do we want to load the plugin file?
* ``pluginfilename``  (default ``'my_plugins.py'``) The name of our plugin file.
* ``logfile`` (default: ``'False'``): Should we output to a log file in the
  filesystem?
* ``loglevel`` (default: ``'20'``): What is the threshold (0 to 50) for outputting
  log files?
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