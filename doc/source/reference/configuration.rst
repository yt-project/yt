Customizing yt: The Configuration and Plugin Files
==================================================

yt features ways to customize it to your personal preferences in terms of
how much output it displays, loading custom fields, loading custom colormaps,
accessing test datasets regardless of where you are in the file system, etc.
This customization is done through :ref:`configuration-file` and
:ref:`plugin-file` both of which exist in your ``$HOME/.config/yt`` directory.

.. _configuration-file:

The Configuration
-----------------

The configuration is stored in simple text files (in the `toml <https://github.com/toml-lang/toml>`_ format).
The files allow to set internal yt variables to custom default values to be used in future sessions.
The configuration can either be stored :ref:`globally <global-conf>` or :ref:`locally <local-conf>`.

.. _global-conf:

Global Configuration
^^^^^^^^^^^^^^^^^^^^

If no local configuration file exists, yt will look for and recognize the file
``$HOME/.config/yt/yt.toml`` as a configuration file, containing several options
that can be modified and adjusted to control runtime behavior.  For example, a sample
``$HOME/.config/yt/yt.toml`` file could look
like:

.. code-block:: none

   [yt]
   log_level = 1
   maximum_stored_datasets = 10000

This configuration file would set the logging threshold much lower, enabling
much more voluminous output from yt.  Additionally, it increases the number of
datasets tracked between instantiations of yt. The configuration file can be
managed using the ``yt config --global`` helper. It can list, add, modify and remove
options from the configuration file, e.g.:

.. code-block:: none

   $ yt config -h
   $ yt config list
   $ yt config set yt log_level 1
   $ yt config rm yt maximum_stored_datasets


.. _local-conf:

Local Configuration
^^^^^^^^^^^^^^^^^^^

yt will look for a file named ``yt.toml`` in the current directory. If present, its options
are loaded and the global configuration is not read. Local configuration files
can contain the same options as the global one.

Local configuration files can either be edited manually, or alternatively they
can be managed using ``yt config --local``. It can list, add, modify and remove
options, and display the path to the local configuration file, e.g.:

.. code-block:: none

   $ yt config -h
   $ yt config list --local
   $ yt config set --local yt log_level 1
   $ yt config rm --local yt maximum_stored_datasets
   $ yt config print-path --local

If no local configuration file is present, these commands will create an (empty) one
in the current working directory.

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

   yt.set_log_level(1)

   ds = yt.load("my_data0001")
   ds.print_stats()

This has the same effect as setting ``log_level = 1`` in the configuration
file. Note that a log level of 1 means that all log messages are printed to
stdout.  To disable logging, set the log level to 50.


.. _config-options:

Available Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following external parameters are available.  A number of parameters are
used internally.

* ``colored_logs`` (default: ``False``): Should logs be colored?
* ``default_colormap`` (default: ``cmyt.arbre``): What colormap should be used by
  default for yt-produced images?
* ``plugin_filename``  (default ``my_plugins.py``) The name of our plugin file.
* ``log_level`` (default: ``20``): What is the threshold (0 to 50) for
  outputting log files?
* ``test_data_dir`` (default: ``/does/not/exist``): The default path the
  ``load()`` function searches for datasets when it cannot find a dataset in the
  current directory.
* ``reconstruct_index`` (default: ``True``): If true, grid edges for patch AMR
  datasets will be adjusted such that they fall as close as possible to an
  integer multiple of the local cell width. If you are working with a dataset
  with a large number of grids, setting this to False can speed up loading
  your dataset possibly at the cost of grid-aligned artifacts showing up in
  slice visualizations.
* ``notebook_password`` (default: empty): If set, this will be fed to the
  IPython notebook created by ``yt notebook``.  Note that this should be an
  sha512 hash, not a plaintext password.  Starting ``yt notebook`` with no
  setting will provide instructions for setting this.
* ``requires_ds_strict`` (default: ``True``): If true, answer tests wrapped
  with :func:`~yt.utilities.answer_testing.framework.requires_ds` will raise
  :class:`~yt.utilities.exceptions.YTUnidentifiedDataType` rather than consuming
  it if required dataset is not present.
* ``serialize`` (default: ``False``): If true, perform automatic
  :ref:`object serialization <object-serialization>`
* ``sketchfab_api_key`` (default: empty): API key for https://sketchfab.com/ for
  uploading AMRSurface objects.
* ``suppress_stream_logging`` (default: ``False``): If true, execution mode will be
  quiet.
* ``stdout_stream_logging`` (default: ``False``): If true, logging is directed
  to stdout rather than stderr
* ``skip_dataset_cache`` (default: ``False``): If true, automatic caching of datasets
  is turned off.
* ``supp_data_dir`` (default: ``/does/not/exist``): The default path certain
  submodules of yt look in for supplemental data files.


.. _per-field-config:

Available per-field Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to customize the default behaviour of plots using per-field configuration.
The default options for plotting a given field can be specified in the configuration file
in ``[plot.field_type.field_name]`` blocks. The available keys are

* ``cmap`` (default: ``yt.default_colormap``, see :ref:`config-options`): the colormap to
  use for the field.
* ``log`` (default: ``True``): use a log scale (or symlog if ``linthresh`` is also set).
* ``linthresh`` (default: ``None``): if set to a float different than ``None`` and ``log`` is
  ``True``, use a symlog normalization with the given linear threshold.
* ``units`` (defaults to the units of the field): the units to use to represent the field.
* ``path_length_units`` (default: ``cm``): the unit of the integration length when doing
  e.g. projections. This always has the dimensions of a length. Note that this will only
  be used if ``units`` is also set for the field. The final units will then be
  ``units*path_length_units``.

You can also set defaults for all fields of a given field type by omitting the field name,
as illustrated below in the deposit block.

.. code-block:: toml

  [plot.gas.density]
  cmap = "plasma"
  log = true
  units = "mp/cm**3"

  [plot.gas.velocity_divergence]
  cmap = "bwr"  # use a diverging colormap
  log = false   # and a linear scale

  [plot.deposit]
  path_length_units = "kpc"  # use kpc for deposition projections


.. _plugin-file:

Plugin Files
------------

Plugin files are a means of creating custom fields, quantities, data objects,
colormaps, and other code executable functions or classes to be used in future
yt sessions without modifying the source code directly.

To enable a plugin file, call the function
:func:`~yt.funcs.enable_plugins` at the top of your script.

Global system plugin file
^^^^^^^^^^^^^^^^^^^^^^^^^

yt will look for and recognize the file ``$HOME/.config/yt/my_plugins.py`` as a
plugin file. It is possible to rename this file to ``$HOME/.config/yt/<plugin_filename>.py``
by defining ``plugin_filename`` in your ``yt.toml`` file, as mentioned above.

.. note::

   You can tell that your system plugin file is being parsed by watching for a logging
   message when you import yt. Note that both the ``yt load`` and ``iyt``
   command line entry points parse the plugin file.


Local project plugin file
^^^^^^^^^^^^^^^^^^^^^^^^^

Optionally, :func:`~yt.funcs.enable_plugins` can be passed an argument to specify
a custom location for a plugin file. This can be useful to define project wise customizations.
In that use case, any system-level plugin file will be ignored.

Plugin File Format
^^^^^^^^^^^^^^^^^^

Plugin files should contain pure Python code. If accessing yt functions and classes
they will not require the ``yt.`` prefix, because of how they are loaded.

For example, if one created a plugin file containing:

.. code-block:: python

   def _myfunc(field, data):
       return np.random.random(data["density"].shape)


   add_field(
       "random",
       function=_myfunc,
       sampling_type="cell",
       dimensions="dimensionless",
       units="auto",
   )

then all of my data objects would have access to the field ``random``.

You can also define other convenience functions in your plugin file.  For
instance, you could define some variables or functions, and even import common
modules:

.. code-block:: python

   import os

   HOMEDIR = "/home/username/"
   RUNDIR = "/scratch/runs/"


   def load_run(fn):
       if not os.path.exists(RUNDIR + fn):
           return None
       return load(RUNDIR + fn)

In this case, we've written ``load_run`` to look in a specific directory to see
if it can find an output with the given name.  So now we can write scripts that
use this function:

.. code-block:: python

   import yt

   yt.enable_plugins()

   my_run = yt.load_run("hotgasflow/DD0040/DD0040")

And because we have used ``yt.enable_plugins`` we have access to the
``load_run`` function defined in our plugin file.

.. note::
    if your convenience function's name colliding with an existing object
    within yt's namespace, it will be ignored.

Note that using the plugins file implies that your script is no longer fully
reproducible. If you share your script with someone else and use some of the
functionality if your plugins file, you will also need to share your plugins
file for someone else to re-run your script properly.

Adding Custom Colormaps
^^^^^^^^^^^^^^^^^^^^^^^

To add custom :ref:`colormaps` to your plugin file, you must use the
:func:`~yt.visualization.color_maps.make_colormap` function to generate a
colormap of your choice and then add it to the plugin file. You can see
an example of this in :ref:`custom-colormaps`. Remember that you don't need
to prefix commands in your plugin file with ``yt.``, but you'll only be
able to access the colormaps when you load the ``yt.mods`` module, not simply
``yt``.
