.. _visualizing_particle_datasets_with_firefly:

Visualizing Particle Datasets with Firefly
==========================================
`Firefly <https://github.com/ageller/Firefly>`_
is an interactive, browser-based,
particle visualization platform that allows you to filter, colormap, and fly
through their data. The Python frontend allows users to both load in their
own datasets and customize every aspect of the user interface.
yt offers to ability
to export your data to Firefly's ffly or JSON format through the
:meth:`~yt.data_objects.data_containers.YTDataContainer.create_firefly_object`
method.

You can adjust the interface settings, particle colors, decimation factors, and
other `Firefly settings <https://ageller.github.io/Firefly/docs/build/html/index.html>`_
through the returned ``Firefly.reader`` object. Once the
settings are tuned to your liking, calling the ``reader.writeToDisk()`` method will
produce the final ffly files. Note that ``reader.clean_datadir`` defaults to true
when using
:meth:`~yt.data_objects.data_containers.YTDataContainer.create_firefly_object`
so if you would like to manage multiple datasets make sure to pass different
``datadir`` keyword arguments.

.. image:: _images/firefly_example.png
   :width: 85%
   :align: center
   :alt: Screenshot of a sample Firefly visualization

Exporting an Example Dataset to Firefly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here is an example of how to use yt to export data to Firefly using some
`sample data <https://yt-project.org/data/>`_.

.. code-block:: python

   ramses_ds = yt.load("DICEGalaxyDisk_nonCosmological/output_00002/info_00002.txt")

   region = ramses_ds.sphere(ramses_ds.domain_center, (1000, "kpc"))

   reader = region.create_firefly_object(
       "IsoGalaxyRamses",
       fields_to_include=["particle_extra_field_1", "particle_extra_field_2"],
       fields_units=["dimensionless", "dimensionless"],
   )

   ## adjust some of the options
   reader.settings["color"]["io"] = [1, 1, 0, 1]  ## set default color
   reader.particleGroups[0].decimation_factor = 100  ## increase the decimation factor

   ## dump files to
   ##  ~/IsoGalaxyRamses/Dataio000.ffly
   ##  ~/IsoGalaxyRamses/filenames.json
   ##  ~/IsoGalaxyRamses/DataSettings.json
   reader.writeToDisk()
