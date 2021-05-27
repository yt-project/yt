.. _api-reference:

API Reference
=============

Plots and the Plotting Interface
--------------------------------

SlicePlot and ProjectionPlot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

   ~yt.visualization.plot_window.SlicePlot
   ~yt.visualization.plot_window.AxisAlignedSlicePlot
   ~yt.visualization.plot_window.OffAxisSlicePlot
   ~yt.visualization.plot_window.ProjectionPlot
   ~yt.visualization.plot_window.OffAxisProjectionPlot
   ~yt.visualization.plot_window.WindowPlotMPL
   ~yt.visualization.plot_window.PlotWindow
   ~yt.visualization.plot_window.plot_2d

ProfilePlot and PhasePlot
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

   ~yt.visualization.profile_plotter.ProfilePlot
   ~yt.visualization.profile_plotter.PhasePlot
   ~yt.visualization.profile_plotter.PhasePlotMPL

Particle Plots
^^^^^^^^^^^^^^

.. autosummary::

   ~yt.visualization.particle_plots.ParticleProjectionPlot
   ~yt.visualization.particle_plots.ParticlePhasePlot
   ~yt.visualization.particle_plots.ParticlePlot

Fixed Resolution Pixelization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

   ~yt.visualization.fixed_resolution.FixedResolutionBuffer
   ~yt.visualization.fixed_resolution.ParticleImageBuffer
   ~yt.visualization.fixed_resolution.CylindricalFixedResolutionBuffer
   ~yt.visualization.fixed_resolution.OffAxisProjectionFixedResolutionBuffer

Writing FITS images
^^^^^^^^^^^^^^^^^^^

.. autosummary::

   ~yt.visualization.fits_image.FITSImageData
   ~yt.visualization.fits_image.FITSSlice
   ~yt.visualization.fits_image.FITSProjection
   ~yt.visualization.fits_image.FITSOffAxisSlice
   ~yt.visualization.fits_image.FITSOffAxisProjection

Data Sources
------------

.. _physical-object-api:

Physical Objects
^^^^^^^^^^^^^^^^

These are the objects that act as physical selections of data, describing a
region in space.  These are not typically addressed directly; see
:ref:`available-objects` for more information.

Base Classes
++++++++++++

These will almost never need to be instantiated on their own.

.. autosummary::

   ~yt.data_objects.data_containers.YTDataContainer
   ~yt.data_objects.selection_objects.data_selection_objects.YTSelectionContainer
   ~yt.data_objects.selection_objects.data_selection_objects.YTSelectionContainer0D
   ~yt.data_objects.selection_objects.data_selection_objects.YTSelectionContainer1D
   ~yt.data_objects.selection_objects.data_selection_objects.YTSelectionContainer2D
   ~yt.data_objects.selection_objects.data_selection_objects.YTSelectionContainer3D

Selection Objects
+++++++++++++++++

These objects are defined by some selection method or mechanism.  Most are
geometric.

.. autosummary::

   ~yt.data_objects.selection_objects.point.YTPoint
   ~yt.data_objects.selection_objects.ray.YTOrthoRay
   ~yt.data_objects.selection_objects.ray.YTRay
   ~yt.data_objects.selection_objects.slices.YTSlice
   ~yt.data_objects.selection_objects.slices.YTCuttingPlane
   ~yt.data_objects.selection_objects.disk.YTDisk
   ~yt.data_objects.selection_objects.region.YTRegion
   ~yt.data_objects.selection_objects.object_collection.YTDataCollection
   ~yt.data_objects.selection_objects.spheroids.YTSphere
   ~yt.data_objects.selection_objects.spheroids.YTEllipsoid
   ~yt.data_objects.selection_objects.cur_region.YTCutRegion
   ~yt.data_objects.index_subobjects.grid_patch.AMRGridPatch
   ~yt.data_objects.index_subobjects.octree_subset.OctreeSubset
   ~yt.data_objects.index_subobjects.particle_container.ParticleContainer
   ~yt.data_objects.index_subobjects.unstructured_mesh.UnstructuredMesh
   ~yt.data_objects.index_subobjects.unstructured_mesh.SemiStructuredMesh

Construction Objects
++++++++++++++++++++

These objects typically require some effort to build.  Often this means
integrating through the simulation in some way, or creating some large or
expensive set of intermediate data.

.. autosummary::

   ~yt.data_objects.construction_data_containers.YTStreamline
   ~yt.data_objects.construction_data_containers.YTQuadTreeProj
   ~yt.data_objects.construction_data_containers.YTCoveringGrid
   ~yt.data_objects.construction_data_containers.YTArbitraryGrid
   ~yt.data_objects.construction_data_containers.YTSmoothedCoveringGrid
   ~yt.data_objects.construction_data_containers.YTSurface

Time Series Objects
^^^^^^^^^^^^^^^^^^^

These are objects that either contain and represent or operate on series of
datasets.

.. autosummary::

   ~yt.data_objects.time_series.DatasetSeries
   ~yt.data_objects.time_series.DatasetSeriesObject
   ~yt.data_objects.time_series.SimulationTimeSeries
   ~yt.data_objects.time_series.TimeSeriesQuantitiesContainer
   ~yt.data_objects.time_series.AnalysisTaskProxy
   ~yt.data_objects.particle_trajectories.ParticleTrajectories

Geometry Handlers
-----------------

These objects generate an "index" into multiresolution data.

.. autosummary::

   ~yt.geometry.geometry_handler.Index
   ~yt.geometry.grid_geometry_handler.GridIndex
   ~yt.geometry.oct_geometry_handler.OctreeIndex
   ~yt.geometry.particle_geometry_handler.ParticleIndex
   ~yt.geometry.unstructured_mesh_handler.UnstructuredIndex

Units
-----

These classes and functions enable yt's symbolic unit handling system.

.. autosummary::

   yt.data_objects.static_output.Dataset.arr
   yt.data_objects.static_output.Dataset.quan
   ~yt.units.unit_object.define_unit
   ~yt.units.unit_object.Unit
   ~yt.units.unit_registry.UnitRegistry
   ~yt.units.unit_systems.UnitSystem
   ~yt.units.yt_array.YTArray
   ~yt.units.yt_array.YTQuantity
   ~yt.units.yt_array.uconcatenate
   ~yt.units.yt_array.uintersect1d
   ~yt.units.yt_array.uunion1d
   ~yt.units.yt_array.unorm
   ~yt.units.yt_array.udot
   ~yt.units.yt_array.uvstack
   ~yt.units.yt_array.uhstack
   ~yt.units.yt_array.ustack
   ~yt.units.yt_array.display_ytarray

Frontends
---------

.. autosummary::

AMRVAC
^^^^^^

.. autosummary::

   ~yt.frontends.amrvac.data_structures.AMRVACGrid
   ~yt.frontends.amrvac.data_structures.AMRVACHierarchy
   ~yt.frontends.amrvac.data_structures.AMRVACDataset
   ~yt.frontends.amrvac.fields.AMRVACFieldInfo
   ~yt.frontends.amrvac.io.AMRVACIOHandler
   ~yt.frontends.amrvac.io.read_amrvac_namelist

ARTIO
^^^^^

.. autosummary::

   ~yt.frontends.artio.data_structures.ARTIOIndex
   ~yt.frontends.artio.data_structures.ARTIOOctreeSubset
   ~yt.frontends.artio.data_structures.ARTIORootMeshSubset
   ~yt.frontends.artio.data_structures.ARTIODataset
   ~yt.frontends.artio.definitions.ARTIOconstants
   ~yt.frontends.artio.fields.ARTIOFieldInfo
   ~yt.frontends.artio.io.IOHandlerARTIO


Athena
^^^^^^

.. autosummary::

   ~yt.frontends.athena.data_structures.AthenaGrid
   ~yt.frontends.athena.data_structures.AthenaHierarchy
   ~yt.frontends.athena.data_structures.AthenaDataset
   ~yt.frontends.athena.fields.AthenaFieldInfo
   ~yt.frontends.athena.io.IOHandlerAthena

AMReX/Boxlib
^^^^^^^^^^^^

.. autosummary::

   ~yt.frontends.boxlib.data_structures.BoxlibGrid
   ~yt.frontends.boxlib.data_structures.BoxlibHierarchy
   ~yt.frontends.boxlib.data_structures.BoxlibDataset
   ~yt.frontends.boxlib.data_structures.CastroDataset
   ~yt.frontends.boxlib.data_structures.MaestroDataset
   ~yt.frontends.boxlib.data_structures.NyxHierarchy
   ~yt.frontends.boxlib.data_structures.NyxDataset
   ~yt.frontends.boxlib.data_structures.OrionHierarchy
   ~yt.frontends.boxlib.data_structures.OrionDataset
   ~yt.frontends.boxlib.fields.BoxlibFieldInfo
   ~yt.frontends.boxlib.io.IOHandlerBoxlib
   ~yt.frontends.boxlib.io.IOHandlerOrion

Chombo
^^^^^^

.. autosummary::

   ~yt.frontends.chombo.data_structures.ChomboGrid
   ~yt.frontends.chombo.data_structures.ChomboHierarchy
   ~yt.frontends.chombo.data_structures.ChomboDataset
   ~yt.frontends.chombo.data_structures.Orion2Hierarchy
   ~yt.frontends.chombo.data_structures.Orion2Dataset
   ~yt.frontends.chombo.io.IOHandlerChomboHDF5
   ~yt.frontends.chombo.io.IOHandlerOrion2HDF5

Enzo
^^^^

.. autosummary::

   ~yt.frontends.enzo.answer_testing_support.ShockTubeTest
   ~yt.frontends.enzo.data_structures.EnzoGrid
   ~yt.frontends.enzo.data_structures.EnzoGridGZ
   ~yt.frontends.enzo.data_structures.EnzoGridInMemory
   ~yt.frontends.enzo.data_structures.EnzoHierarchy1D
   ~yt.frontends.enzo.data_structures.EnzoHierarchy2D
   ~yt.frontends.enzo.data_structures.EnzoHierarchy
   ~yt.frontends.enzo.data_structures.EnzoHierarchyInMemory
   ~yt.frontends.enzo.data_structures.EnzoDatasetInMemory
   ~yt.frontends.enzo.data_structures.EnzoDataset
   ~yt.frontends.enzo.fields.EnzoFieldInfo
   ~yt.frontends.enzo.io.IOHandlerInMemory
   ~yt.frontends.enzo.io.IOHandlerPacked1D
   ~yt.frontends.enzo.io.IOHandlerPacked2D
   ~yt.frontends.enzo.io.IOHandlerPackedHDF5
   ~yt.frontends.enzo.io.IOHandlerPackedHDF5GhostZones
   ~yt.frontends.enzo.simulation_handling.EnzoCosmology
   ~yt.frontends.enzo.simulation_handling.EnzoSimulation

FITS
^^^^

.. autosummary::

   ~yt.frontends.fits.data_structures.FITSGrid
   ~yt.frontends.fits.data_structures.FITSHierarchy
   ~yt.frontends.fits.data_structures.FITSDataset
   ~yt.frontends.fits.fields.FITSFieldInfo
   ~yt.frontends.fits.io.IOHandlerFITS

FLASH
^^^^^

.. autosummary::

   ~yt.frontends.flash.data_structures.FLASHGrid
   ~yt.frontends.flash.data_structures.FLASHHierarchy
   ~yt.frontends.flash.data_structures.FLASHDataset
   ~yt.frontends.flash.fields.FLASHFieldInfo
   ~yt.frontends.flash.io.IOHandlerFLASH

GDF
^^^

.. autosummary::

   ~yt.frontends.gdf.data_structures.GDFGrid
   ~yt.frontends.gdf.data_structures.GDFHierarchy
   ~yt.frontends.gdf.data_structures.GDFDataset
   ~yt.frontends.gdf.io.IOHandlerGDFHDF5

Halo Catalogs
^^^^^^^^^^^^^

.. autosummary::

   ~yt.frontends.ahf.data_structures.AHFHalosDataset
   ~yt.frontends.ahf.fields.AHFHalosFieldInfo
   ~yt.frontends.ahf.io.IOHandlerAHFHalos
   ~yt.frontends.gadget_fof.data_structures.GadgetFOFDataset
   ~yt.frontends.gadget_fof.data_structures.GadgetFOFHDF5File
   ~yt.frontends.gadget_fof.data_structures.GadgetFOFHaloDataset
   ~yt.frontends.gadget_fof.io.IOHandlerGadgetFOFHDF5
   ~yt.frontends.gadget_fof.io.IOHandlerGadgetFOFHaloHDF5
   ~yt.frontends.gadget_fof.fields.GadgetFOFFieldInfo
   ~yt.frontends.gadget_fof.fields.GadgetFOFHaloFieldInfo
   ~yt.frontends.halo_catalog.data_structures.YTHaloCatalogFile
   ~yt.frontends.halo_catalog.data_structures.YTHaloCatalogDataset
   ~yt.frontends.halo_catalog.fields.YTHaloCatalogFieldInfo
   ~yt.frontends.halo_catalog.io.IOHandlerYTHaloCatalog
   ~yt.frontends.owls_subfind.data_structures.OWLSSubfindParticleIndex
   ~yt.frontends.owls_subfind.data_structures.OWLSSubfindHDF5File
   ~yt.frontends.owls_subfind.data_structures.OWLSSubfindDataset
   ~yt.frontends.owls_subfind.fields.OWLSSubfindFieldInfo
   ~yt.frontends.owls_subfind.io.IOHandlerOWLSSubfindHDF5
   ~yt.frontends.rockstar.data_structures.RockstarBinaryFile
   ~yt.frontends.rockstar.data_structures.RockstarDataset
   ~yt.frontends.rockstar.fields.RockstarFieldInfo
   ~yt.frontends.rockstar.io.IOHandlerRockstarBinary

MOAB
^^^^

.. autosummary::

   ~yt.frontends.moab.data_structures.MoabHex8Hierarchy
   ~yt.frontends.moab.data_structures.MoabHex8Mesh
   ~yt.frontends.moab.data_structures.MoabHex8Dataset
   ~yt.frontends.moab.data_structures.PyneHex8Mesh
   ~yt.frontends.moab.data_structures.PyneMeshHex8Hierarchy
   ~yt.frontends.moab.data_structures.PyneMoabHex8Dataset
   ~yt.frontends.moab.io.IOHandlerMoabH5MHex8
   ~yt.frontends.moab.io.IOHandlerMoabPyneHex8

OpenPMD
^^^^^^^

.. autosummary::

   ~yt.frontends.open_pmd.data_structures.OpenPMDGrid
   ~yt.frontends.open_pmd.data_structures.OpenPMDHierarchy
   ~yt.frontends.open_pmd.data_structures.OpenPMDDataset
   ~yt.frontends.open_pmd.fields.OpenPMDFieldInfo
   ~yt.frontends.open_pmd.io.IOHandlerOpenPMDHDF5
   ~yt.frontends.open_pmd.misc.parse_unit_dimension
   ~yt.frontends.open_pmd.misc.is_const_component
   ~yt.frontends.open_pmd.misc.get_component

RAMSES
^^^^^^

.. autosummary::

   ~yt.frontends.ramses.data_structures.RAMSESDomainFile
   ~yt.frontends.ramses.data_structures.RAMSESDomainSubset
   ~yt.frontends.ramses.data_structures.RAMSESIndex
   ~yt.frontends.ramses.data_structures.RAMSESDataset
   ~yt.frontends.ramses.fields.RAMSESFieldInfo
   ~yt.frontends.ramses.io.IOHandlerRAMSES

SPH and Particle Codes
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

   ~yt.frontends.gadget.data_structures.GadgetBinaryFile
   ~yt.frontends.gadget.data_structures.GadgetHDF5Dataset
   ~yt.frontends.gadget.data_structures.GadgetDataset
   ~yt.frontends.http_stream.data_structures.HTTPParticleFile
   ~yt.frontends.http_stream.data_structures.HTTPStreamDataset
   ~yt.frontends.owls.data_structures.OWLSDataset
   ~yt.frontends.sph.data_structures.ParticleDataset
   ~yt.frontends.tipsy.data_structures.TipsyFile
   ~yt.frontends.tipsy.data_structures.TipsyDataset
   ~yt.frontends.sph.fields.SPHFieldInfo
   ~yt.frontends.gadget.io.IOHandlerGadgetBinary
   ~yt.frontends.gadget.io.IOHandlerGadgetHDF5
   ~yt.frontends.http_stream.io.IOHandlerHTTPStream
   ~yt.frontends.owls.io.IOHandlerOWLS
   ~yt.frontends.tipsy.io.IOHandlerTipsyBinary

Stream
^^^^^^

.. autosummary::

   ~yt.frontends.stream.data_structures.StreamDictFieldHandler
   ~yt.frontends.stream.data_structures.StreamGrid
   ~yt.frontends.stream.data_structures.StreamHandler
   ~yt.frontends.stream.data_structures.StreamHexahedralHierarchy
   ~yt.frontends.stream.data_structures.StreamHexahedralMesh
   ~yt.frontends.stream.data_structures.StreamHexahedralDataset
   ~yt.frontends.stream.data_structures.StreamHierarchy
   ~yt.frontends.stream.data_structures.StreamOctreeHandler
   ~yt.frontends.stream.data_structures.StreamOctreeDataset
   ~yt.frontends.stream.data_structures.StreamOctreeSubset
   ~yt.frontends.stream.data_structures.StreamParticleFile
   ~yt.frontends.stream.data_structures.StreamParticleIndex
   ~yt.frontends.stream.data_structures.StreamParticlesDataset
   ~yt.frontends.stream.data_structures.StreamDataset
   ~yt.frontends.stream.fields.StreamFieldInfo
   ~yt.frontends.stream.io.IOHandlerStream
   ~yt.frontends.stream.io.IOHandlerStreamHexahedral
   ~yt.frontends.stream.io.IOHandlerStreamOctree
   ~yt.frontends.stream.io.StreamParticleIOHandler

ytdata
^^^^^^

.. autosummary::

   ~yt.frontends.ytdata.data_structures.YTDataContainerDataset
   ~yt.frontends.ytdata.data_structures.YTSpatialPlotDataset
   ~yt.frontends.ytdata.data_structures.YTGridDataset
   ~yt.frontends.ytdata.data_structures.YTGridHierarchy
   ~yt.frontends.ytdata.data_structures.YTGrid
   ~yt.frontends.ytdata.data_structures.YTNonspatialDataset
   ~yt.frontends.ytdata.data_structures.YTNonspatialHierarchy
   ~yt.frontends.ytdata.data_structures.YTNonspatialGrid
   ~yt.frontends.ytdata.data_structures.YTProfileDataset
   ~yt.frontends.ytdata.data_structures.YTClumpTreeDataset
   ~yt.frontends.ytdata.data_structures.YTClumpContainer
   ~yt.frontends.ytdata.fields.YTDataContainerFieldInfo
   ~yt.frontends.ytdata.fields.YTGridFieldInfo
   ~yt.frontends.ytdata.io.IOHandlerYTDataContainerHDF5
   ~yt.frontends.ytdata.io.IOHandlerYTGridHDF5
   ~yt.frontends.ytdata.io.IOHandlerYTSpatialPlotHDF5
   ~yt.frontends.ytdata.io.IOHandlerYTNonspatialhdf5

Loading Data
------------

.. autosummary::

   ~yt.loaders.load
   ~yt.loaders.simulation
   ~yt.loaders.load_uniform_grid
   ~yt.loaders.load_amr_grids
   ~yt.loaders.load_particles
   ~yt.loaders.load_octree
   ~yt.loaders.load_hexahedral_mesh
   ~yt.loaders.load_unstructured_mesh
   ~yt.loaders.load_sample

Derived Datatypes
-----------------

Profiles and Histograms
^^^^^^^^^^^^^^^^^^^^^^^

These types are used to sum data up and either return that sum or return an
average.  Typically they are more easily used through the ``ProfilePlot``
``PhasePlot`` interface. We also provide the ``create_profile`` function
to create these objects in a uniform manner.


.. autosummary::

   ~yt.data_objects.profiles.ProfileND
   ~yt.data_objects.profiles.Profile1D
   ~yt.data_objects.profiles.Profile2D
   ~yt.data_objects.profiles.Profile3D
   ~yt.data_objects.profiles.ParticleProfile
   ~yt.data_objects.profiles.create_profile

.. _clump_finding_ref:

Clump Finding
^^^^^^^^^^^^^

The ``Clump`` object and associated functions can be used for identification
of topologically disconnected structures, i.e., clump finding.

.. autosummary::

   ~yt.data_objects.level_sets.clump_handling.Clump
   ~yt.data_objects.level_sets.clump_handling.Clump.add_info_item
   ~yt.data_objects.level_sets.clump_handling.Clump.add_validator
   ~yt.data_objects.level_sets.clump_handling.Clump.save_as_dataset
   ~yt.data_objects.level_sets.clump_handling.find_clumps
   ~yt.data_objects.level_sets.clump_info_items.add_clump_info
   ~yt.data_objects.level_sets.clump_validators.add_validator

X-ray Emission Fields
^^^^^^^^^^^^^^^^^^^^^

This can be used to create derived fields of X-ray emission in
different energy bands.

.. autosummary::

   ~yt.fields.xray_emission_fields.XrayEmissivityIntegrator
   ~yt.fields.xray_emission_fields.add_xray_emissivity_field

Field Types
-----------

.. autosummary::

   ~yt.fields.field_info_container.FieldInfoContainer
   ~yt.fields.derived_field.DerivedField
   ~yt.fields.derived_field.ValidateDataField
   ~yt.fields.derived_field.ValidateGridType
   ~yt.fields.derived_field.ValidateParameter
   ~yt.fields.derived_field.ValidateProperty
   ~yt.fields.derived_field.ValidateSpatial

Field Functions
---------------

.. autosummary::

   ~yt.fields.field_info_container.FieldInfoContainer.add_field
   ~yt.data_objects.static_output.Dataset.add_field
   ~yt.data_objects.static_outputs.add_deposited_particle_field
   ~yt.data_objects.static_outputs.add_mesh_sampling_particle_field
   ~yt.data_objects.static_outputs.add_smoothed_particle_field
   ~yt.data_objects.static_outputs.add_gradient_fields
   ~yt.frontends.stream.data_structures.add_SPH_fields

Particle Filters
----------------

.. autosummary::

   ~yt.data_objects.particle_filters.add_particle_filter
   ~yt.data_objects.particle_filters.particle_filter

Image Handling
--------------

For volume renderings and fixed resolution buffers the image object returned is
an ``ImageArray`` object, which has useful functions for image saving and
writing to bitmaps.

.. autosummary::

   ~yt.data_objects.image_array.ImageArray

Volume Rendering
^^^^^^^^^^^^^^^^

See also :ref:`volume_rendering`.

Here are the primary entry points and the main classes involved in the
Scene infrastructure:

.. autosummary::

   ~yt.visualization.volume_rendering.volume_rendering.volume_render
   ~yt.visualization.volume_rendering.volume_rendering.create_scene
   ~yt.visualization.volume_rendering.off_axis_projection.off_axis_projection
   ~yt.visualization.volume_rendering.scene.Scene
   ~yt.visualization.volume_rendering.camera.Camera
   ~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree

The different kinds of sources:

.. autosummary::

   ~yt.visualization.volume_rendering.render_source.RenderSource
   ~yt.visualization.volume_rendering.render_source.VolumeSource
   ~yt.visualization.volume_rendering.render_source.PointSource
   ~yt.visualization.volume_rendering.render_source.LineSource
   ~yt.visualization.volume_rendering.render_source.BoxSource
   ~yt.visualization.volume_rendering.render_source.GridSource
   ~yt.visualization.volume_rendering.render_source.CoordinateVectorSource
   ~yt.visualization.volume_rendering.render_source.MeshSource

The different kinds of transfer functions:

.. autosummary::

   ~yt.visualization.volume_rendering.transfer_functions.TransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.PlanckTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.MultiVariateTransferFunction
   ~yt.visualization.volume_rendering.transfer_function_helper.TransferFunctionHelper

The different kinds of lenses:

.. autosummary::

   ~yt.visualization.volume_rendering.lens.Lens
   ~yt.visualization.volume_rendering.lens.PlaneParallelLens
   ~yt.visualization.volume_rendering.lens.PerspectiveLens
   ~yt.visualization.volume_rendering.lens.StereoPerspectiveLens
   ~yt.visualization.volume_rendering.lens.FisheyeLens
   ~yt.visualization.volume_rendering.lens.SphericalLens
   ~yt.visualization.volume_rendering.lens.StereoSphericalLens

Streamlining
^^^^^^^^^^^^

See also :ref:`streamlines`.


.. autosummary::

   ~yt.visualization.streamlines.Streamlines

Image Writing
^^^^^^^^^^^^^

These functions are all used for fast writing of images directly to disk,
without calling matplotlib.  This can be very useful for high-cadence outputs
where colorbars are unnecessary or for volume rendering.


.. autosummary::

   ~yt.visualization.image_writer.multi_image_composite
   ~yt.visualization.image_writer.write_bitmap
   ~yt.visualization.image_writer.write_projection
   ~yt.visualization.image_writer.write_image
   ~yt.visualization.image_writer.map_to_colors
   ~yt.visualization.image_writer.strip_colormap_data
   ~yt.visualization.image_writer.splat_points
   ~yt.visualization.image_writer.scale_image

We also provide a module that is very good for generating EPS figures,
particularly with complicated layouts.

.. autosummary::

   ~yt.visualization.eps_writer.DualEPS
   ~yt.visualization.eps_writer.single_plot
   ~yt.visualization.eps_writer.multiplot
   ~yt.visualization.eps_writer.multiplot_yt
   ~yt.visualization.eps_writer.return_cmap
   ~yt.visualization.eps_writer.return_colormap

.. _derived-quantities-api:

Derived Quantities
------------------

See :ref:`derived-quantities`.


.. autosummary::

   ~yt.data_objects.derived_quantities.DerivedQuantity
   ~yt.data_objects.derived_quantities.DerivedQuantityCollection
   ~yt.data_objects.derived_quantities.WeightedAverageQuantity
   ~yt.data_objects.derived_quantities.AngularMomentumVector
   ~yt.data_objects.derived_quantities.BulkVelocity
   ~yt.data_objects.derived_quantities.CenterOfMass
   ~yt.data_objects.derived_quantities.Extrema
   ~yt.data_objects.derived_quantities.MaxLocation
   ~yt.data_objects.derived_quantities.MinLocation
   ~yt.data_objects.derived_quantities.SpinParameter
   ~yt.data_objects.derived_quantities.TotalMass
   ~yt.data_objects.derived_quantities.TotalQuantity
   ~yt.data_objects.derived_quantities.WeightedAverageQuantity
   ~yt.data_objects.derived_quantities.WeightedVariance

.. _callback-api:

Callback List
-------------


See also :ref:`callbacks`.

.. autosummary::

   ~yt.visualization.plot_window.PWViewerMPL.annotate_clear
   ~yt.visualization.plot_window.PWViewerMPL.clear_annotations
   ~yt.visualization.plot_modifications.ArrowCallback
   ~yt.visualization.plot_modifications.CellEdgesCallback
   ~yt.visualization.plot_modifications.ClumpContourCallback
   ~yt.visualization.plot_modifications.ContourCallback
   ~yt.visualization.plot_modifications.CuttingQuiverCallback
   ~yt.visualization.plot_modifications.GridBoundaryCallback
   ~yt.visualization.plot_modifications.HaloCatalogCallback
   ~yt.visualization.plot_modifications.ImageLineCallback
   ~yt.visualization.plot_modifications.LinePlotCallback
   ~yt.visualization.plot_modifications.MagFieldCallback
   ~yt.visualization.plot_modifications.MarkerAnnotateCallback
   ~yt.visualization.plot_modifications.ParticleCallback
   ~yt.visualization.plot_modifications.PointAnnotateCallback
   ~yt.visualization.plot_modifications.QuiverCallback
   ~yt.visualization.plot_modifications.RayCallback
   ~yt.visualization.plot_modifications.ScaleCallback
   ~yt.visualization.plot_modifications.SphereCallback
   ~yt.visualization.plot_modifications.StreamlineCallback
   ~yt.visualization.plot_modifications.TextLabelCallback
   ~yt.visualization.plot_modifications.TimestampCallback
   ~yt.visualization.plot_modifications.TitleCallback
   ~yt.visualization.plot_modifications.TriangleFacetsCallback
   ~yt.visualization.plot_modifications.VelocityCallback

Colormap Functions
------------------


See also :ref:`colormaps`.

.. autosummary::

   ~yt.visualization.color_maps.add_cmap
   ~yt.visualization.color_maps.add_colormap
   ~yt.visualization.color_maps.make_colormap
   ~yt.visualization.color_maps.show_colormaps

Function List
-------------


.. autosummary::

   ~yt.frontends.ytdata.utilities.save_as_dataset
   ~yt.data_objects.data_containers.YTDataContainer.save_as_dataset
   ~yt.data_objects.static_output.Dataset.all_data
   ~yt.data_objects.static_output.Dataset.box
   ~yt.funcs.enable_plugins
   ~yt.funcs.get_pbar
   ~yt.funcs.humanize_time
   ~yt.funcs.insert_ipython
   ~yt.funcs.is_root
   ~yt.funcs.is_sequence
   ~yt.funcs.iter_fields
   ~yt.funcs.just_one
   ~yt.funcs.only_on_root
   ~yt.funcs.paste_traceback
   ~yt.funcs.pdb_run
   ~yt.funcs.print_tb
   ~yt.funcs.rootonly
   ~yt.funcs.time_execution
   ~yt.data_objects.level_sets.contour_finder.identify_contours
   ~yt.utilities.parallel_tools.parallel_analysis_interface.enable_parallelism
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_blocking_call
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_objects
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_passthrough
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_root_only
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_simple_proxy
   ~yt.data_objects.data_containers.YTDataContainer.get_field_parameter
   ~yt.data_objects.data_containers.YTDataContainer.set_field_parameter

Math Utilities
--------------


.. autosummary::

   ~yt.utilities.math_utils.periodic_position
   ~yt.utilities.math_utils.periodic_dist
   ~yt.utilities.math_utils.euclidean_dist
   ~yt.utilities.math_utils.rotate_vector_3D
   ~yt.utilities.math_utils.modify_reference_frame
   ~yt.utilities.math_utils.compute_rotational_velocity
   ~yt.utilities.math_utils.compute_parallel_velocity
   ~yt.utilities.math_utils.compute_radial_velocity
   ~yt.utilities.math_utils.compute_cylindrical_radius
   ~yt.utilities.math_utils.ortho_find
   ~yt.utilities.math_utils.quartiles
   ~yt.utilities.math_utils.get_rotation_matrix
   ~yt.utilities.math_utils.get_sph_r
   ~yt.utilities.math_utils.resize_vector
   ~yt.utilities.math_utils.get_sph_theta
   ~yt.utilities.math_utils.get_sph_phi
   ~yt.utilities.math_utils.get_cyl_r
   ~yt.utilities.math_utils.get_cyl_z
   ~yt.utilities.math_utils.get_cyl_theta
   ~yt.utilities.math_utils.get_cyl_r_component
   ~yt.utilities.math_utils.get_cyl_theta_component
   ~yt.utilities.math_utils.get_cyl_z_component
   ~yt.utilities.math_utils.get_sph_r_component
   ~yt.utilities.math_utils.get_sph_phi_component
   ~yt.utilities.math_utils.get_sph_theta_component


Miscellaneous Types
-------------------


.. autosummary::

   ~yt.config.YTConfigParser
   ~yt.utilities.parameter_file_storage.ParameterFileStore
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ObjectIterator
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelAnalysisInterface
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelObjectIterator

.. _cosmology-calculator-ref:

Cosmology Calculator
--------------------

.. autosummary::

   ~yt.utilities.cosmology.Cosmology
   ~yt.utilities.cosmology.Cosmology.hubble_distance
   ~yt.utilities.cosmology.Cosmology.comoving_radial_distance
   ~yt.utilities.cosmology.Cosmology.comoving_transverse_distance
   ~yt.utilities.cosmology.Cosmology.comoving_volume
   ~yt.utilities.cosmology.Cosmology.angular_diameter_distance
   ~yt.utilities.cosmology.Cosmology.angular_scale
   ~yt.utilities.cosmology.Cosmology.luminosity_distance
   ~yt.utilities.cosmology.Cosmology.lookback_time
   ~yt.utilities.cosmology.Cosmology.hubble_time
   ~yt.utilities.cosmology.Cosmology.critical_density
   ~yt.utilities.cosmology.Cosmology.hubble_parameter
   ~yt.utilities.cosmology.Cosmology.expansion_factor
   ~yt.utilities.cosmology.Cosmology.z_from_t
   ~yt.utilities.cosmology.Cosmology.t_from_z
   ~yt.utilities.cosmology.Cosmology.get_dark_factor

Testing Infrastructure
----------------------

The first set of functions are all provided by NumPy.

.. autosummary::

   ~yt.testing.assert_array_equal
   ~yt.testing.assert_almost_equal
   ~yt.testing.assert_approx_equal
   ~yt.testing.assert_array_almost_equal
   ~yt.testing.assert_equal
   ~yt.testing.assert_array_less
   ~yt.testing.assert_string_equal
   ~yt.testing.assert_array_almost_equal_nulp
   ~yt.testing.assert_allclose
   ~yt.testing.assert_raises

These are yt-provided functions:

.. autosummary::

   ~yt.testing.assert_rel_equal
   ~yt.testing.amrspace
   ~yt.testing.fake_random_ds
   ~yt.testing.expand_keywords

These are for the pytest infrastructure:

.. autosummary::

    ~conftest.tempdir
    ~conftest.answer_file
    ~conftest.hashing
    ~yt.utilities.answer_testing.answer_tests.grid_hierarchy
    ~yt.utilities.answer_testing.answer_tests.parentage_relationships
    ~yt.utilities.answer_testing.answer_tests.grid_values
    ~yt.utilities.answer_testing.answer_tests.projection_values
    ~yt.utilities.answer_testing.answer_tests.field_values
    ~yt.utilities.answer_testing.answer_tests.pixelized_projection_values
    ~yt.utilities.answer_testing.answer_tests.simulated_halo_mass_function
    ~yt.utilities.answer_testing.answer_tests.analytic_halo_mass_function
    ~yt.utilities.answer_testing.answer_tests.small_patch_amr
    ~yt.utilities.answer_testing.answer_tests.big_patch_amr
    ~yt.utilities.answer_testing.answer_tests.generic_array
    ~yt.utilities.answer_testing.answer_tests.sph_answer
    ~yt.utilities.answer_testing.answer_tests.get_field_size_and_mean
    ~yt.utilities.answer_testing.answer_tests.plot_window_attribute
    ~yt.utilities.answer_testing.answer_tests.phase_plot_attribute
    ~yt.utilities.answer_testing.answer_tests.generic_image
    ~yt.utilities.answer_testing.answer_tests.axial_pixelization
    ~yt.utilities.answer_testing.answer_tests.light_cone_projection
    ~yt.utilities.answer_testing.answer_tests.extract_connected_sets
    ~yt.utilities.answer_testing.answer_tests.VR_image_comparison
