.. _api-reference:

API Reference
=============

Plots and the Plotting Interface
--------------------------------

SlicePlot and ProjectionPlot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.visualization.plot_window.SlicePlot
   ~yt.visualization.plot_window.AxisAlignedSlicePlot
   ~yt.visualization.plot_window.OffAxisSlicePlot
   ~yt.visualization.plot_window.ProjectionPlot
   ~yt.visualization.plot_window.OffAxisProjectionPlot
   ~yt.visualization.plot_window.WindowPlotMPL
   ~yt.visualization.plot_window.PlotWindow

ProfilePlot and PhasePlot
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.visualization.profile_plotter.ProfilePlot
   ~yt.visualization.profile_plotter.PhasePlot
   ~yt.visualization.profile_plotter.PhasePlotMPL

Fixed Resolution Pixelization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.visualization.fixed_resolution.FixedResolutionBuffer
   ~yt.visualization.fixed_resolution.CylindricalFixedResolutionBuffer
   ~yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer
   ~yt.visualization.fixed_resolution.OffAxisProjectionFixedResolutionBuffer

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
   :toctree: generated/

   ~yt.data_objects.data_containers.YTSelectionContainer
   ~yt.data_objects.data_containers.YTSelectionContainer0D
   ~yt.data_objects.data_containers.YTSelectionContainer1D
   ~yt.data_objects.data_containers.YTSelectionContainer2D
   ~yt.data_objects.data_containers.YTSelectionContainer3D

Selection Objects
+++++++++++++++++

These objects are defined by some selection method or mechanism.  Most are
geometric.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.selection_data_containers.YTPointBase
   ~yt.data_objects.selection_data_containers.YTOrthoRayBase
   ~yt.data_objects.selection_data_containers.YTRayBase
   ~yt.data_objects.selection_data_containers.YTSliceBase
   ~yt.data_objects.selection_data_containers.YTCuttingPlaneBase
   ~yt.data_objects.selection_data_containers.YTDiskBase
   ~yt.data_objects.selection_data_containers.YTRegionBase
   ~yt.data_objects.selection_data_containers.YTDataCollectionBase
   ~yt.data_objects.selection_data_containers.YTSphereBase
   ~yt.data_objects.selection_data_containers.YTEllipsoidBase
   ~yt.data_objects.selection_data_containers.YTCutRegionBase
   ~yt.data_objects.grid_patch.AMRGridPatch

Construction Objects
++++++++++++++++++++

These objects typically require some effort to build.  Often this means
integrating through the simulation in some way, or creating some large or
expensive set of intermediate data.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.construction_data_containers.YTStreamlineBase
   ~yt.data_objects.construction_data_containers.YTQuadTreeProjBase
   ~yt.data_objects.construction_data_containers.YTCoveringGridBase
   ~yt.data_objects.construction_data_containers.YTArbitraryGridBase
   ~yt.data_objects.construction_data_containers.YTSmoothedCoveringGridBase
   ~yt.data_objects.construction_data_containers.YTSurfaceBase

Time Series Objects
^^^^^^^^^^^^^^^^^^^

These are objects that either contain and represent or operate on series of
datasets.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.time_series.DatasetSeries
   ~yt.data_objects.time_series.DatasetSeriesObject
   ~yt.data_objects.time_series.TimeSeriesQuantitiesContainer
   ~yt.data_objects.time_series.AnalysisTaskProxy

Geometry Handlers
-----------------

These objects generate an "index" into multiresolution data.

.. autosummary::
   :toctree: generated/

   ~yt.geometry.geometry_handler.Index
   ~yt.geometry.grid_geometry_handler.GridIndex
   ~yt.geometry.oct_geometry_handler.OctreeIndex
   ~yt.geometry.particle_geometry_handler.ParticleIndex
   ~yt.geometry.unstructured_mesh_handler.UnstructuredIndex

Units
-----

These classes and functions enable yt's symbolic unit handling system.

.. autosummary::
   :toctree: generated/

   yt.data_objects.static_output.Dataset.arr
   yt.data_objects.static_output.Dataset.quan
   ~yt.units.unit_object.Unit
   ~yt.units.unit_registry.UnitRegistry
   ~yt.units.yt_array.YTArray
   ~yt.units.yt_array.YTQuantity


Frontends
---------

.. autosummary::
   :toctree: generated/

ARTIO
^^^^^

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.frontends.athena.data_structures.AthenaGrid
   ~yt.frontends.athena.data_structures.AthenaHierarchy
   ~yt.frontends.athena.data_structures.AthenaDataset
   ~yt.frontends.athena.fields.AthenaFieldInfo
   ~yt.frontends.athena.io.IOHandlerAthena

Boxlib
^^^^^^

.. autosummary::
   :toctree: generated/

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
   ~yt.frontends.boxlib.io.IOHandlerCastro
   ~yt.frontends.boxlib.io.IOHandlerNyx
   ~yt.frontends.boxlib.io.IOHandlerOrion

Chombo
^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.frontends.chombo.data_structures.ChomboGrid
   ~yt.frontends.chombo.data_structures.ChomboHierarchy
   ~yt.frontends.chombo.data_structures.ChomboDataset
   ~yt.frontends.chombo.data_structures.Orion2Hierarchy
   ~yt.frontends.chombo.data_structures.Orion2Dataset
   ~yt.frontends.chombo.io.IOHandlerChomboHDF5
   ~yt.frontends.chombo.io.IOHandlerChombo2DHDF5
   ~yt.frontends.chombo.io.IOHandlerChombo1DHDF5
   ~yt.frontends.chombo.io.IOHandlerOrion2HDF5

Enzo
^^^^

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.frontends.fits.data_structures.FITSGrid
   ~yt.frontends.fits.data_structures.FITSHierarchy
   ~yt.frontends.fits.data_structures.FITSDataset
   ~yt.frontends.fits.fields.FITSFieldInfo
   ~yt.frontends.fits.io.IOHandlerFITS

FLASH
^^^^^

.. autosummary::
   :toctree: generated/
   
   ~yt.frontends.flash.data_structures.FLASHGrid
   ~yt.frontends.flash.data_structures.FLASHHierarchy
   ~yt.frontends.flash.data_structures.FLASHDataset
   ~yt.frontends.flash.fields.FLASHFieldInfo
   ~yt.frontends.flash.io.IOHandlerFLASH

GDF
^^^

.. autosummary::
   :toctree: generated/

   ~yt.frontends.gdf.data_structures.GDFGrid
   ~yt.frontends.gdf.data_structures.GDFHierarchy
   ~yt.frontends.gdf.data_structures.GDFDataset
   ~yt.frontends.gdf.io.IOHandlerGDFHDF5

Halo Catalogs
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.frontends.halo_catalog.data_structures.HaloCatalogHDF5File
   ~yt.frontends.halo_catalog.data_structures.HaloCatalogDataset
   ~yt.frontends.halo_catalog.fields.HaloCatalogFieldInfo
   ~yt.frontends.halo_catalog.io.IOHandlerHaloCatalogHDF5
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
   :toctree: generated/

   ~yt.frontends.moab.data_structures.MoabHex8Hierarchy
   ~yt.frontends.moab.data_structures.MoabHex8Mesh
   ~yt.frontends.moab.data_structures.MoabHex8Dataset
   ~yt.frontends.moab.data_structures.PyneHex8Mesh
   ~yt.frontends.moab.data_structures.PyneMeshHex8Hierarchy
   ~yt.frontends.moab.data_structures.PyneMoabHex8Dataset
   ~yt.frontends.moab.io.IOHandlerMoabH5MHex8
   ~yt.frontends.moab.io.IOHandlerMoabPyneHex8

RAMSES
^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.frontends.ramses.data_structures.RAMSESDomainFile
   ~yt.frontends.ramses.data_structures.RAMSESDomainSubset
   ~yt.frontends.ramses.data_structures.RAMSESIndex
   ~yt.frontends.ramses.data_structures.RAMSESDataset
   ~yt.frontends.ramses.fields.RAMSESFieldInfo
   ~yt.frontends.ramses.io.IOHandlerRAMSES

SPH and Particle Codes
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

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

Loading Data
------------

.. autosummary::
   :toctree: generated/

   ~yt.convenience.load
   ~yt.convenience.simulation
   ~yt.frontends.stream.data_structures.load_uniform_grid
   ~yt.frontends.stream.data_structures.load_amr_grids
   ~yt.frontends.stream.data_structures.load_particles
   ~yt.frontends.stream.data_structures.load_hexahedral_mesh

Derived Datatypes
-----------------

Profiles and Histograms
^^^^^^^^^^^^^^^^^^^^^^^

These types are used to sum data up and either return that sum or return an
average.  Typically they are more easily used through the ``ProfilePlot``
``PhasePlot`` interface. We also provide the ``create_profile`` function
to create these objects in a uniform manner.


.. autosummary::
   :toctree: generated/

   ~yt.data_objects.profiles.ProfileND
   ~yt.data_objects.profiles.Profile1D
   ~yt.data_objects.profiles.Profile2D
   ~yt.data_objects.profiles.Profile3D
   ~yt.data_objects.profiles.create_profile

.. _halo_analysis_ref:

Halo Analysis
^^^^^^^^^^^^^

The ``HaloCatalog`` object is the primary means for performing custom analysis 
on cosmological halos.  It is also the primary interface for halo finding.

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_analysis.halo_catalog.HaloCatalog
   ~yt.analysis_modules.halo_analysis.halo_finding_methods.HaloFindingMethod
   ~yt.analysis_modules.halo_analysis.halo_callbacks.HaloCallback
   ~yt.analysis_modules.halo_analysis.halo_callbacks.delete_attribute
   ~yt.analysis_modules.halo_analysis.halo_callbacks.halo_sphere
   ~yt.analysis_modules.halo_analysis.halo_callbacks.iterative_center_of_mass
   ~yt.analysis_modules.halo_analysis.halo_callbacks.load_profiles
   ~yt.analysis_modules.halo_analysis.halo_callbacks.phase_plot
   ~yt.analysis_modules.halo_analysis.halo_callbacks.profile
   ~yt.analysis_modules.halo_analysis.halo_callbacks.save_profiles
   ~yt.analysis_modules.halo_analysis.halo_callbacks.sphere_bulk_velocity
   ~yt.analysis_modules.halo_analysis.halo_callbacks.sphere_field_max_recenter
   ~yt.analysis_modules.halo_analysis.halo_callbacks.virial_quantities
   ~yt.analysis_modules.halo_analysis.halo_filters.HaloFilter
   ~yt.analysis_modules.halo_analysis.halo_filters.not_subhalo
   ~yt.analysis_modules.halo_analysis.halo_filters.quantity_value
   ~yt.analysis_modules.halo_analysis.halo_quantities.HaloQuantity
   ~yt.analysis_modules.halo_analysis.halo_quantities.bulk_velocity
   ~yt.analysis_modules.halo_analysis.halo_quantities.center_of_mass

Halo Finding
^^^^^^^^^^^^

These provide direct access to the halo finders.  However, it is strongly recommended 
to use the ``HaloCatalog``.

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_finding.halo_objects.FOFHaloFinder
   ~yt.analysis_modules.halo_finding.halo_objects.HOPHaloFinder
   ~yt.analysis_modules.halo_finding.rockstar.rockstar.RockstarHaloFinder

Two Point Functions
^^^^^^^^^^^^^^^^^^^

These functions are designed to create correlations or other results of
operations acting on two spatially-distinct points in a data source.  See also
:ref:`two_point_functions`.


.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.two_point_functions.two_point_functions.TwoPointFunctions
   ~yt.analysis_modules.two_point_functions.two_point_functions.FcnSet

Field Types
-----------

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.fields.field_info_container.FieldInfoContainer.add_field
   ~yt.data_objects.static_output.Dataset.add_field

Image Handling
--------------

For volume renderings and fixed resolution buffers the image object returned is
an ``ImageArray`` object, which has useful functions for image saving and 
writing to bitmaps.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.image_array.ImageArray

Extension Types
---------------

Cosmology, Star Particle Analysis, and Simulated Observations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the generation of stellar SEDs.  (See also :ref:`star_analysis`.)


.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.star_analysis.sfr_spectrum.StarFormationRate
   ~yt.analysis_modules.star_analysis.sfr_spectrum.SpectrumBuilder

Light cone generation and simulation analysis.  (See also
:ref:`light-cone-generator`.)


.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.cosmological_observation.light_cone.light_cone.LightCone
   ~yt.analysis_modules.cosmological_observation.light_ray.light_ray.LightRay

Absorption and X-ray spectra and spectral lines:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.absorption_spectrum.absorption_spectrum.AbsorptionSpectrum
   ~yt.analysis_modules.spectral_integrator.spectral_frequency_integrator.EmissivityIntegrator
   ~yt.analysis_modules.spectral_integrator.spectral_frequency_integrator.add_xray_emissivity_field

Absorption spectra fitting:

.. autosummary:: 
   :toctree: generated/

   ~yt.analysis_modules.absorption_spectrum.absorption_spectrum_fit.generate_total_fit

Sunrise exporting:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.sunrise_export.sunrise_exporter.export_to_sunrise
   ~yt.analysis_modules.sunrise_export.sunrise_exporter.export_to_sunrise_from_halolist

RADMC-3D exporting:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.radmc3d_export.RadMC3DInterface.RadMC3DLayer
   ~yt.analysis_modules.radmc3d_export.RadMC3DInterface.RadMC3DWriter

Volume Rendering
^^^^^^^^^^^^^^^^

See also :ref:`volume_rendering`.

Here are the primary entry points:

.. autosummary::
   :toctree: generated/

   ~yt.visualization.volume_rendering.camera.Camera
   ~yt.visualization.volume_rendering.camera.off_axis_projection
   ~yt.visualization.volume_rendering.camera.allsky_projection

These objects set up the way the image looks:

.. autosummary::
   :toctree: generated/

   ~yt.visualization.volume_rendering.transfer_functions.ColorTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.MultiVariateTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.PlanckTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.ProjectionTransferFunction
   ~yt.visualization.volume_rendering.transfer_functions.TransferFunction

There are also advanced objects for particular use cases:

.. autosummary::
   :toctree: generated/

   ~yt.visualization.volume_rendering.camera.MosaicFisheyeCamera
   ~yt.visualization.volume_rendering.camera.FisheyeCamera
   ~yt.visualization.volume_rendering.camera.MosaicCamera
   ~yt.visualization.volume_rendering.camera.plot_allsky_healpix
   ~yt.visualization.volume_rendering.camera.PerspectiveCamera
   ~yt.utilities.amr_kdtree.amr_kdtree.AMRKDTree
   ~yt.visualization.volume_rendering.camera.StereoPairCamera

Streamlining
^^^^^^^^^^^^

See also :ref:`streamlines`.


.. autosummary::
   :toctree: generated/

   ~yt.visualization.streamlines.Streamlines

Image Writing
^^^^^^^^^^^^^

These functions are all used for fast writing of images directly to disk,
without calling matplotlib.  This can be very useful for high-cadence outputs
where colorbars are unnecessary or for volume rendering.


.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.visualization.eps_writer.DualEPS
   ~yt.visualization.eps_writer.single_plot
   ~yt.visualization.eps_writer.multiplot
   ~yt.visualization.eps_writer.multiplot_yt
   ~yt.visualization.eps_writer.return_cmap

.. _derived-quantities-api:

Derived Quantities
------------------

See :ref:`derived-quantities`.


.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.visualization.plot_window.PWViewerMPL.annotate_clear
   ~yt.visualization.plot_modifications.ArrowCallback
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
   :toctree: generated/
   
   ~yt.visualization.color_maps.add_cmap
   ~yt.visualization.color_maps.make_colormap
   ~yt.visualization.color_maps.show_colormaps

Function List
-------------


.. autosummary::
   :toctree: generated/

   ~yt.convenience.load
   ~yt.data_objects.static_output.Dataset.all_data
   ~yt.data_objects.static_output.Dataset.box
   ~yt.funcs.deprecate
   ~yt.funcs.ensure_list
   ~yt.funcs.get_pbar
   ~yt.funcs.humanize_time
   ~yt.funcs.insert_ipython
   ~yt.funcs.is_root
   ~yt.funcs.iterable
   ~yt.funcs.just_one
   ~yt.funcs.only_on_root
   ~yt.funcs.paste_traceback
   ~yt.funcs.pdb_run
   ~yt.funcs.print_tb
   ~yt.funcs.rootonly
   ~yt.funcs.time_execution
   ~yt.analysis_modules.level_sets.contour_finder.identify_contours
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
   :toctree: generated/

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
   ~yt.utilities.math_utils.get_ortho_basis
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
   :toctree: generated/

   ~yt.config.YTConfigParser
   ~yt.utilities.parameter_file_storage.ParameterFileStore
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ObjectIterator
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelAnalysisInterface
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelObjectIterator


Testing Infrastructure
----------------------

The first set of functions are all provided by NumPy.

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   ~yt.testing.assert_rel_equal
   ~yt.testing.amrspace
   ~yt.testing.fake_random_ds
   ~yt.testing.expand_keywords
