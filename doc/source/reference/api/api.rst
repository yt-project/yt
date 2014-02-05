API Reference
=============

Plots and the Plotting Interface
--------------------------------

SlicePlot and ProjectionPlot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.visualization.plot_window.SlicePlot
   ~yt.visualization.plot_window.OffAxisSlicePlot
   ~yt.visualization.plot_window.ProjectionPlot
   ~yt.visualization.plot_window.OffAxisProjectionPlot

PlotCollection
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated/

   ~yt.visualization.plot_collection.PlotCollection
   ~yt.visualization.plot_collection.PlotCollectionInteractive
   ~yt.visualization.fixed_resolution.FixedResolutionBuffer
   ~yt.visualization.fixed_resolution.ObliqueFixedResolutionBuffer
   ~yt.visualization.base_plot_types.get_multi_plot

Data Sources
------------

.. _physical-object-api:

Physical Objects
^^^^^^^^^^^^^^^^

These are the objects that act as physical selections of data, describing a
region in space.  These are not typically addressed directly; see
:ref:`available-objects` for more information.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.data_containers.AMRCoveringGridBase
   ~yt.data_objects.data_containers.AMRCuttingPlaneBase
   ~yt.data_objects.data_containers.AMRCylinderBase
   ~yt.data_objects.data_containers.AMRGridCollectionBase
   ~yt.data_objects.data_containers.AMRRayBase
   ~yt.data_objects.data_containers.AMROrthoRayBase
   ~yt.data_objects.data_containers.AMRStreamlineBase
   ~yt.data_objects.data_containers.AMRQuadTreeProjBase
   ~yt.data_objects.data_containers.AMRProjBase
   ~yt.data_objects.data_containers.AMRRegionBase
   ~yt.data_objects.data_containers.AMRPeriodicRegionBase
   ~yt.data_objects.data_containers.AMRRegionStrictBase
   ~yt.data_objects.data_containers.AMRSliceBase
   ~yt.data_objects.data_containers.AMRSmoothedCoveringGridBase
   ~yt.data_objects.data_containers.AMRSphereBase
   ~yt.data_objects.data_containers.AMRSurfaceBase
   ~yt.data_objects.data_containers.AMRBooleanRegionBase
   ~yt.data_objects.data_containers.ExtractedRegionBase
   ~yt.data_objects.data_containers.AMRSurfaceBase

Time Series Objects
^^^^^^^^^^^^^^^^^^^

These are objects that either contain and represent or operate on series of
datasets.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.time_series.TimeSeriesData
   ~yt.data_objects.time_series.TimeSeriesDataObject
   ~yt.data_objects.time_series.TimeSeriesQuantitiesContainer
   ~yt.data_objects.time_series.AnalysisTaskProxy

Frontends
---------

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.grid_patch.AMRGridPatch
   ~yt.data_objects.hierarchy.AMRHierarchy
   ~yt.data_objects.static_output.StaticOutput

Enzo
^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.enzo.data_structures.EnzoGrid
   ~yt.frontends.enzo.data_structures.EnzoHierarchy
   ~yt.frontends.enzo.data_structures.EnzoStaticOutput

Orion
^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.orion.data_structures.OrionGrid
   ~yt.frontends.orion.data_structures.OrionHierarchy
   ~yt.frontends.orion.data_structures.OrionStaticOutput

FLASH
^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.flash.data_structures.FLASHGrid
   ~yt.frontends.flash.data_structures.FLASHHierarchy
   ~yt.frontends.flash.data_structures.FLASHStaticOutput

Chombo
^^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.chombo.data_structures.ChomboGrid
   ~yt.frontends.chombo.data_structures.ChomboHierarchy
   ~yt.frontends.chombo.data_structures.ChomboStaticOutput

Pluto
^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.pluto.data_structures.PlutoGrid
   ~yt.frontends.pluto.data_structures.PlutoHierarchy
   ~yt.frontends.pluto.data_structures.PlutoStaticOutput

Stream
^^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.stream.data_structures.StreamGrid
   ~yt.frontends.stream.data_structures.StreamHierarchy
   ~yt.frontends.stream.data_structures.StreamStaticOutput

Nyx
^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.nyx.data_structures.NyxGrid
   ~yt.frontends.nyx.data_structures.NyxHierarchy
   ~yt.frontends.nyx.data_structures.NyxStaticOutput

Athena
^^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.athena.data_structures.AthenaGrid
   ~yt.frontends.athena.data_structures.AthenaHierarchy
   ~yt.frontends.athena.data_structures.AthenaStaticOutput

GDF
^^^


.. autosummary::
   :toctree: generated/

   ~yt.frontends.gdf.data_structures.GDFGrid
   ~yt.frontends.gdf.data_structures.GDFHierarchy
   ~yt.frontends.gdf.data_structures.GDFStaticOutput

Derived Datatypes
-----------------

Profiles and Histograms
^^^^^^^^^^^^^^^^^^^^^^^

These types are used to sum data up and either return that sum or return an
average.  Typically they are more easily used through the
`yt.visualization.plot_collection` interface.


.. autosummary::
   :toctree: generated/

   ~yt.data_objects.profiles.BinnedProfile1D
   ~yt.data_objects.profiles.BinnedProfile2D
   ~yt.data_objects.profiles.BinnedProfile3D

Halo Finding and Particle Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Halo finding can be executed using these types.  Here we list the main halo
finders as well as a few other supplemental objects.

.. rubric:: Halo Finders

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_finding.halo_objects.FOFHaloFinder
   ~yt.analysis_modules.halo_finding.halo_objects.HOPHaloFinder
   ~yt.analysis_modules.halo_finding.halo_objects.parallelHF
   ~yt.analysis_modules.halo_finding.rockstar.api.RockstarHaloFinder

You can also operate on the Halo and HAloList objects themselves:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_finding.halo_objects.Halo
   ~yt.analysis_modules.halo_finding.halo_objects.HaloList
   ~yt.analysis_modules.halo_finding.halo_objects.HOPHalo
   ~yt.analysis_modules.halo_finding.halo_objects.RockstarHalo
   ~yt.analysis_modules.halo_finding.halo_objects.parallelHOPHalo
   ~yt.analysis_modules.halo_finding.halo_objects.FOFHalo
   ~yt.analysis_modules.halo_finding.halo_objects.LoadedHalo
   ~yt.analysis_modules.halo_finding.halo_objects.TextHalo
   ~yt.analysis_modules.halo_finding.halo_objects.RockstarHaloList
   ~yt.analysis_modules.halo_finding.halo_objects.HOPHaloList
   ~yt.analysis_modules.halo_finding.halo_objects.FOFHaloList
   ~yt.analysis_modules.halo_finding.halo_objects.LoadedHaloList
   ~yt.analysis_modules.halo_finding.halo_objects.TextHaloList
   ~yt.analysis_modules.halo_finding.halo_objects.parallelHOPHaloList

There are also functions for loading halos from disk:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_finding.halo_objects.LoadHaloes
   ~yt.analysis_modules.halo_finding.halo_objects.LoadTextHaloes
   ~yt.analysis_modules.halo_finding.halo_objects.LoadRockstarHalos

We have several methods that work to create merger trees:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_merger_tree.merger_tree.MergerTree
   ~yt.analysis_modules.halo_merger_tree.merger_tree.MergerTreeConnect
   ~yt.analysis_modules.halo_merger_tree.merger_tree.MergerTreeDotOutput
   ~yt.analysis_modules.halo_merger_tree.merger_tree.MergerTreeTextOutput

You can use Halo catalogs generatedl externally as well:

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_merger_tree.enzofof_merger_tree.HaloCatalog
   ~yt.analysis_modules.halo_merger_tree.enzofof_merger_tree.EnzoFOFMergerTree
   ~yt.analysis_modules.halo_merger_tree.enzofof_merger_tree.plot_halo_evolution

Halo Profiling
^^^^^^^^^^^^^^

yt provides a comprehensive halo profiler that can filter, center, and analyze
halos en masse.

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.halo_profiler.multi_halo_profiler.HaloProfiler
   ~yt.analysis_modules.halo_profiler.multi_halo_profiler.VirialFilter


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

   ~yt.data_objects.field_info_container.DerivedField
   ~yt.data_objects.field_info_container.FieldInfoContainer
   ~yt.data_objects.field_info_container.ValidateDataField
   ~yt.data_objects.field_info_container.ValidateGridType
   ~yt.data_objects.field_info_container.ValidateParameter
   ~yt.data_objects.field_info_container.ValidateProperty
   ~yt.data_objects.field_info_container.ValidateSpatial

Image Handling
--------------

For volume renderings and fixed resolution buffers the image object returned is
an ``ImageArray`` object, which has useful functions for image saving and 
writing to bitmaps.

.. autosummary::
   :toctree: generated/

   ~yt.data_objects.image_array.ImageArray
   ~yt.data_objects.image_array.ImageArray.write_png
   ~yt.data_objects.image_array.ImageArray.write_hdf5

Extension Types
---------------

Coordinate Transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^


.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.coordinate_transformation.transforms.arbitrary_regrid
   ~yt.analysis_modules.coordinate_transformation.transforms.spherical_regrid

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
   ~yt.analysis_modules.spectral_integrator.spectral_frequency_integrator.add_xray_luminosity_field
   ~yt.analysis_modules.spectral_integrator.spectral_frequency_integrator.add_xray_photon_emissivity_field

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

Radial Column Density
^^^^^^^^^^^^^^^^^^^^^

If you'd like to calculate the column density out to a given point, from a
specified center, yt can provide that information.

.. autosummary::
   :toctree: generated/

   ~yt.analysis_modules.radial_column_density.radial_column_density.RadialColumnDensity

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
   ~yt.visualization.image_writer.write_fits
   ~yt.visualization.image_writer.write_image
   ~yt.visualization.image_writer.map_to_colors
   ~yt.visualization.image_writer.strip_colormap_data
   ~yt.visualization.image_writer.splat_points
   ~yt.visualization.image_writer.annotate_image
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

   ~yt.data_objects.derived_quantities._AngularMomentumVector
   ~yt.data_objects.derived_quantities._BaryonSpinParameter
   ~yt.data_objects.derived_quantities._BulkVelocity
   ~yt.data_objects.derived_quantities._CenterOfMass
   ~yt.data_objects.derived_quantities._Extrema
   ~yt.data_objects.derived_quantities._IsBound
   ~yt.data_objects.derived_quantities._MaxLocation
   ~yt.data_objects.derived_quantities._ParticleSpinParameter
   ~yt.data_objects.derived_quantities._TotalMass
   ~yt.data_objects.derived_quantities._TotalQuantity
   ~yt.data_objects.derived_quantities._WeightedAverageQuantity

.. _callback-api:

Callback List
-------------


See also :ref:`callbacks`.

.. autosummary::
   :toctree: generated/

   ~yt.visualization.plot_modifications.ArrowCallback
   ~yt.visualization.plot_modifications.ClumpContourCallback
   ~yt.visualization.plot_modifications.ContourCallback
   ~yt.visualization.plot_modifications.CoordAxesCallback
   ~yt.visualization.plot_modifications.CuttingQuiverCallback
   ~yt.visualization.plot_modifications.GridBoundaryCallback
   ~yt.visualization.plot_modifications.HopCircleCallback
   ~yt.visualization.plot_modifications.HopParticleCallback
   ~yt.visualization.plot_modifications.LabelCallback
   ~yt.visualization.plot_modifications.LinePlotCallback
   ~yt.visualization.plot_modifications.MarkerAnnotateCallback
   ~yt.visualization.plot_modifications.ParticleCallback
   ~yt.visualization.plot_modifications.PointAnnotateCallback
   ~yt.visualization.plot_modifications.QuiverCallback
   ~yt.visualization.plot_modifications.SphereCallback
   ~yt.visualization.plot_modifications.TextLabelCallback
   ~yt.visualization.plot_modifications.TitleCallback
   ~yt.visualization.plot_modifications.UnitBoundaryCallback
   ~yt.visualization.plot_modifications.VelocityCallback

Function List
-------------


.. autosummary::
   :toctree: generated/

   ~yt.convenience.load
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
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_blocking_call
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_passthrough
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_root_only
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_simple_proxy
   ~yt.utilities.parallel_tools.parallel_analysis_interface.parallel_splitter

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
   ~yt.data_objects.data_containers.FakeGridForParticles
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ObjectIterator
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelAnalysisInterface
   ~yt.utilities.parallel_tools.parallel_analysis_interface.ParallelObjectIterator
   ~yt.analysis_modules.hierarchy_subset.hierarchy_subset.ConstructedRootGrid
   ~yt.analysis_modules.hierarchy_subset.hierarchy_subset.ExtractedHierarchy


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
   ~yt.testing.fake_random_pf
   ~yt.testing.expand_keywords
