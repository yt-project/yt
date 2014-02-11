

.. class:: boolean(self, regions, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRBooleanRegionBase`.)


.. class:: covering_grid(self, level, left_edge, dims, fields=None, pf=None, num_ghost_zones=0, use_pbar=True, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRCoveringGridBase`.)


.. class:: cut_region(self, base_region, field_cuts, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.InLineExtractedRegionBase`.)


.. class:: cutting(self, normal, center, fields=None, node_name=None, north_vector=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRCuttingPlaneBase`.)


.. class:: disk(self, center, normal, radius, height, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRCylinderBase`.)


.. class:: ellipsoid(self, center, A, B, C, e0, tilt, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMREllipsoidBase`.)


.. class:: extracted_region(self, base_region, indices, force_refresh=True, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.ExtractedRegionBase`.)


.. class:: fixed_res_cutting(self, normal, center, width, dims, fields=None, node_name=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRFixedResCuttingPlaneBase`.)


.. class:: fixed_res_proj(self, axis, level, left_edge, dims, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRFixedResProjectionBase`.)


.. class:: grid_collection(self, center, grid_list, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRGridCollectionBase`.)


.. class:: grid_collection_max_level(self, center, max_level, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRMaxLevelCollectionBase`.)


.. class:: inclined_box(self, origin, box_vectors, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRInclinedBoxBase`.)


.. class:: ortho_ray(self, axis, coords, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMROrthoRayBase`.)


.. class:: overlap_proj(self, axis, field, weight_field=None, max_level=None, center=None, pf=None, source=None, node_name=None, field_cuts=None, preload_style='level', serialize=True, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRProjBase`.)


.. class:: periodic_region(self, center, left_edge, right_edge, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRPeriodicRegionBase`.)


.. class:: periodic_region_strict(self, center, left_edge, right_edge, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRPeriodicRegionStrictBase`.)


.. class:: proj(self, axis, field, weight_field=None, max_level=None, center=None, pf=None, source=None, node_name=None, field_cuts=None, preload_style=None, serialize=True, style='integrate', **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRQuadTreeProjBase`.)


.. class:: ray(self, start_point, end_point, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRRayBase`.)


.. class:: region(self, center, left_edge, right_edge, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRRegionBase`.)


.. class:: region_strict(self, center, left_edge, right_edge, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRRegionStrictBase`.)


.. class:: slice(self, axis, coord, fields=None, center=None, pf=None, node_name=False, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRSliceBase`.)


.. class:: smoothed_covering_grid(self, *args, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRSmoothedCoveringGridBase`.)


.. class:: sphere(self, center, radius, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRSphereBase`.)


.. class:: streamline(self, positions, length=1.0, fields=None, pf=None, **field_parameters):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRStreamlineBase`.)


.. class:: surface(self, data_source, surface_field, field_value):

   For more information, see :ref:`physical-object-api`
   (This is a proxy for :class:`~yt.data_objects.data_containers.AMRSurfaceBase`.)
