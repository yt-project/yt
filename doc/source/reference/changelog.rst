.. _changelog:

ChangeLog
=========

This is a non-comprehensive log of changes to yt over its many releases.

Contributors
------------

The `CREDITS file <https://github.com/yt-project/yt/blob/master/CREDITS>`_ 
contains the most up-to-date list of everyone who has contributed to the yt 
source code.

Version 3.5.0
-------------

Version 3.5.0 is the first major release of yt since August 2017. It includes 
328 pull requests from 41 contributors, including 22 new contributors.

Major Changes
^^^^^^^^^^^^^

- ``yt.analysis_modules`` has been deprecated in favor of the new
  ``yt_astro_analysis`` package. New features and new astronomy-specific
  analysis modules will go into ``yt_astro_analysis`` and importing from
  ``yt.analysis_modules`` will raise a noisy warning. We will remove
  ``yt.analysis_modules`` in a future release. See `PR 1938
  <https://github.com/yt-project/yt/pull/1938>`__.
- Vector fields and derived fields depending on vector fields have been
  systematically updated to account for a bulk correction field parameter. For
  example, for the velocity field, all derived fields that depend on velocity
  will now account for the ``"bulk_velocity"`` field parameter. In addition, we
  have defined ``"relative_velocity"`` and ``"relative_magnetic_field"`` fields
  that include the bulk correction. Both of these are vector fields, to access
  the components, use e.g. ``"relative_velocity_x"``. The
  ``"particle_position_relative"`` and ``"particle_velocity_relative"`` fields
  have been deprecated. See `PR 1693
  <https://github.com/yt-project/yt/pull/1693>`__ and `PR 2022
  <https://github.com/yt-project/yt/pull/2022>`__.
- Aliases to spatial fields with the ``"gas"`` field type will now be returned
  in the default unit system for the dataset. As an example the ``"x"`` field
  might resolve to the field tuples ``("index", "x")`` or ``("gas",
  "x")``. Accessing the former will return data in code units while the latter
  will return data in whatever unit system the dataset is configured to use
  (CGS, by default). This means that to ensure the units of a spatial field will
  always be consistent, one must access the field as a tuple, explicitly
  specifying the field type. Accessing a spatial field using a string field name
  may return data in either code units or the dataset's default unit system
  depending on the history of field accesses prior to accessing that field. In
  the future accessing fields using an ambiguous field name will raise an
  error. See `PR 1799 <https://github.com/yt-project/yt/pull/1799>`__ and `PR
  1850 <https://github.com/yt-project/yt/pull/1850>`__.
- The ``max_level`` and ``min_level`` attributes of yt data objects now
  correctly update the state of the underlying data objects when set. In
  addition we have added an example to the cookbook that shows how to downsample
  AMR data using this functionality. See `PR 1737
  <https://github.com/yt-project/yt/pull/1737>`__.
- It is now possible to customize the formatting of labels for ion species
  fields. Rather than using the default spectroscopic notation, one can call
  ``ds.set_field_label_format("ionization_label", "plus_minus")`` to use the
  more traditional notation where ionization state is indicated with ``+`` and
  ``-`` symbols. See `PR 1867 <https://github.com/yt-project/yt/pull/1867>`__.

Improvements to the RAMSES frontend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We would particularly like to recognize Corentin Cadiou for his tireless work over the past year on improving support for RAMSES and octree AMR data in yt.
  
- Added support for reading RAMSES sink particles. See `PR 1548
  <https://github.com/yt-project/yt/pull/1548>`__.
- Add support for the new self-describing Ramses particle output format. See `PR
  1616 <https://github.com/yt-project/yt/pull/1616>`__.
- It is now possible to restrict the domain of a loaded Ramses dataset by
  passing a ``bbox`` keyword argument to ``yt.load()``. If passed this
  corresponds to the coordinates of the top-left and bottom-right hand corner of
  the subvolume to load. Data outside the bounding box will be ignored. This is
  useful for loading very large Ramses datasets where yt currently has poor
  scaling. See `PR 1637 <https://github.com/yt-project/yt/pull/1637>`__.
- The Ramses ``"particle_birth_time"`` field now contains the time when star
  particles form in a simulation in CGS units, formerly these times were only
  accessible via the incorrectly named ``"particle_age"`` field in conformal
  units. Correspondingly the ``"particle_age"`` field has been deprecated. The
  conformal birth time is not available via the ``"conformal_birth_time``"
  field. See `PR 1649 <https://github.com/yt-project/yt/pull/1649>`__.
- Substantial performance improvement for reading RAMSES AMR data. See `PR 1671
  <https://github.com/yt-project/yt/pull/1671>`__.
- The RAMSES frontend will now produce less voluminous logging feedback when
  loading the dataset or reading data. This is particularly noticeable for very
  large datasets with many CPU files. See `PR 1738
  <https://github.com/yt-project/yt/pull/1738>`__.
- Avoid repeated parsing of RAMSES particle and RT descriptors. See `PR 1739
  <https://github.com/yt-project/yt/pull/1739>`__.
- Added support for reading the RAMSES gravitational potential field. See `PR
  1751 <https://github.com/yt-project/yt/pull/1751>`__.
- Add support for RAMSES datasets that use the ``groupsize`` feature. See `PR
  1769 <https://github.com/yt-project/yt/pull/1769>`__.
- Dramatically improve the overall performance of the RAMSES frontend. See `PR
  1771 <https://github.com/yt-project/yt/pull/1771>`__.

Additional Improvements
^^^^^^^^^^^^^^^^^^^^^^^  
  
- Added support for particle data in the Enzo-P frontend. See `PR 1490
  <https://github.com/yt-project/yt/pull/1490>`__.
- Added an ``equivalence`` keyword argument to ``YTArray.in_units()`` and
  ``YTArray.to()``. This makes it possible to specify an equivalence when
  converting data to a new unit. Also added ``YTArray.to_value()`` which allows
  convering to a new unit, then stripping off the units to return a plain numpy
  array. See `PR 1563 <https://github.com/yt-project/yt/pull/1563>`__.
- Rather than crashing, yt will now assume default values for cosmology
  parameters in Gadget HDF5 data if it cannot find the relevant header
  information. See `PR 1578
  <https://github.com/yt-project/yt/pull/1578>`__.
- Improve detection for OpenMP support at compile-time, including adding support
  for detecting OpenMP on Windows. See `PR 1591
  <https://github.com/yt-project/yt/pull/1591>`__, `PR 1695
  <https://github.com/yt-project/yt/pull/1695>`__ and `PR 1696
  <https://github.com/yt-project/yt/pull/1696>`__.
- Add support for 2D cylindrical data for most plot callbacks. See `PR 1598
  <https://github.com/yt-project/yt/pull/1598>`__.
- Particles outside the domain are now ignored by ``load_uniform_grid()`` and
  ``load_amr_grids()``. See `PR 1602
  <https://github.com/yt-project/yt/pull/1602>`__.
- Fix incorrect units for the Gadget internal energy field in cosmology
  simulations. See `PR 1611
  <https://github.com/yt-project/yt/pull/1611>`__.
- Add support for calculating covering grids in parallel. See `PR 1612
  <https://github.com/yt-project/yt/pull/1612>`__.
- The number of particles in a dataset loaded by the stream frontend (e.g. via
  ``load_uniform_grid``) no longer needs to be explicitly provided via the
  ``number_of_particles`` keyword argument, using the ``number_of_particles``
  keyword will now generate a deprecation warning. See `PR 1620
  <https://github.com/yt-project/yt/pull/1620>`__.
- Add support for non-cartesian GAMER data. See `PR 1622
  <https://github.com/yt-project/yt/pull/1622>`__.
- If a particle filter depends on another particle filter, both particle filters
  will be registered for a dataset if the dependent particle filter is
  registered with a dataset. See `PR 1624
  <https://github.com/yt-project/yt/pull/1624>`__.
- The ``save()`` method of the various yt plot objets now optionally can accept
  a tuple of strings instead of a string. If a tuple is supplied, the elments
  are joined with ``os.sep`` to form a path. See `PR 1630
  <https://github.com/yt-project/yt/pull/1630>`__.
- The quiver callback now accepts a ``plot_args`` keyword argument that allows
  passing keyword arguments to matplotlib to allow for customization of the
  quiver plot. See `PR 1636 <https://github.com/yt-project/yt/pull/1636>`__.
- Updates and improvements for the OpenPMD frontend. See `PR 1645
  <https://github.com/yt-project/yt/pull/1645>`__.
- The mapserver now works correctly under Python3 and has new features like a
  colormap selector and plotting multiple fields via layers. See `PR 1654
  <https://github.com/yt-project/yt/pull/1654>`__ and `PR 1668
  <https://github.com/yt-project/yt/pull/1668>`__.
- Substantial performance improvement for calculating the gravitational
  potential in the clump finder. See `PR 1684
  <https://github.com/yt-project/yt/pull/1684>`__.
- Added new methods to ``ProfilePlot``: ``set_xlabel()``, ``set_ylabel()``,
  ``annotate_title()``, and ``annotate_text()``. See `PR 1700
  <https://github.com/yt-project/yt/pull/1700>`__ and `PR 1705
  <https://github.com/yt-project/yt/pull/1705>`__.
- Speedup for parallel halo finding operation for the FOF and HOP halo
  finders. See `PR 1724 <https://github.com/yt-project/yt/pull/1724>`__.
- Add support for halo finding using the rockstar halo finder on Python3. See
  `PR 1740 <https://github.com/yt-project/yt/pull/1740>`__.
- The ``ValidateParameter`` field validator has gained the ability for users to
  explicitly specify the values of field parameters during field detection. This
  makes it possible to write fields that access different sets of fields
  depending on the value of the field parameter. For example, a field might
  define an ``'axis'`` field parameter that can be either ``'x'``, ``'y'`` or
  ``'z'``. One can now explicitly tell the field detection system to access the
  field using all three values of ``'axis'``. This improvement avoids errors one
  would see now where only one value or an invalid value of the field parameter
  will be tested by yt. See `PR 1741
  <https://github.com/yt-project/yt/pull/1741>`__.
- It is now legal to pass a dataset instance as the first argument to
  ``ProfilePlot`` and ``PhasePlot``. This is equivalent to passing
  ``ds.all_data()``.
- Functions that accept a ``(length, unit)`` tuple (e.g. ``(3, 'km')`` for 3
  kilometers) will not raise an error if ``length`` is a ``YTQuantity`` instance
  with units attached. See `PR 1749
  <https://github.com/yt-project/yt/pull/1749>`__.
- The ``annotate_timestamp`` plot annotation now optionally accepts a
  ``time_offset`` keyword argument that sets the zero point of the time
  scale. Additionally, the ``annotate_scale`` plot annotation now accepts a
  ``format`` keyword argument, allowing custom formatting of the scale
  annotation. See `PR 1755 <https://github.com/yt-project/yt/pull/1755>`__.
- Add support for magnetic field variables and creation time fields in the GIZMO
  frontend. See `PR 1756 <https://github.com/yt-project/yt/pull/1756>`__ and `PR
  1914 <https://github.com/yt-project/yt/pull/1914>`__.
- ``ParticleProjectionPlot`` now supports the ``annotate_particles`` plot
  callback. See `PR 1765 <https://github.com/yt-project/yt/pull/1765>`__.
- Optmized the performance of off-axis projections for octree AMR data. See `PR
  1766 <https://github.com/yt-project/yt/pull/1766>`__.
- Added support for several radiative transfer fields in the ARTIO frontend. See
  `PR 1804 <https://github.com/yt-project/yt/pull/1804>`__.
- Performance improvement for Boxlib datasets that don't use AMR. See `PR 1834
  <https://github.com/yt-project/yt/pull/1834>`__.
- It is now possible to set custom profile bin edges. See `PR 1837
  <https://github.com/yt-project/yt/pull/1837>`__.
- Dropped support for Python3.4. See `PR 1840
  <https://github.com/yt-project/yt/pull/1840>`__.
- Add support for reading RAMSES cooling fields. See `PR 1853
  <https://github.com/yt-project/yt/pull/1853>`__.
- Add support for NumPy 1.15. See `PR 1854
  <https://github.com/yt-project/yt/pull/1854>`__.
- Ensure that functions defined in the plugins file are available in the yt
  namespace. See `PR 1855 <https://github.com/yt-project/yt/pull/1855>`__.
- Creating a profiles with log-scaled bins but where the bin edges are negative
  or zero now raises an error instead of silently generating a corrupt,
  incorrect answer. See `PR 1856
  <https://github.com/yt-project/yt/pull/1856>`__.
- Systematically added validation for inputs to data object initializers. See
  `PR 1871 <https://github.com/yt-project/yt/pull>`__.
- It is now possible to select only a specific particle type in the particle
  trajectories analysis module. See `PR 1887
  <https://github.com/yt-project/yt/pull/1887>`__.
- Substantially improve the performance of selecting particle fields with a
  ``cut_region`` data object. See `PR 1892
  <https://github.com/yt-project/yt/pull/1892>`__.
- The ``iyt`` command-line entry-point into IPython now installs yt-specific
  tab-completions. See `PR 1900 <https://github.com/yt-project/yt/pull/1900>`__.
- Derived quantities have been systematically updated to accept a
  ``particle_type`` keyword argument, allowing easier analysis of only a single
  particle type. See `PR 1902 <https://github.com/yt-project/yt/pull/1902>`__
  and `PR 1922 <https://github.com/yt-project/yt/pull/1922>`__.
- The ``annotate_streamlines()`` function now accepts a ``display_threshold``
  keyword argument. This suppresses drawing streamlines over any region of a
  dataset where the field being displayed is less than the threshold. See `PR
  1922 <https://github.com/yt-project/yt/pull/1922>`__.
- Add support for 2D nodal data. See `PR 1923
  <https://github.com/yt-project/yt/pull/1923>`__.
- Add support for GAMER outputs that use patch groups. This substantially
  reduces the memory requirements for loading large GAMER datasets. See `PR 1935
  <https://github.com/yt-project/yt/pull/1935>`__.
- Add a ``data_source`` keyword argument to the ``annotate_particles`` plot
  callback. See `PR 1937 <https://github.com/yt-project/yt/pull/1937>`__.
- Define species fields in the NMSU Art frontend. See `PR 1981
  <https://github.com/yt-project/yt/pull/1981>`__.
- Added a ``__format__`` implementation for ``YTArray``. See `PR 1985
  <https://github.com/yt-project/yt/pull/1985>`__.
- Derived fields that use a particle filter now only need to be derived for the
  particle filter type, not for the particle types used to define the particle
  filter. See `PR 1993 <https://github.com/yt-project/yt/pull/1993>`__.
- Added support for periodic visualizations using
  ``ParticleProjectionPlot``. See `PR 1996
  <https://github.com/yt-project/yt/pull/1996>`__.
- Added ``YTArray.argsort()``. See `PR 2002
  <https://github.com/yt-project/yt/pull/2002>`__.
- Calculate the header size from the header specification in the Gadget frontend
  to allow reading from Gadget binary datasets with nonstandard headers. See `PR
  2005 <https://github.com/yt-project/yt/pull/2005>`__ and `PR 2036
  <https://github.com/yt-project/yt/pull/2036>`__.
- Save the standard deviation in ``profile.save_as_dataset()``. See `PR 2008
  <https://github.com/yt-project/yt/pull/2008>`__.
- Allow the ``color`` keyword argument to be passed to matplotlib in the
  ``annotate_clumps`` callback to control the color of the clump annotation. See
  `PR 2019 <https://github.com/yt-project/yt/pull/2019>`__.
- Raise an exception when profiling fields of unequal shape. See `PR 2025
  <https://github.com/yt-project/yt/pull/2025>`__.
- The clump info dictionary is now populated as clumps get created instead of
  during ``clump.save_as_dataset()``. See `PR 2053
  <https://github.com/yt-project/yt/pull/2053>`__.
- Avoid segmentation fault in slice selector by clipping slice integer
  coordinates. See `PR 2055 <https://github.com/yt-project/yt/pull/2055>`__.

  
Minor Enhancements and Bugfixes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
- Fix incorrect use of floating point division in the parallel analysis framework.
  See `PR 1538 <https://github.com/yt-project/yt/pull/1538>`__.
- Fix integration with that matplotlib QT backend for interactive plotting.
  See `PR 1540 <https://github.com/yt-project/yt/pull/1540>`__.
- Add support for the particle creation time field in the GAMER frontend.
  See `PR 1546 <https://github.com/yt-project/yt/pull/1546>`__.
- Various minor improvements to the docs. See `PR 1542
  <https://github.com/yt-project/yt/pull/1542>`__. and `PR 1547
  <https://github.com/yt-project/yt/pull/1547>`__.
- Add better error handling for invalid tipsy aux files. See `PR 1549
  <https://github.com/yt-project/yt/pull/1549>`__.
- Fix typo in default Gadget header specification. See `PR 1550
  <https://github.com/yt-project/yt/pull/1550>`__.
- Use the git version in the get_yt_version function. See `PR 1551
  <https://github.com/yt-project/yt/pull/1551>`__.
- Assume dimensionless units for fields from FITS datasets when we can't infer
  the units. See `PR 1553 <https://github.com/yt-project/yt/pull/1553>`__.
- Autodetect ramses extra particle fields. See `PR 1555
  <https://github.com/yt-project/yt/pull/1555>`__.
- Fix issue with handling unitless halo quantities in HaloCatalog. See `PR 1558
  <https://github.com/yt-project/yt/pull/1558>`__.
- Track the halo catalog creation process using a parallel-safe progress bar. 
  See `PR 1559 <https://github.com/yt-project/yt/pull/1559>`__.
- The PPV Cube functionality no longer crashes if there is no temperature field
  in the dataset. See `PR 1562
  <https://github.com/yt-project/yt/pull/1562>`__.
- Fix crash caused by saving the ``'x'``, ``'y'``, or ``'z'`` fields in
  clump.save_as_dataset().  See `PR 1567
  <https://github.com/yt-project/yt/pull/1567>`__.
- Accept both string and tuple field names in ``ProfilePlot.set_unit()`` and
  ``PhasePlot.set_unit()``. See `PR 1568
  <https://github.com/yt-project/yt/pull/1568>`__.
- Fix issues with some arbitrary grid attributes not being reloaded properly
  after being saved with ``save_as_dataset()``. See `PR 1569
  <https://github.com/yt-project/yt/pull/1569>`__.
- Fix units issue in the light cone projection operation. See `PR 1574
  <https://github.com/yt-project/yt/pull/1574>`__.
- Use ``astropy.wcsaxes`` instead of the independent ``wcsaxes`` project.  See
  `PR 1577 <https://github.com/yt-project/yt/pull/1577>`__.
- Correct typo in WarpX field definitions. See `PR 1583
  <https://github.com/yt-project/yt/pull/1583>`__.
- Avoid crashing when loading an Enzo dataset with a parameter file that has
  commented out parameters. See `PR 1584
  <https://github.com/yt-project/yt/pull/1584>`__.
- Fix a corner case in the clump finding machinery where the reference to the
  parent clump is invalid after pruning a child clump that has no siblings. See
  `PR 1587 <https://github.com/yt-project/yt/pull/1587>`__.
- Fix issues with setting up yt fields for the magnetic and velocity field
  components and associated derived fields in curvilinear coordinate
  systems. See `PR 1588 <https://github.com/yt-project/yt/pull/1588>`__ and `PR
  1687 <https://github.com/yt-project/yt/pull/1687>`__.
- Fix incorrect profile values when the profile weight field has values equal to
  zero. See `PR 1590 <https://github.com/yt-project/yt/pull/1590>`__.
- Fix issues with making matplotlib animations of a
  ``ParticleProjectionPlot``. See `PR 1594
  <https://github.com/yt-project/yt/pull/1594>`__.
- The ``Scene.annotate_axes()`` function will now use the correct colors for
  drawing the axes annotation. See `PR 1596
  <https://github.com/yt-project/yt/pull/1596>`__.
- Fix incorrect default plot bounds for a zoomed-in slice plot of a 2D
  cylindrical dataset. See `PR 1597
  <https://github.com/yt-project/yt/pull/1597>`__.
- Fix issue where field accesses on 2D grids would return data with incorrect
  shapes. See `PR 1603 <https://github.com/yt-project/yt/pull/1603>`__.
- Added a cookbook example for a multipanel phase plot. See `PR 1605
  <https://github.com/yt-project/yt/pull/1605>`__.
- Boolean simulation parameters in the Boxlib frontend will now be interpreted
  correctly. See `PR 1618 <https://github.com/yt-project/yt/pull/1618>`__.
- The ``ds.particle_type_counts`` attribute will now be populated correctly for
  AMReX data.
- The ``"rad"`` unit (added for compatibility with astropy) now has the correct
  dimensions of angle instead of solid angle. See `PR 1628
  <https://github.com/yt-project/yt/pull/1628>`__.
- Fix units issues in several plot callbacks. See `PR 1633
  <https://github.com/yt-project/yt/pull/1633>`__ and `PR 1674
  <https://github.com/yt-project/yt/pull/1674>`__.
- Various fixes for how WarpX fields are interpreted. See `PR 1634
  <https://github.com/yt-project/yt/pull/1634>`__.
- Fix incorrect units in the automatically deposited particle fields. See `PR
  1638 <https://github.com/yt-project/yt/pull/1638>`__.
- It is now possible to set the axes background color after calling
  ``plot.hide_axes()``. See `PR 1662
  <https://github.com/yt-project/yt/pull/1662>`__.
- Fix a typo in the name of the ``colors`` keyword argument passed to matplotlib
  for the contour callback. See `PR 1664
  <https://github.com/yt-project/yt/pull/1664>`__.
- Add support for Enzo Active Particle fields that arrays. See `PR 1665
  <https://github.com/yt-project/yt/pull/1665>`__.
- Avoid crash when generating halo catalogs from the rockstar halo finder for
  small simulation domains. See `PR 1679
  <https://github.com/yt-project/yt/pull/1679>`__.
- The clump callback now functions correctly for a reloaded clump dataset. See
  `PR 1683 <https://github.com/yt-project/yt/pull/1683>`__.
- Fix incorrect calculation for tangential components of vector fields. See `PR
  1688 <https://github.com/yt-project/yt/pull/1688>`__.
- Allow halo finders to run in parallel on Python3. See `PR 1690
  <https://github.com/yt-project/yt/pull/1690>`__.
- Fix issues with Gadget particle IDs for simulations with large numbers of
  particles being incorrectly rounded. See `PR 1692
  <https://github.com/yt-project/yt/pull/1692>`__.
- ``ParticlePlot`` no longer needs to be passed spatial fields in a particular
  order to ensure that a ``ParticleProjectionPlot`` is returned. See `PR 1697
  <https://github.com/yt-project/yt/pull/1697>`__.
- Accessing data from a FLASH grid directly now returns float64 data. See `PR
  1708 <https://github.com/yt-project/yt/pull/1708>`__.
- Fix periodicity check in ``YTPoint`` data object. See `PR 1712
  <https://github.com/yt-project/yt/pull/1712>`__.
- Avoid crash on matplotlib 2.2.0 when generating yt plots with symlog
  colorbars. See `PR 1720 <https://github.com/yt-project/yt/pull/1720>`__.
- Avoid crash when FLASH ``"unitsystem"`` parameter is quoted in the HDF5
  file. See `PR 1722 <https://github.com/yt-project/yt/pull/1722>`__.
- Avoid issues with creating custom particle filters for OWLS/EAGLE
  datasets. See `PR 1723 <https://github.com/yt-project/yt/pull/1723>`__.
- Adapt to behavior change in matplotlib that caused plot inset boxes for
  annotated text to be drawn when none was requested. See `PR 1731
  <https://github.com/yt-project/yt/pull/1731>`__ and `PR 1827
  <https://github.com/yt-project/yt/pull/1827>`__.
- Fix clump finder ignoring field parameters. See `PR 1732
  <https://github.com/yt-project/yt/pull/1732>`__.
- Avoid generating NaNs in x-ray emission fields. See `PR 1742
  <https://github.com/yt-project/yt/pull/1742>`__.
- Fix compatibility with Sphinx 1.7 when building the docs. See `PR 1743
  <https://github.com/yt-project/yt/pull/1743>`__.
- Eliminate usage of deprecated ``"clobber"`` keyword argument for various
  usages of astropy in yt. See `PR 1744
  <https://github.com/yt-project/yt/pull/1744>`__.
- Fix incorrect definition of the ``"d"`` unit (an alias of ``"day"``). See `PR
  1746 <https://github.com/yt-project/yt/pull/1746>`__.
- ``PhasePlot.set_log()`` now correctly handles tuple field names as well as
  string field names. See `PR 1787
  <https://github.com/yt-project/yt/pull/1787>`__.
- Fix incorrect axis order in aitoff pixelizer. See `PR 1791
  <https://github.com/yt-project/yt/pull/1791>`__.
- Fix crash in when exporting a surface as a ply model. See `PR 1792
  <https://github.com/yt-project/yt/pull/1792>`__ and `PR 1817
  <https://github.com/yt-project/yt/pull/1817>`__.
- Fix crash in scene.save_annotated() in newer numpy versions. See `PR 1793
  <https://github.com/yt-project/yt/pull/1793>`__.
- Many tests no longer depend on real datasets. See `PR 1801
  <https://github.com/yt-project/yt/pull/1801>`__, `PR 1805
  <https://github.com/yt-project/yt/pull/1805>`__, `PR 1809
  <https://github.com/yt-project/yt/pull/1809>`__, `PR 1883
  <https://github.com/yt-project/yt/pull/1883>`__, and `PR 1941
  <https://github.com/yt-project/yt/pull/1941>`__
- New tests were added to improve test coverage or the performance of the
  tests. See `PR 1820 <https://github.com/yt-project/yt/pull/1820>`__, `PR 1831
  <https://github.com/yt-project/yt/pull/1831>`__, `PR 1833
  <https://github.com/yt-project/yt/pull/1833>`__, `PR 1841
  <https://github.com/yt-project/yt/pull/1841>`__, `PR 1842
  <https://github.com/yt-project/yt/pull/1842>`__, `PR 1885
  <https://github.com/yt-project/yt/pull/1885>`__, `PR 1886
  <https://github.com/yt-project/yt/pull/1886>`__, `PR 1952
  <https://github.com/yt-project/yt/pull/1952>`__, `PR 1953
  <https://github.com/yt-project/yt/pull/1953>`__, `PR 1955
  <https://github.com/yt-project/yt/pull/1955>`__, and `PR 1957
  <https://github.com/yt-project/yt/pull/1957>`__.
- The particle trajectories machinery will raise an error if it is asked to
  analyze a set of particles with duplicated particle IDs. See `PR 1818
  <https://github.com/yt-project/yt/pull/1818>`__.
- Fix incorrect velocity unit int he ``gadget_fof`` frontend. See `PR 1829
  <https://github.com/yt-project/yt/pull/1829>`__.
- Making an off-axis projection of a cut_region data object with an octree AMR
  dataset now works correctly. See `PR 1858
  <https://github.com/yt-project/yt/pull/1858>`__.
- Replace hard-coded constants in Enzo frontend with calculations to improve
  agreement with Enzo's internal constants and improve clarity. See `PR 1873
  <https://github.com/yt-project/yt/pull/1873>`__.
- Correct issues with Enzo magnetic units in cosmology simulations. See `PR 1876
  <https://github.com/yt-project/yt/pull/1876>`__.
- Use the species names from the dataset rather than hardcoding species names in
  the WarpX frontend. See `PR 1884
  <https://github.com/yt-project/yt/pull/1884>`__.
- Fix issue with masked I/O for unstructured mesh data. See `PR 1918
  <https://github.com/yt-project/yt/pull/1918>`__.
- Fix crash when reading DM-only Enzo datasets where some grids have no particles. See `PR 1919 <https://github.com/yt-project/yt/pull/1919>`__.
- Fix crash when loading pure-hydro Nyx dataset. See `PR 1950
  <https://github.com/yt-project/yt/pull/1950>`__.
- Avoid crashes when plotting fields that contain NaN. See `PR 1951
  <https://github.com/yt-project/yt/pull/1951>`__.
- Avoid crashes when loading NMSU ART data. See `PR 1960
  <https://github.com/yt-project/yt/pull/1960>`__.
- Avoid crash when loading WarpX dataset with no particles. See `PR 1979
  <https://github.com/yt-project/yt/pull/1979>`__.
- Adapt to API change in glue to fix the ``to_glue()`` method on yt data
  objects. See `PR 1991 <https://github.com/yt-project/yt/pull/1991>`__.
- Fix incorrect width calculation in the ``annotate_halos()`` plot callback. See
  `PR 1995 <https://github.com/yt-project/yt/pull/1995>`__.
- Don't try to read from files containing zero halos in the ``gadget_fof``
  frontend. See `PR 2001 <https://github.com/yt-project/yt/pull/2001>`__.
- Fix incorrect calculation in ``get_ortho_base()``. See `PR 2013
  <https://github.com/yt-project/yt/pull/2013>`__.
- Avoid issues with the axes background color being inconsistently set. See `PR
  2018 <https://github.com/yt-project/yt/pull/2018>`__.
- Fix issue with reading multiple fields at once for octree AMR data sometimes
  returning data for another field for one of the requested fields. See `PR 2020
  <https://github.com/yt-project/yt/pull/2020>`__.
- Fix incorrect domain annotation for ``Scene.annotate_domain()`` when using the
  plane-parallel camera. See `PR 2024
  <https://github.com/yt-project/yt/pull/2024>`__.
- Avoid crash when particles are on the domain edges for ``gadget_fof``
  data. See `PR 2034 <https://github.com/yt-project/yt/pull/2034>`__.
- Avoid stripping code units when processing units through a dataset's unit
  system. See `PR 2035 <https://github.com/yt-project/yt/pull/2035>`__.
- Avoid incorrectly rescaling units of metalicity fields. See `PR 2038
  <https://github.com/yt-project/yt/pull/2038>`__.
- Fix incorrect units for FLASH ``"divb"`` field. See `PR 2062
  <https://github.com/yt-project/yt/pull/2062>`__.

Version 3.4
-----------

Version 3.4 is the first major release of yt since July 2016. It includes 450
pull requests from 44 contributors including 18 new contributors.

-  yt now supports displaying plots using the interactive matplotlib
   backends. To enable this functionality call
   ``yt.toggle_interactivity()``. This is currently supported at an
   experimental level, please let us know if you come across issues
   using it. See `Bitbucket PR
   2294 <https://bitbucket.org/yt_analysis/yt/pull-requests/2294>`__.
-  The yt configuration file should now be located in a location
   following the XDG\_CONFIG convention (usually ``~/.config/yt/ytrc``)
   rather than the old default location (usually ``~/.yt/config``). You
   can use ``yt config migrate`` at the bash command line to migrate
   your configuration file to the new location. See `Bitbucket PR
   2343 <https://bitbucket.org/yt_analysis/yt/pull-requests/2343>`__.
-  Added ``yt.LinePlot``, a new plotting class for creating 1D plots
   along lines through a dataset. See `Github PR
   1509 <https://github.com/yt-project/yt/pull/1509>`__ and `Github PR
   1440 <https://github.com/yt-project/yt/pull/1440>`__.
-  Added ``yt.define_unit`` to easily define new units in yt's unit
   system. See `Bitbucket PR
   2485 <https://bitbucket.org/yt_analysis/yt/pull-requests/2485>`__.
-  Added ``yt.plot_2d``, a wrapper around SlicePlot for plotting 2D
   datasets. See `Github PR
   1476 <https://github.com/yt-project/yt/pull/1476>`__.
-  We have restored support for boolean data objects. Boolean objects
   are data objects that are defined in terms of boolean operations on
   other data objects. See `Bitbucket PR
   2257 <https://bitbucket.org/yt_analysis/yt/pull-requests/2257>`__.
-  Datasets now have a ``fields`` attribute that allows access to fields
   via a python object. For example, instead of using a tuple field name
   like ``('gas', 'density')``, one can now use
   ``ds.fields.gas.density``. See `Bitbucket PR
   2459 <https://bitbucket.org/yt_analysis/yt/pull-requests/2459>`__.
-  It is now possible to create a wider variety of data objects via
   ``ds.r``, including rays, fixed resolution rays, points, and images.
   See `Github PR 1518 <https://github.com/yt-project/yt/pull/1518>`__
   and `Github PR 1393 <https://github.com/yt-project/yt/pull/1393>`__.
-  ``add_field`` and ``ds.add_field`` must now be called with a
   ``sampling_type`` keyword argument. Possible values are currently
   ``cell`` and ``particle``. We have also deprecated the
   ``particle_type`` keyword argument in favor of
   ``sampling_type='cell'``. For now a ``'cell'`` ``sampling_type`` is
   assumed if ``sampling_type`` is not specified but in the future
   ``sampling_type`` will always need to be specified.
-  Added support for the ``Athena++`` code. See `Bitbucket PR
   2149 <https://bitbucket.org/yt_analysis/yt/pull-requests/2149>`__.
-  Added support for the ``Enzo-p`` code. See `Github PR
   1447 <https://github.com/yt-project/yt/pull/1447>`__, `Github PR
   1443 <https://github.com/yt-project/yt/pull/1443>`__ and `Github PR
   1439 <https://github.com/yt-project/yt/pull/1439>`__.
-  Added support for the ``AMReX`` code. See `Bitbucket PR
   2530 <https://bitbucket.org/yt_analysis/yt/pull-requests/2530>`__.
-  Added support for the ``openPMD`` output format. See `Bitbucket PR
   2376 <https://bitbucket.org/yt_analysis/yt/pull-requests/2376>`__.
-  Added support for reading face-centered and vertex-centered fields
   for block AMR codes. See `Bitbucket PR
   2575 <https://bitbucket.org/yt_analysis/yt/pull-requests/2575>`__.
-  Added support for loading outputs from the Amiga Halo Finder. See
   `Github PR 1477 <https://github.com/yt-project/yt/pull/1477>`__.
-  Added support for particle fields for Boxlib data. See `Bitbucket PR
   2510 <https://bitbucket.org/yt_analysis/yt/pull-requests/2510>`__ and
   `Bitbucket PR
   2497 <https://bitbucket.org/yt_analysis/yt/pull-requests/2497>`__.
-  Added support for custom RAMSES particle fields. See `Github PR
   1470 <https://github.com/yt-project/yt/pull/1470>`__.
-  Added support for RAMSES-RT data. See `Github PR
   1456 <https://github.com/yt-project/yt/pull/1456>`__ and `Github PR
   1449 <https://github.com/yt-project/yt/pull/1449>`__.
-  Added support for Enzo MHDCT fields. See `Github PR
   1438 <https://github.com/yt-project/yt/pull/1438>`__.
-  Added support for units and particle fields to the GAMER frontend.
   See `Bitbucket PR
   2366 <https://bitbucket.org/yt_analysis/yt/pull-requests/2366>`__ and
   `Bitbucket PR
   2408 <https://bitbucket.org/yt_analysis/yt/pull-requests/2408>`__.
-  Added support for type 2 Gadget binary outputs. See `Bitbucket PR
   2355 <https://bitbucket.org/yt_analysis/yt/pull-requests/2355>`__.
-  Added the ability to detect and read double precision Gadget data.
   See `Bitbucket PR
   2537 <https://bitbucket.org/yt_analysis/yt/pull-requests/2537>`__.
-  Added the ability to detect and read in big endian Gadget data. See
   `Github PR 1353 <https://github.com/yt-project/yt/pull/1353>`__.
-  Added support for Nyx datasets that do not contain particles. See
   `Bitbucket PR
   2571 <https://bitbucket.org/yt_analysis/yt/pull-requests/2571>`__
-  A number of untested and unmaintained modules have been deprecated
   and moved to the `yt attic
   repository <https://github.com/yt-project/yt_attic>`__. This includes
   the functionality for calculating two point functions, the Sunrise
   exporter, the star analysis module, and the functionality for
   calculating halo mass functions. If you are interested in working on
   restoring the functionality in these modules, we welcome
   contributions. Please contact us on the mailing list or by opening an
   issue on GitHub if you have questions.
-  The particle trajectories functionality has been removed from the
   analysis modules API and added as a method of the ``DatasetSeries``
   object. You can now create a ``ParticleTrajectories`` object using
   ``ts.particle_trajectories()`` where ``ts`` is a time series of
   datasets.
-  The ``spectral_integrator`` analysis module is now available via
   ``yt.fields.xray_emission_fields``. See `Bitbucket PR
   2465 <https://bitbucket.org/yt_analysis/yt/pull-requests/2465>`__.
-  The ``photon_simulator`` analysis module has been deprecated in favor
   of the ``pyXSIM`` package, available separately from ``yt``. See
   `Bitbucket PR
   2441 <https://bitbucket.org/yt_analysis/yt/pull-requests/2441>`__.
-  ``yt.utilities.fits_image`` is now available as
   ``yt.visualization.fits_image``. In addition classes that were in the
   ``yt.utilities.fits_image`` namespace are now available in the main
   ``yt`` namespace.
-  The ``profile.variance`` attribute has been deprecated in favor of
   ``profile.standard_deviation``.
-  The ``number_of_particles`` key no longer needs to be defined when
   loading data via the stream frontend. See `Github PR
   1428 <https://github.com/yt-project/yt/pull/1428>`__.
-  The install script now only supports installing via miniconda. We
   have removed support for compiling python and yt's dependencies from
   source. See `Github PR
   1459 <https://github.com/yt-project/yt/pull/1459>`__.
-  Added ``plot.set_background_color`` for ``PlotWindow`` and
   ``PhasePlot`` plots. This lets users specify a color to fill in the
   background of a plot instead of the default color, white. See
   `Bitbucket PR
   2513 <https://bitbucket.org/yt_analysis/yt/pull-requests/2513>`__.
-  ``PlotWindow`` plots can now optionally use a right-handed coordinate
   system. See `Bitbucket PR
   2318 <https://bitbucket.org/yt_analysis/yt/pull-requests/2318>`__.
-  The isocontour API has been overhauled to make use of units. See
   `Bitbucket PR
   2453 <https://bitbucket.org/yt_analysis/yt/pull-requests/2453>`__.
-  ``Dataset`` instances now have a ``checksum`` property, which can be
   accessed via ``ds.checksum``. This provides a unique identifier that
   is guaranteed to be the same from session to session. See `Bitbucket
   PR 2503 <https://bitbucket.org/yt_analysis/yt/pull-requests/2503>`__.
-  Added a ``data_source`` keyword argument to
   ``OffAxisProjectionPlot``. See `Bitbucket PR
   2490 <https://bitbucket.org/yt_analysis/yt/pull-requests/2490>`__.
-  Added a ``yt download`` command-line helper to download test data
   from http://yt-project.org/data. For more information see
   ``yt download --help`` at the bash command line. See `Bitbucket PR
   2495 <https://bitbucket.org/yt_analysis/yt/pull-requests/2495>`__ and
   `Bitbucket PR
   2471 <https://bitbucket.org/yt_analysis/yt/pull-requests/2471>`__.
-  Added a ``yt upload`` command-line helper to upload files to the `yt
   curldrop <https://docs.hub.yt/services.html#curldrop>`__ at the bash
   command line. See `Github PR
   1471 <https://github.com/yt-project/yt/pull/1471>`__.
-  If it's installed, colormaps from the `cmocean
   package <http://matplotlib.org/cmocean/>`__ will be made available as
   yt colormaps. See `Bitbucket PR
   2439 <https://bitbucket.org/yt_analysis/yt/pull-requests/2439>`__.
-  It is now possible to visualize unstructured mesh fields defined on
   multiple mesh blocks. See `Bitbucket PR
   2487 <https://bitbucket.org/yt_analysis/yt/pull-requests/2487>`__.
-  Add support for second-order interpolation when slicing tetrahedral
   unstructured meshes. See `Bitbucket PR
   2550 <https://bitbucket.org/yt_analysis/yt/pull-requests/2550>`__.
-  Add support for volume rendering second-order tetrahedral meshes. See
   `Bitbucket PR
   2401 <https://bitbucket.org/yt_analysis/yt/pull-requests/2401>`__.
-  Add support for QUAD9 mesh elements. See `Bitbucket PR
   2549 <https://bitbucket.org/yt_analysis/yt/pull-requests/2549>`__.
-  Add support for second-order triangle mesh elements. See `Bitbucket
   PR 2378 <https://bitbucket.org/yt_analysis/yt/pull-requests/2378>`__.
-  Added support for dynamical dark energy parameterizations to the
   ``Cosmology`` object. See `Bitbucket PR
   2572 <https://bitbucket.org/yt_analysis/yt/pull-requests/2572>`__.
-  ``ParticleProfile`` can now handle log-scaled bins and data with
   negative values. See `Bitbucket PR
   2564 <https://bitbucket.org/yt_analysis/yt/pull-requests/2564>`__ and
   `Github PR 1510 <https://github.com/yt-project/yt/pull/1510>`__.
-  Cut region data objects can now be saved as reloadable datasets using
   ``save_as_dataset``. See `Bitbucket PR
   2541 <https://bitbucket.org/yt_analysis/yt/pull-requests/2541>`__.
-  Clump objects can now be saved as reloadable datasets using
   ``save_as_dataset``. See `Bitbucket PR
   2326 <https://bitbucket.org/yt_analysis/yt/pull-requests/2326>`__.
-  It is now possible to specify the field to use for the size of the
   circles in the ``annotate_halos`` plot modifying function. See
   `Bitbucket PR
   2493 <https://bitbucket.org/yt_analysis/yt/pull-requests/2493>`__.
-  The ``ds.max_level`` attribute is now a property that is computed on
   demand. The more verbose ``ds.index.max_level`` will continue to
   work. See `Bitbucket PR
   2461 <https://bitbucket.org/yt_analysis/yt/pull-requests/2461>`__.
-  The ``PointSource`` volume rendering source now optionally accepts a
   ``radius`` keyword argument to draw spatially extended points. See
   `Bitbucket PR
   2404 <https://bitbucket.org/yt_analysis/yt/pull-requests/2404>`__.
-  It is now possible to save volume rendering images in eps, ps, and
   pdf format. See `Github PR
   1504 <https://github.com/yt-project/yt/pull/1504>`__.

Minor Enhancements and Bugfixes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Fixed issue selecting and visualizing data at very high AMR levels.
   See `Github PR 1521 <https://github.com/yt-project/yt/pulls/1521>`__
   and `Github PR 1433 <https://github.com/yt-project/yt/pull/1433>`__.
-  Print a more descriptive error message when defining a particle
   filter fails with missing fields See `Github PR
   1517 <https://github.com/yt-project/yt/pull/1517>`__.
-  Removed grid edge rounding from the FLASH frontend. This fixes a
   number of pernicious visualization artifacts for FLASH data. See
   `Github PR 1493 <https://github.com/yt-project/yt/pull/1493>`__.
-  Parallel projections no longer error if there are less io chunks than
   MPI tasks. See `Github PR
   1488 <https://github.com/yt-project/yt/pull/1488>`__.
-  A memory leak in the volume renderer has been fixed. See `Github PR
   1485 <https://github.com/yt-project/yt/pull/1485>`__ and `Github PR
   1435 <https://github.com/yt-project/yt/pull/1435>`__.
-  The ``force_override`` keyword argument now raises an error when used
   with on-disk fields. See `Github PR
   1516 <https://github.com/yt-project/yt/pull/1516>`__.
-  Restore support for making plots from reloaded plots. See `Github PR
   1514 <https://github.com/yt-project/yt/pull/1514>`__
-  Don't ever try to read inputs or probin files for Castro and Maestro.
   See `Github PR 1445 <https://github.com/yt-project/yt/pull/1445>`__.
-  Fixed issue that caused visualization artifacts when creating an
   off-axis projection for particle or octree AMR data. See `Github PR
   1434 <https://github.com/yt-project/yt/pull/1434>`__.
-  Fix i/o for the Enzo ``'Dark_Matter_Density'`` field. See `Github PR
   1360 <https://github.com/yt-project/yt/pull/1360>`__.
-  Create the ``'particle_ones'`` field even if we don't have a particle
   mass field. See `Github PR
   1424 <https://github.com/yt-project/yt/pull/1424>`__.
-  Fixed issues with minor colorbar ticks with symlog colorbar scaling.
   See `Github PR 1423 <https://github.com/yt-project/yt/pull/1423>`__.
-  Using the rockstar halo finder is now supported under Python3. See
   `Github PR 1414 <https://github.com/yt-project/yt/pull/1414>`__.
-  Fixed issues with orientations of volume renderings when compositing
   multiple sources. See `Github PR
   1411 <https://github.com/yt-project/yt/pull/1411>`__.
-  Added a check for valid AMR structure in ``load_amr_grids``. See
   `Github PR 1408 <https://github.com/yt-project/yt/pull/1408>`__.
-  Fix bug in handling of periodic boundary conditions in the
   ``annotate_halos`` plot modifying function. See `Github PR
   1351 <https://github.com/yt-project/yt/pull/1351>`__.
-  Add support for plots with non-unit aspect ratios to the
   ``annotate_scale`` plot modifying function. See `Bitbucket PR
   2551 <https://bitbucket.org/yt_analysis/yt/pull-requests/2551>`__.
-  Fixed issue with saving light ray datasets. See `Bitbucket PR
   2589 <https://bitbucket.org/yt_analysis/yt/pull-requests/2589>`__.
-  Added support for 2D WarpX data. ee `Bitbucket PR
   2583 <https://bitbucket.org/yt_analysis/yt/pull-requests/2583>`__.
-  Ensure the ``particle_radius`` field is always accessed with the
   correct field type. See `Bitbucket PR
   2562 <https://bitbucket.org/yt_analysis/yt/pull-requests/2562>`__.
-  It is now possible to use a covering grid to access particle filter
   fields. See `Bitbucket PR
   2569 <https://bitbucket.org/yt_analysis/yt/pull-requests/2569>`__.
-  The x limits of a ``ProfilePlot`` will now snap exactly to the limits
   specified in calls to ``ProfilePlot.set_xlim``. See `Bitbucket PR
   2546 <https://bitbucket.org/yt_analysis/yt/pull-requests/2546>`__.
-  Added a cookbook example showing how to make movies using
   matplotlib's animation framework. See `Bitbucket PR
   2544 <https://bitbucket.org/yt_analysis/yt/pull-requests/2544>`__.
-  Use a parallel-safe wrapper around mkdir when creating new
   directories. See `Bitbucket PR
   2570 <https://bitbucket.org/yt_analysis/yt/pull-requests/2570>`__.
-  Removed ``yt.utilities.spatial``. This was a forked version of
   ``scipy.spatial`` with support for a periodic KD-tree. Scipy now has
   a periodic KD-tree, so we have removed the forked version from yt.
   Please use ``scipy.spatial`` if you were relying on
   ``yt.utilities.spatial``. See `Bitbucket PR
   2576 <https://bitbucket.org/yt_analysis/yt/pull-requests/2576>`__.
-  Improvements for the ``HaloCatalog``. See `Bitbucket PR
   2536 <https://bitbucket.org/yt_analysis/yt/pull-requests/2536>`__ and
   `Bitbucket PR
   2535 <https://bitbucket.org/yt_analysis/yt/pull-requests/2535>`__.
-  Removed ``'log'`` in colorbar label in annotated volume rendering.
   See `Bitbucket PR
   2548 <https://bitbucket.org/yt_analysis/yt/pull-requests/2548>`__
-  Fixed a crash triggered by depositing particle data onto a covering
   grid. See `Bitbucket PR
   2545 <https://bitbucket.org/yt_analysis/yt/pull-requests/2545>`__.
-  Ensure field type guessing is deterministic on Python3. See
   `Bitbucket PR
   2559 <https://bitbucket.org/yt_analysis/yt/pull-requests/2559>`__.
-  Removed unused yt.utilities.exodusII\_reader module. See `Bitbucket
   PR 2533 <https://bitbucket.org/yt_analysis/yt/pull-requests/2533>`__.
-  The ``cell_volume`` field in curvilinear coordinates now uses an
   exact rather than an approximate definition. See `Bitbucket PR
   2466 <https://bitbucket.org/yt_analysis/yt/pull-requests/2466>`__.

Version 3.3
-----------

Version 3.3 is the first major release of yt since July 2015. It includes more
than 3000 commits from 41 contributors, including 12 new contributors.

Major enhancements
^^^^^^^^^^^^^^^^^^

* Raw and processed data from selections, projections, profiles and so forth can
  now be saved in a ytdata format and loaded back in by yt. See 
  :ref:`saving_data`.
* Totally re-worked volume rendering API. The old API is still available for users
  who prefer it, however. See :ref:`volume_rendering`.
* Support for unstructured mesh visualization. See 
  :ref:`unstructured-mesh-slices` and :ref:`unstructured_mesh_rendering`.
* Interactive Data Visualization for AMR and unstructured mesh datasets. See
  :ref:`interactive_data_visualization`.
* Several new colormaps, including a new default, 'arbre'. The other new
  colormaps are named 'octarine', 'kelp', and 'dusk'. All these new colormaps
  were generated using the `viscm package
  <https://github.com/matplotlib/viscm>`_ and should do a better job of
  representing the data for colorblind viewers and when printed out in
  grayscale. See :ref:`colormaps` for more detail.
* New frontends for the :ref:`ExodusII <loading-exodusii-data>`, 
  :ref:`GAMER <loading-gamer-data>`, and :ref:`Gizmo <loading-gizmo-data>` data 
  formats.
* The unit system associated with a dataset is now customizable, defaulting to
  CGS. See :ref:`unit_systems`.
* Enhancements and usability improvements for analysis modules, especially the
  ``absorption_spectrum``, ``photon_simulator``, and ``light_ray`` modules. See
  :ref:`synthetic-observations`.
* Data objects can now be created via an alternative Numpy-like API. See
  :ref:`quickly-selecting-data`.
* A line integral convolution plot modification. See
  :ref:`annotate-line-integral-convolution`.
* Many speed optimizations, including to the volume rendering, units, tests,
  covering grids, the absorption spectrum and photon simulator analysis modules,
  and ghost zone generation.
* Packaging and release-related improvements: better install and setup scripts,
  automated PR backporting.
* Readability improvements to the codebase, including linting, removing dead
  code, and refactoring much of the Cython.
* Improvements to the CI infrastructure, including more extensible answer tests
  and automated testing for Python 3 and Windows.
* Numerous documentation improvements, including formatting tweaks, bugfixes,
  and many new cookbook recipes.
* Support for geographic (lat/lon) coordinates.
* Several improvements for SPH codes, including alternative smoothing kernels,
  an ``add_smoothed_particle_field`` function, and particle type-aware octree
  construction for Gadget data.
* Roundtrip conversions between Pint and yt units.
* Added halo data containers for gadget_fof frontend.
* Enabled support for spherical datasets in the BoxLib frontend.
* Many new tests have been added.
* Better hashing for Selector objects.

Minor enhancements and bugfixes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Fixed many bugs related to Python 3 compatibility
* Fixed bugs related to compatibility issues with newer versions of numpy
* Added the ability to export data objects to a Pandas dataframe
* Added support for the fabs ufunc to YTArray
* Fixed two licensing issues
* Fixed a number of bugs related to Windows compatibility.
* We now avoid hard-to-decipher tracebacks when loading empty files or
  directories
* Fixed a bug related to ART star particle creation time field
* Fixed a bug caused by using the wrong int type for indexing in particle deposit
* Fixed a NameError bug in comparing temperature units with offsets
* Fixed an API bug in YTArray casting during coercion from YTQuantity
* Added loadtxt and savetxt convenience functions for ``YTArray``
* Fixed an issue caused by not sort species names with Enzo
* Fixed a units bug for RAMSES when ``boxlen > 1``.
* Fixed ``process_chunk`` function for non-cartesian geometry.
* Added ``scale_factor`` attribute to cosmological simulation datasets
* Fixed a bug where "center" vectors are used instead of "normal" vectors in
  get_sph_phi(), etc.
* Fixed issues involving invalid FRBs when uses called _setup_plots in their
  scripts
* Added a ``text_args`` keyword to ``annotate_scale()`` callback
* Added a print_stats function for RAMSES
* Fixed a number of bugs in the Photon Simulator
* Added support for particle fields to the [Min,Max]Location derived quantities
* Fixed some units bugs for Gadget cosmology simulations
* Fixed a bug with Gadget/GIZMO StarFormationRate units
* Fixed an issue in TimeSeriesData where all the filenames were getting passed
  to ``load`` on each processor.
* Fixed a units bug in the Tipsy frontend
* Ensured that ARTIOIndex.get_smallest_dx() returns a quantity with units
* Ensured that plots are valid after invalidating the figure
* Fixed a bug regarding code unit labels
* Fixed a bug with reading Tipsy Aux files
* Added an effective redshift field to the Light Ray analysis module for use in
  AbsorptionSpectrum
* Fixed a bug with the redshift calculation in LightRay analysis module
* Fixed a bug in the Orion frontend when you had more than 10 on-disk particle
  fields in the file
* Detect more types of ART files
* Update derived_field_list in add_volume_weighted_smoothed_field
* Fixed casting issues for 1D and 2D Enzo simulations
* Avoid type indirection when setting up data object entry points
* Fixed issues with SIMPUT files
* Fixed loading athena data in python3 with provided parameters
* Tipsy cosmology unit fixes
* Fixed bad unit labels for compound units
* Making the xlim and ylim of the PhasePlot plot axes controllable
* Adding grid_arrays to grid_container
* An Athena and a GDF bugfix
* A small bugfix and some small enhancements for sunyaev_zeldovich
* Defer to coordinate handlers for width
* Make array_like_field return same units as get_data
* Fixing bug in ray "dts" and "t" fields
* Check against string_types not str
* Closed a loophole that allowed improper LightRay use
* Enabling AbsorptionSpectrum to deposit unresolved spectral lines
* Fixed an ART byte/string/array issue
* Changing AbsorptionSpectrum attribute lambda_bins to be lambda_field for
  consistency
* No longer require user to save to disk when generating an AbsorptionSpectrum
* ParticlePlot FRBs can now use save_as_dataset and save attributes properly
* Added checks to assure ARTIO creates a metal_density field from existing metal
  fields.
* Added mask to LightRay to assure output elements have non-zero density (a
  problem in some SPH datasets)
* Added a "fields" attribute to datasets
* Updated the TransferFunctionHelper to work with new profiles
* Fixed a bug where the field_units kwarg to load_amr_grids didn't do anything
* Changed photon_simulator's output file structure
* Fixed a bug related to setting output_units.
* Implemented ptp operation.
* Added effects of transverse doppler redshift to LightRay
* Fixed a casting error for float and int64 multiplication in sdf class
* Added ability to read and write YTArrays to and from groups within HDF5 files
* Made ftype of "on-disk" stream fields "stream"
* Fixed a strings decoding issue in the photon simulator
* Fixed an incorrect docstring in load_uniform_grid
* Made PlotWindow show/hide helpers for axes and colorbar return self
* Made Profile objects store field metadata.
* Ensured GDF unit names are strings
* Taught off_axis_projection about its resolution keyword.
* Reintroduced sanitize_width for polar/cyl coordinates.
* We now fail early when load_uniform_grid is passed data with an incorrect shape
* Replaced progress bar with tqdm
* Fixed redshift scaling of "Overdensity" field in yt-2.x
* Fixed several bugs in the eps_writer
* Fixed bug affecting 2D BoxLib simulations.
* Implemented to_json and from_json for the UnitRegistry object
* Fixed a number of issues with ds.find_field_values_at_point[s]
* Fixed a bug where sunrise_exporter was using wrong imports
* Import HUGE from utilities.physical_ratios
* Fixed bug in ARTIO table look ups
* Adding support for longitude and latitude
* Adding halo data containers for gadget_fof frontend.
* Can now compare YTArrays without copying them
* Fixed several bugs related to active particle datasets
* Angular_momentum_vector now only includes space for particle fields if they
  exist.
* Image comparison tests now print a meaningful error message if they fail.
* Fixed numpy 1.11 compatibility issues.
* Changed _skip_cache to be True by default.
* Enable support for spherical datasets in the BoxLib frontend.
* Fixed a bug in add_deposited_particle_field.
* Fixed issues with input sanitization in the point data object.
* Fixed a copy/paste error introduced by refactoring WeightedMenParticleField
* Fixed many formatting issues in the docs build
* Now avoid creating particle unions for particle types that have no common
  fields
* Patched ParticlePlot to work with filtered particle fields.
* Fixed a couple corner cases in gadget_fof frontend
* We now properly normalise all normal vectors in functions that take a normal
  vector (for e.g. get_sph_theta)
* Fixed a bug where the transfer function features were not always getting
  cleared properly.
* Made the Chombo frontend is_valid method smarter.
* Added a get_hash() function to yt/funcs.py which returns a hash for a file
* Added Sievert to the default unit symbol table
* Corrected an issue with periodic "wiggle" in AbsorptionSpectrum instances
* Made ``ds.field_list`` sorted by default
* Bug fixes for the Nyx frontend
* Fixed a bug where the index needed to be created before calling derived
  quantities
* Made latex_repr a property, computed on-demand
* Fixed a bug in off-axis slice deposition
* Fixed a bug with some types of octree block traversal
* Ensured that mpi operations retain ImageArray type instead of downgrading to
  YTArray parent class
* Added a call to _setup_plots in the custom colorbar tickmark example
* Fixed two minor bugs in save_annotated
* Added ability to specify that DatasetSeries is not a mixed data type
* Fixed a memory leak in ARTIO
* Fixed copy/paste error in to_frb method.
* Ensured that particle dataset max_level is consistent with the index max_level
* Fixed an issue where fields were getting added multiple times to
  field_info.field_list
* Enhanced annotate_ray and annotate_arrow callbacks
* Added GDF answer tests
* Made the YTFieldTypeNotFound exception more informative
* Added a new function, fake_vr_orientation_test_ds(), for use in testing
* Ensured that instances of subclasses of YTArray have the correct type
* Re-enabled max_level for projections, ProjectionPlot, and OffAxisProjectionPlot
* Fixed a bug in the Orion 2 field definitions
* Fixed a bug caused by matplotlib not being added to install_requires
* Edited PhasePlot class to have an annotate_title method
* Implemented annotate_cell_edges
* Handled KeyboardInterrupt in volume rendering Cython loop
* Made old halo finders now accept ptype
* Updated the latex commands in yt cheatsheet
* Fixed a circular dependency loop bug in abar field definition for FLASH
  datasets
* Added neutral species aliases as described in YTEP 0003
* Fixed a logging issue: don't create a StreamHandler unless we will use it
* Correcting how theta and phi are calculated in
  ``_particle_velocity_spherical_radius``,
  ``_particle_velocity_spherical_theta``,
  ``_particle_velocity_cylindrical_radius``, and
  ``_particle_velocity_cylindrical_theta``
* Fixed a bug related to the field dictionary in ``load_particles``
* Allowed for the special case of supplying width as a tuple of tuples
* Made yt compile with MSVC on Windows
* Fixed a bug involving mask for dt in octree
* Merged the get_yt.sh and install_script.sh into one
* Added tests for the install script
* Allowed use axis names instead of dimensions for spherical pixelization
* Fixed a bug where close() wasn't being called in HDF5FileHandler
* Enhanced commandline image upload/delete
* Added get_brewer_cmap to get brewer colormaps without importing palettable at
  the top level
* Fixed a bug where a parallel_root_only function was getting called inside
  another parallel_root_only function
* Exit the install script early if python can't import '_ssl' module
* Make PlotWindow's annotate_clear method invalidate the plot
* Adding int wrapper to avoid deprecation warning from numpy
* Automatically create vector fields for magnetic_field
* Allow users to completely specify the filename of a 1D profile
* Force nose to produce meaningful traceback for cookbook recipes' tests
* Fixed x-ray display_name and documentation
* Try to guess and load particle file for FLASH dataset
* Sped up top-level yt import
* Set the field type correctly for fields added as particle fields
* Added a position location method for octrees
* Fixed a copy/paste error in uhstack function
* Made trig functions give correct results when supplied data with dimensions of
  angle but units that aren't radian
* Print out some useful diagnostic information if check_for_openmp() fails
* Give user-added derived fields a default field type
* Added support for periodicity in annotate_particles.
* Added a check for whether returned field has units in volume-weighted smoothed
  fields
* Casting array indices as ints in colormaps infrastructure
* Fixed a bug where the standard particle fields weren't getting set up
  correctly for the Orion frontends
* Enabled LightRay to accept loaded datasets instead of just filenames
* Allowed for adding or subtracting arrays filled with zeros without checking
  units.
* Fixed a bug in selection for semistructured meshes.
* Removed 'io' from enzo particle types for active particle datasets
* Added support for FLASH particle datasets.
* Silenced a deprecation warning from IPython
* Eliminated segfaults in KDTree construction
* Fixed add_field handling when passed a tuple
* Ensure field parameters are correct for fields that need ghost zones
* Made it possible to use DerivedField instances to access data
* Added ds.particle_type_counts
* Bug fix and improvement for generating Google Cardboard VR in
  StereoSphericalLens
* Made DarkMatterARTDataset more robust in its _is_valid
* Added Earth radius to units
* Deposit hydrogen fields to grid in gizmo frontend
* Switch to index values being int64
* ValidateParameter ensures parameter values are used during field detection
* Switched to using cythonize to manage dependencies in the setup script
* ProfilePlot style changes and refactoring
* Cancel terms with identical LaTeX representations in a LaTeX representation of
  a unit
* Only return early from comparison validation if base values are equal
* Enabled particle fields for clump objects
* Added validation checks for data types in callbacks
* Enabled modification of image axis names in coordinate handlers
* Only add OWLS/EAGLE ion fields if they are present
* Ensured that PlotWindow plots continue to look the same under matplotlib 2.0
* Fixed bug in quiver callbacks for off-axis slice plots
* Only visit octree children if going to next level
* Check that CIC always gets at least two cells
* Fixed compatibility with matplotlib 1.4.3 and earlier
* Fixed two EnzoSimulation bugs
* Moved extraction code from YTSearchCmd to its own utility module
* Changed amr_kdtree functions to be Node class methods
* Sort block indices in order of ascending levels to match order of grid patches
* MKS code unit system fixes
* Disabled bounds checking on pixelize_element_mesh
* Updated light_ray.py for domain width != 1
* Implemented a DOAP file generator
* Fixed bugs for 2D and 1D enzo IO
* Converted mutable Dataset attributes to be properties that return copies
* Allowing LightRay segments to extend further than one box length
* Fixed a divide-by-zero error that occasionally happens in
  triangle_plane_intersect
* Make sure we have an index in subclassed derived quantities
* Added an initial draft of an extensions document
* Made it possible to pass field tuples to command-line plotting
* Ensured the positions of coordinate vector lines are in code units
* Added a minus sign to definition of sz_kinetic field
* Added grid_levels and grid_indices fields to octrees
* Added a morton_index derived field
* Added Exception to AMRKDTree in the case of particle of oct-based data



Version 3.2
-----------

Major enhancements
^^^^^^^^^^^^^^^^^^

* Particle-Only Plots - a series of new plotting functions for visualizing
  particle data.  See here for more information.
* Late-stage beta support for Python 3 - unit tests and answer tests pass for
  all the major frontends under python 3.4, and yt should now be mostly if not
  fully usable.  Because many of the yt developers are still on Python 2 at
  this point, this should be considered a "late stage beta" as there may be
  remaining issues yet to be identified or worked out.
* Now supporting Gadget Friend-of-Friends/Subfind catalogs - see here to learn
  how to load halo catalogs as regular yt datasets.
* Custom colormaps can now be easily defined and added - see here to learn how!
* Now supporting Fargo3D data
* Performance improvements throughout the code base for memory and speed

Minor enhancements
^^^^^^^^^^^^^^^^^^

* Various updates to the following frontends: ART, Athena, Castro, Chombo,
  Gadget, GDF, Maestro, Pluto, RAMSES, Rockstar, SDF, Tipsy
* Numerous documentation updates
* Generic hexahedral mesh pixelizer
* Adding annotate_ray() callback for plots
* AbsorptionSpectrum returned to full functionality and now using faster SciPy
  Voigt profile
* Add a color_field argument to annotate_streamline
* Smoothing lengths auto-calculated for Tipsy Datasets
* Adding SimulationTimeSeries support for Gadget and OWLS.
* Generalizing derived quantity outputs to all be YTArrays or lists of
  YTArrays as appropriate
* Star analysis returned to full functionality
* FITS image writing refactor
* Adding gradient fields on the fly
* Adding support for Gadget Nx4 metallicity fields
* Updating value of solar metal mass fraction to be consistent with Cloudy.
* Gadget raw binary snapshot handling & non-cosmological simulation units
* Adding support for LightRay class to work with Gadget+Tipsy
* Add support for subclasses of frontends
* Dependencies updated
* Serialization for projections using minimal representation
* Adding Grid visitors in Cython
* Improved semantics for derived field units
* Add a yaw() method for the PerspectiveCamera + switch back to LHS
* Adding annotate_clear() function to remove previous callbacks from a plot
* Added documentation for hexahedral mesh on website
* Speed up nearest neighbor evaluation
* Add a convenience method to create deposited particle fields
* UI and docs updates for 3D streamlines
* Ensure particle fields are tested in the field unit tests
* Allow a suffix to be specified to save()
* Add profiling using airspeed velocity
* Various plotting enhancements and bugfixes
* Use hglib to update
* Various minor updates to halo_analysis toolkit
* Docker-based tests for install_script.sh
* Adding support for single and non-cosmological datasets to LightRay
* Adding the Pascal unit
* Add weight_field to PPVCube
* FITS reader: allow HDU in auxiliary
* Fixing electromagnetic units
* Specific Angular Momentum [xyz] computed relative to a normal vector

Bugfixes
^^^^^^^^

* Adding ability to create union fields from alias fields
* Small fix to allow enzo AP datasets to load in parallel when no APs present
* Use proper cell dimension in gradient function.
* Minor memory optimization for smoothed particle fields
* Fix thermal_energy for Enzo HydroMethod==6
* Make sure annotate_particles handles unitful widths properly
* Improvements for add_particle_filter and particle_filter
* Specify registry in off_axis_projection's image finalization
* Apply fix for particle momentum units to the boxlib frontend
* Avoid traceback in "yt version" when python-hglib is not installed
* Expose no_ghost from export_sketchfab down to _extract_isocontours_from_grid
* Fix broken magnetic_unit attribute
* Fixing an off-by-one error in the set x/y lim methods for profile plots
* Providing better error messages to PlotWindow callbacks
* Updating annotate_timestamp to avoid auto-override
* Updating callbacks to consistently define coordinate system
* Fixing species fields for OWLS and tipsy
* Fix extrapolation for vertex-centered data
* Fix periodicity check in FRBs
* Rewrote project_to_plane() in PerspectiveCamera for draw_domain()
* Fix intermittent failure in test_add_deposited_particle_field
* Improve minorticks for a symlog plot with one-sided data
* Fix smoothed covering grid cell computation
* Absorption spectrum generator now 3.0 compliant
* Fix off-by-one-or-more in particle smallest dx
* Fix dimensionality mismatch error in covering grid
* Fix curvature term in cosmology calculator
* Fix geographic axes and pixelization
* Ensure axes aspect ratios respect the user-selected plot aspect ratio
* Avoid clobbering field_map when calling profile.add_fields
* Fixing the arbitrary grid deposit code
* Fix spherical plotting centering
* Make the behavior of to_frb consistent with the docstring
* Ensure projected units are initialized when there are no chunks.
* Removing "field already exists" warnings from the Owls and Gadget frontends
* Various photon simulator bugs
* Fixed use of LaTeX math mode
* Fix upload_image
* Enforce plot width in CSS when displayed in a notebook
* Fix cStringIO.StringIO -> cStringIO in png_writer
* Add some input sanitizing and error checking to covering_grid initializer
* Fix for geographic plotting
* Use the correct filename template for single-file OWLS datasets.
* Fix Enzo IO performance for 32 bit datasets
* Adding a number density field for Enzo MultiSpecies=0 datasets.
* Fix RAMSES block ordering
* Updating ragged array tests for NumPy 1.9.1
* Force returning lists for HDF5FileHandler

Version 3.1
-----------

This is a scheduled feature release.  Below are the itemized, aggregate changes
since version 3.0.


Major changes:
^^^^^^^^^^^^^^

* The RADMC-3D export analysis module has been updated. `PR 1358 <https://bitbucket.org/yt_analysis/yt/pull-requests/1358>`_, `PR 1332 <https://bitbucket.org/yt_analysis/yt/pull-requests/1332>`_.

* Performance improvements for grid frontends. `PR 1350 <https://bitbucket.org/yt_analysis/yt/pull-requests/1350>`_. `PR 1382 <https://bitbucket.org/yt_analysis/yt/pull-requests/1382>`_, `PR 1322 <https://bitbucket.org/yt_analysis/yt/pull-requests/1322>`_.

* Added a frontend for Dark Matter-only NMSU Art simulations. `PR 1258 <https://bitbucket.org/yt_analysis/yt/pull-requests/1258>`_.

* The absorption spectrum generator has been updated. `PR 1356 <https://bitbucket.org/yt_analysis/yt/pull-requests/1356>`_.

* The PerspectiveCamera has been updated and a new SphericalCamera has been
  added. `PR 1346 <https://bitbucket.org/yt_analysis/yt/pull-requests/1346>`_, `PR 1299 <https://bitbucket.org/yt_analysis/yt/pull-requests/1299>`_.

* The unit system now supports unit equivalencies and has improved support for MKS units.  See :ref:`unit_equivalencies`. `PR 1291 <https://bitbucket.org/yt_analysis/yt/pull-requests/1291>`_, `PR 1286 <https://bitbucket.org/yt_analysis/yt/pull-requests/1286>`_.

* Data object selection can now be chained, allowing selecting based on multiple constraints. `PR 1264 <https://bitbucket.org/yt_analysis/yt/pull-requests/1264>`_.

* Added the ability to manually override the simulation unit system. `PR 1236 <https://bitbucket.org/yt_analysis/yt/pull-requests/1236>`_.

* The documentation has been reorganized and has seen substantial improvements. `PR 1383 <https://bitbucket.org/yt_analysis/yt/pull-requests/1383>`_, `PR 1373 <https://bitbucket.org/yt_analysis/yt/pull-requests/1373>`_, `PR 1364 <https://bitbucket.org/yt_analysis/yt/pull-requests/1364>`_, `PR 1351 <https://bitbucket.org/yt_analysis/yt/pull-requests/1351>`_, `PR 1345 <https://bitbucket.org/yt_analysis/yt/pull-requests/1345>`_. `PR 1333 <https://bitbucket.org/yt_analysis/yt/pull-requests/1333>`_, `PR 1342 <https://bitbucket.org/yt_analysis/yt/pull-requests/1342>`_, `PR 1338 <https://bitbucket.org/yt_analysis/yt/pull-requests/1338>`_, `PR 1330 <https://bitbucket.org/yt_analysis/yt/pull-requests/1330>`_, `PR 1326 <https://bitbucket.org/yt_analysis/yt/pull-requests/1326>`_, `PR 1323 <https://bitbucket.org/yt_analysis/yt/pull-requests/1323>`_, `PR 1315 <https://bitbucket.org/yt_analysis/yt/pull-requests/1315>`_, `PR 1305 <https://bitbucket.org/yt_analysis/yt/pull-requests/1305>`_, `PR 1289 <https://bitbucket.org/yt_analysis/yt/pull-requests/1289>`_, `PR 1276 <https://bitbucket.org/yt_analysis/yt/pull-requests/1276>`_.

Minor or bugfix changes:
^^^^^^^^^^^^^^^^^^^^^^^^

* The Ampere unit now accepts SI prefixes.  `PR 1393 <https://bitbucket.org/yt_analysis/yt/pull-requests/1393>`_.

* The Gadget InternalEnergy and StarFormationRate fields are now read in with the correct units.  `PR 1392 <https://bitbucket.org/yt_analysis/yt/pull-requests/1392>`_, `PR 1379 <https://bitbucket.org/yt_analysis/yt/pull-requests/1379>`_.

* Substantial improvements for the PPVCube analysis module and support for FITS dataset. `PR 1390 <https://bitbucket.org/yt_analysis/yt/pull-requests/1390>`_, `PR 1367 <https://bitbucket.org/yt_analysis/yt/pull-requests/1367>`_, `PR 1347 <https://bitbucket.org/yt_analysis/yt/pull-requests/1347>`_, `PR 1326 <https://bitbucket.org/yt_analysis/yt/pull-requests/1326>`_, `PR 1280 <https://bitbucket.org/yt_analysis/yt/pull-requests/1280>`_, `PR 1336 <https://bitbucket.org/yt_analysis/yt/pull-requests/1336>`_.

* The center of a PlotWindow plot can now be set to the maximum or minimum of any field. `PR 1280 <https://bitbucket.org/yt_analysis/yt/pull-requests/1280>`_.

* Fixes for yt testing infrastructure. `PR 1388 <https://bitbucket.org/yt_analysis/yt/pull-requests/1388>`_, `PR 1348 <https://bitbucket.org/yt_analysis/yt/pull-requests/1348>`_.

* Projections are now performed using an explicit path length field for all
  coordinate systems. `PR 1307 <https://bitbucket.org/yt_analysis/yt/pull-requests/1307>`_.

* An example notebook for simulations using the OWLS data format has been added
  to the documentation. `PR 1386 <https://bitbucket.org/yt_analysis/yt/pull-requests/1386>`_.

* Fix for the camera.draw_line function. `PR 1380 <https://bitbucket.org/yt_analysis/yt/pull-requests/1380>`_.

* Minor fixes and improvements for yt plots. `PR 1376 <https://bitbucket.org/yt_analysis/yt/pull-requests/1376>`_, `PR 1374 <https://bitbucket.org/yt_analysis/yt/pull-requests/1374>`_, `PR 1288 <https://bitbucket.org/yt_analysis/yt/pull-requests/1288>`_, `PR 1290 <https://bitbucket.org/yt_analysis/yt/pull-requests/1290>`_.

* Significant documentation reorganization and improvement. `PR 1375 <https://bitbucket.org/yt_analysis/yt/pull-requests/1375>`_, `PR 1359 <https://bitbucket.org/yt_analysis/yt/pull-requests/1359>`_.

* Fixed a conflict in the CFITSIO library used by the x-ray analysis module. `PR 1365 <https://bitbucket.org/yt_analysis/yt/pull-requests/1365>`_.

* Miscellaneous code cleanup. `PR 1371 <https://bitbucket.org/yt_analysis/yt/pull-requests/1371>`_, `PR 1361 <https://bitbucket.org/yt_analysis/yt/pull-requests/1361>`_.

* yt now hooks up to the python logging infrastructure in a more standard
  fashion, avoiding issues with yt logging showing up with using other
  libraries. `PR 1355 <https://bitbucket.org/yt_analysis/yt/pull-requests/1355>`_, `PR 1362 <https://bitbucket.org/yt_analysis/yt/pull-requests/1362>`_, `PR 1360 <https://bitbucket.org/yt_analysis/yt/pull-requests/1360>`_.

* The docstring for the projection data object has been corrected. `PR 1366 <https://bitbucket.org/yt_analysis/yt/pull-requests/1366>`_

* A bug in the calculation of the plot bounds for off-axis slice plots has been fixed. `PR 1357 <https://bitbucket.org/yt_analysis/yt/pull-requests/1357>`_.

* Improvements for the yt-rockstar interface. `PR 1352 <https://bitbucket.org/yt_analysis/yt/pull-requests/1352>`_, `PR 1317 <https://bitbucket.org/yt_analysis/yt/pull-requests/1317>`_.

* Fix issues with plot positioning with saving to postscript or encapsulated postscript. `PR 1353 <https://bitbucket.org/yt_analysis/yt/pull-requests/1353>`_.

* It is now possible to supply a default value for get_field_parameter. `PR 1343 <https://bitbucket.org/yt_analysis/yt/pull-requests/1343>`_.

* A bug in the interpretation of the units of RAMSES simulations has been fixed. `PR 1335 <https://bitbucket.org/yt_analysis/yt/pull-requests/1335>`_.

* Plot callbacks are now only executed once before the plot is saved. `PR 1328 <https://bitbucket.org/yt_analysis/yt/pull-requests/1328>`_.

* Performance improvements for smoothed covering grid alias fields. `PR 1331 <https://bitbucket.org/yt_analysis/yt/pull-requests/1331>`_.

* Improvements and bugfixes for the halo analysis framework. `PR 1349 <https://bitbucket.org/yt_analysis/yt/pull-requests/1349>`_, `PR 1325 <https://bitbucket.org/yt_analysis/yt/pull-requests/1325>`_.

* Fix issues with the default setting for the ``center`` field parameter. `PR 1327 <https://bitbucket.org/yt_analysis/yt/pull-requests/1327>`_.

* Avoid triggering warnings in numpy and matplotlib. `PR 1334 <https://bitbucket.org/yt_analysis/yt/pull-requests/1334>`_, `PR 1300 <https://bitbucket.org/yt_analysis/yt/pull-requests/1300>`_.

* Updates for the field list reference. `PR 1344 <https://bitbucket.org/yt_analysis/yt/pull-requests/1344>`_, `PR 1321 <https://bitbucket.org/yt_analysis/yt/pull-requests/1321>`_, `PR 1318 <https://bitbucket.org/yt_analysis/yt/pull-requests/1318>`_.

* yt can now be run in parallel on a subset of available processors using an MPI subcommunicator. `PR 1340 <https://bitbucket.org/yt_analysis/yt/pull-requests/1340>`_

* Fix for incorrect units when loading an Athena simulation as a time series. `PR 1341 <https://bitbucket.org/yt_analysis/yt/pull-requests/1341>`_.

* Improved support for Enzo 3.0 simulations that have not produced any active particles. `PR 1329 <https://bitbucket.org/yt_analysis/yt/pull-requests/1329>`_.

* Fix for parsing OWLS outputs with periods in the file path.  `PR 1320 <https://bitbucket.org/yt_analysis/yt/pull-requests/1320>`_.

* Fix for periodic radius vector calculation. `PR 1311 <https://bitbucket.org/yt_analysis/yt/pull-requests/1311>`_.

* Improvements for the Maestro and Castro frontends. `PR 1319 <https://bitbucket.org/yt_analysis/yt/pull-requests/1319>`_.

* Clump finding is now supported for more generic types of data. `PR 1314 <https://bitbucket.org/yt_analysis/yt/pull-requests/1314>`_

* Fix unit consistency issue when mixing dimensionless unit symbols. `PR 1300 <https://bitbucket.org/yt_analysis/yt/pull-requests/1300>`_.

* Improved memory footprint in the photon_simulator. `PR 1304 <https://bitbucket.org/yt_analysis/yt/pull-requests/1304>`_.

* Large grids in Athena datasets produced by the join_vtk script can now be optionally split, improving parallel performance.  `PR 1304 <https://bitbucket.org/yt_analysis/yt/pull-requests/1304>`_.

* Slice plots now accept a ``data_source`` keyword argument. `PR 1310 <https://bitbucket.org/yt_analysis/yt/pull-requests/1310>`_.

* Corrected inconsistent octrees in the RAMSES frontend. `PR 1302 <https://bitbucket.org/yt_analysis/yt/pull-requests/1302>`_

* Nearest neighbor distance field added.  `PR 1138 <https://bitbucket.org/yt_analysis/yt/pull-requests/1138>`_.

* Improvements for the ORION2 frontend. `PR 1303 <https://bitbucket.org/yt_analysis/yt/pull-requests/1303>`_

* Enzo 3.0 frontend can now read active particle attributes that are arrays of any shape. `PR 1248 <https://bitbucket.org/yt_analysis/yt/pull-requests/1248>`_.

* Answer tests added for halo finders. `PR 1253 <https://bitbucket.org/yt_analysis/yt/pull-requests/1253>`_

* A ``setup_function`` has been added to the LightRay initializer. `PR 1295 <https://bitbucket.org/yt_analysis/yt/pull-requests/1295>`_.

* The SPH code frontends have been reorganized into separate frontend directories. `PR 1281 <https://bitbucket.org/yt_analysis/yt/pull-requests/1281>`_.

* Fixes for accessing deposit fields for FLASH data. `PR 1294 <https://bitbucket.org/yt_analysis/yt/pull-requests/1294>`_

* Added tests for ORION datasets containing sink and star particles. `PR 1252 <https://bitbucket.org/yt_analysis/yt/pull-requests/1252>`_

* Fix for field names in the particle generator. `PR 1278 <https://bitbucket.org/yt_analysis/yt/pull-requests/1278>`_.

* Added wrapper functions for numpy array manipulation functions.  `PR 1287 <https://bitbucket.org/yt_analysis/yt/pull-requests/1287>`_.

* Added support for packed HDF5 Enzo datasets. `PR 1282 <https://bitbucket.org/yt_analysis/yt/pull-requests/1282>`_.

Version 3.0
-----------

This release of yt features an entirely rewritten infrastructure for
data ingestion, indexing, and representation.  While past versions of
yt were focused on analysis and visualization of data structured as
regular grids, this release features full support for particle
(discrete point) data such as N-body and SPH data, irregular
hexahedral mesh data, and data organized via octrees.  This
infrastructure will be extended in future versions for high-fidelity
representation of unstructured mesh datasets.

Highlighted changes in yt 3.0:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 * Units now permeate the code base, enabling self-consistent unit
   transformations of all arrays and quantities returned by yt.
 * Particle data is now supported using a lightweight octree.  SPH
   data can be smoothed onto an adaptively-defined mesh using standard
   SPH smoothing
 * Support for octree AMR codes
 * Preliminary Support for non-Cartesian data, such as cylindrical,
   spherical, and geographical
 * Revamped analysis framework for halos and halo catalogs, including
   direct ingestion and analysis of halo catalogs of several different
   formats
 * Support for multi-fluid datasets and datasets containing multiple
   particle types
 * Flexible support for dynamically defining new particle types using
   filters on existing particle types or by combining different particle
   types.
 * Vastly improved support for loading generic grid, AMR, hexahedral
   mesh, and particle without hand-coding a frontend for a particular
   data format.
 * New frontends for ART, ARTIO, Boxlib, Chombo, FITS, GDF, Subfind,
   Rockstar, Pluto, RAMSES, SDF, Gadget, OWLS, PyNE, Tipsy, as well as
   rewritten frontends for Enzo, FLASH, Athena, and generic data.
 * First release to support installation of yt on Windows
 * Extended capabilities for construction of simulated observations,
   and new facilities for analyzing and visualizing FITS images and cube
   data
 * Many performance improvements

This release is the first of several; while most functionality from
the previous generation of yt has been updated to work with yt 3.0, it
does not yet have feature parity in all respects.  While the core of
yt is stable, we suggest the support for analysis modules and volume
rendering be viewed as a late-stage beta, with a series of additional
releases (3.1, 3.2, etc) appearing over the course of the next year to
improve support in these areas.

For a description of how to bring your 2.x scripts up to date to 3.0,
and a summary of common gotchas in this transition, please see
:ref:`yt3differences`.

Version 2.6
-----------

This is a scheduled release, bringing to a close the development in the 2.x
series.  Below are the itemized, aggregate changes since version 2.5.

Major changes:
^^^^^^^^^^^^^^

  * yt is now licensed under the 3-clause BSD license.
  * HEALPix has been removed for the time being, as a result of licensing
    incompatibility.
  * The addition of a frontend for the Pluto code
  * The addition of an OBJ exporter to enable transparent and multi-surface
    exports of surfaces to Blender and Sketchfab
  * New absorption spectrum analysis module with documentation
  * Adding ability to draw lines with Grey Opacity in volume rendering
  * Updated physical constants to reflect 2010 CODATA data
  * Dependency updates (including IPython 1.0)
  * Better notebook support for yt plots
  * Considerably (10x+) faster kD-tree building for volume rendering
  * yt can now export to RADMC3D
  * Athena frontend now supports Static Mesh Refinement and units (
    http://hub.yt-project.org/nb/7l1zua )
  * Fix long-standing bug for plotting arrays with range of zero
  * Adding option to have interpolation based on non-uniform bins in
    interpolator code
  * Upgrades to most of the dependencies in the install script
  * ProjectionPlot now accepts a data_source keyword argument

Minor or bugfix changes:
^^^^^^^^^^^^^^^^^^^^^^^^

  * Fix for volume rendering on the command line
  * map_to_colormap will no longer return out-of-bounds errors
  * Fixes for dds in covering grid calculations
  * Library searching for build process is now more reliable
  * Unit fix for "VorticityGrowthTimescale" field
  * Pyflakes stylistic fixes
  * Number density added to FLASH
  * Many fixes for Athena frontend
  * Radius and ParticleRadius now work for reduced-dimensionality datasets
  * Source distributions now work again!
  * Athena data now 64 bits everywhere
  * Grids displays on plots are now shaded to reflect the level of refinement
  * show_colormaps() is a new function for displaying all known colormaps
  * PhasePlotter by default now adds a colormap.
  * System build fix for POSIX systems
  * Fixing domain offsets for halo centers-of-mass
  * Removing some Enzo-specific terminology in the Halo Mass Function
  * Addition of coordinate vectors on volume render
  * Pickling fix for extracted regions
  * Addition of some tracer particle annotation functions
  * Better error message for "yt" command
  * Fix for radial vs poloidal fields
  * Piernik 2D data handling fix
  * Fixes for FLASH current redshift
  * PlotWindows now have a set_font function and a new default font setting
  * Colorbars less likely to extend off the edge of a PlotWindow
  * Clumps overplotted on PlotWindows are now correctly contoured
  * Many fixes to light ray and profiles for integrated cosmological analysis
  * Improvements to OpenMP compilation
  * Typo in value for km_per_pc (not used elsewhere in the code base) has been
    fixed
  * Enable parallel IPython notebook sessions (
    http://hub.yt-project.org/nb/qgn19h )
  * Change (~1e-6) to particle_density deposition, enabling it to be used by
    FLASH and other frontends
  * Addition of is_root function for convenience in parallel analysis sessions
  * Additions to Orion particle reader
  * Fixing TotalMass for case when particles not present
  * Fixing the density threshold or HOP and pHOP to match the merger tree
  * Reason can now plot with latest plot window
  * Issues with VelocityMagnitude and aliases with velo have been corrected in
    the FLASH frontend
  * Halo radii are calculated correctly for domains that do not start at 0,0,0.
  * Halo mass function now works for non-Enzo frontends.
  * Bug fixes for directory creation, typos in docstrings
  * Speed improvements to ellipsoidal particle detection
  * Updates to FLASH fields
  * CASTRO frontend bug fixes
  * Fisheye camera bug fixes
  * Answer testing now includes plot window answer testing
  * Athena data serialization
  * load_uniform_grid can now decompose dims >= 1024.  (#537)
  * Axis unit setting works correctly for unit names  (#534)
  * ThermalEnergy is now calculated correctly for Enzo MHD simulations (#535)
  * Radius fields had an asymmetry in periodicity calculation (#531)
  * Boolean regions can now be pickled (#517)

Version 2.5
-----------

Many below-the-surface changes happened in yt 2.5 to improve reliability,
fidelity of the answers, and streamlined user interface.  The major change in
this release has been the immense expansion in testing of yt.  We now have over
2000 unit tests (run on every commit, thanks to both Kacper Kowalik and Shining
Panda) as well as answer testing for FLASH, Enzo, Chombo and Orion data.

The Stream frontend, which can construct datasets in memory, has been improved
considerably.  It's now easier than ever to load data from disk.  If you know
how to get volumetric data into Python, you can use either the
``load_uniform_grid`` function or the ``load_amr_grid`` function to create an
in-memory dataset that yt can analyze.

yt now supports the Athena code.

yt is now focusing on providing first class support for the IPython notebook.
In this release, plots can be displayed inline.  The Reason HTML5 GUI will be
merged with the IPython notebook in a future release.

Install Script Changes:
^^^^^^^^^^^^^^^^^^^^^^^

 * SciPy can now be installed
 * Rockstar can now be installed
 * Dependencies can be updated with "yt update --all"
 * Cython has been upgraded to 0.17.1
 * Python has been upgraded to 2.7.3
 * h5py has been upgraded to 2.1.0
 * hdf5 has been upgraded to 1.8.9
 * matplotlib has been upgraded to 1.2.0
 * IPython has been upgraded to 0.13.1
 * Forthon has been upgraded to 0.8.10
 * nose has been added
 * sympy has been added
 * python-hglib has been added

We've also improved support for installing on OSX, Ubuntu and OpenSUSE.

Most Visible Improvements
^^^^^^^^^^^^^^^^^^^^^^^^^

 * Nearly 200 pull requests and over 1000 changesets have been merged since yt
   2.4 was release on August 2nd, 2012.
 * numpy is now imported as np, not na.  na will continue to work for the
   foreseeable future.
 * You can now get a `yt cheat sheet <http://yt-project.org/docs/2.5/cheatsheet.pdf>`!
 * yt can now load simulation data created by Athena.
 * The Rockstar halo finder can now be installed by the install script
 * SciPy can now be installed by the install script
 * Data can now be written out in two ways:

   * Sidecar files containing expensive derived fields can be written and
     implicitly loaded from.
   * GDF files, which are portable yt-specific representations of full
     simulations, can be created from any dataset.  Work is underway on
     a pure C library that can be linked against to load these files into
     simulations.

 * The "Stream" frontend, for loading raw data in memory, has been greatly
   expanded and now includes initial conditions generation functionality,
   particle fields, and simple loading of AMR grids with ``load_amr_grids``.
 * Spherical and Cylindrical fields have been sped up and made to have a
   uniform interface.  These fields can be the building blocks of more advanced
   fields.
 * Coordinate transformations have been sped up and streamlined. It is now
   possible to convert any scalar or vector field to a new cartesian, spherical,
   or cylindrical coordinate system with an arbitrary orientation. This makes it
   possible to do novel analyses like profiling the toroidal and poloidal
   velocity as a function of radius in an inclined disk.
 * Many improvements to the EnzoSimulation class, which can now find many
   different types of data.
 * Image data is now encapsulated in an ImageArray class, which carries with it
   provenance information about its trajectory through yt.
 * Streamlines now query at every step along the streamline, not just at every
   cell.
 * Surfaces can now be extracted and examined, as well as uploaded to
   Sketchfab.com for interactive visualization in a web browser.
 * allsky_projection can now accept a datasource, making it easier to cut out
   regions to examine.
 * Many, many improvements to PlotWindow.  If you're still using
   PlotCollection, check out ``ProjectionPlot``, ``SlicePlot``,
   ``OffAxisProjectionPlot`` and ``OffAxisSlicePlot``.
 * PlotWindow can now accept a timeseries instead of a dataset.
 * Many fixes for 1D and 2D data, especially in FLASH datasets.
 * Vast improvements to the particle file handling for FLASH datasets.
 * Particles can now be created ex nihilo with CICSample_3.
 * Rockstar halo finding is now a targeted goal.  Support for using Rockstar
   has improved dramatically.
 * Increased support for tracking halos across time using the FOF halo finder.
 * The command ``yt notebook`` has been added to spawn an IPython notebook
   server, and the ``yt.imods`` module can replace ``yt.mods`` in the IPython
   Notebook to enable better integration.
 * Metallicity-dependent X-ray fields have now been added.
 * Grid lines can now be added to volume renderings.
 * Volume rendering backend has been updated to use an alpha channel, fixing
   parallel opaque volume renderings.  This also enables easier blending of
   multiple images and annotations to the rendering. Users are encouraged
   to look at the capabilities of the ``ImageArray`` for writing out renders,
   as updated in the cookbook examples. Volume renders can now be saved with
   an arbitrary background color.
 * Periodicity, or alternately non-periodicity, is now a part of radius
   calculations.
 * The AMRKDTree has been rewritten.  This allows parallelism with other than
   power-of-2 MPI processes, arbitrary sets of grids, and splitting of
   unigrids.
 * Fixed Resolution Buffers and volume rendering images now utilize a new
   ImageArray class that stores information such as data source, field names,
   and other information in a .info dictionary. See the ``ImageArray``
   docstrings for more information on how they can be used to save to a bitmap
   or hdf5 file.

Version 2.4
-----------

The 2.4 release was particularly large, encompassing nearly a thousand
changesets and a number of new features.

To help you get up to speed, we've made an IPython notebook file demonstrating
a few of the changes to the scripting API.  You can
`download it here <http://yt-project.org/files/yt24.ipynb>`_.

Most Visible Improvements
^^^^^^^^^^^^^^^^^^^^^^^^^

 * Threaded volume renderer, completely refactored from the ground up for
   speed and parallelism.
 * The Plot Window (see :ref:`simple-inspection`) is now fully functional!  No
   more PlotCollections, and full, easy access to Matplotlib axes objects.
 * Many improvements to Time Series analysis:
    * EnzoSimulation now integrates with TimeSeries analysis!
    * Auto-parallelization of analysis and parallel iteration
    * Memory usage when iterating over datasets reduced substantially
 * Many improvements to Reason, the yt GUI
    * Addition of "yt reason" as a startup command
    * Keyboard shortcuts in projection & slice mode: z, Z, x, X for zooms,
      hjkl, HJKL for motion
    * Drag to move in projection & slice mode
    * Contours and vector fields in projection & slice mode
    * Color map selection in projection & slice mode
    * 3D Scene
 * Integration with the all new yt Hub ( http://hub.yt-project.org/ ): upload
   variable resolution projections, slices, project information, vertices and
   plot collections right from the yt command line!

Other Changes
^^^^^^^^^^^^^

 * :class:`~yt.visualization.plot_window.ProjectionPlot` and
   :class:`~yt.visualization.plot_window.SlicePlot` supplant the functionality
   of PlotCollection.
 * Camera path creation from keyframes and splines
 * Ellipsoidal data containers and ellipsoidal parameter calculation for halos
 * PyX and ZeroMQ now available in the install script
 * Consolidation of unit handling
 * HDF5 updated to 1.8.7, Mercurial updated to 2.2, IPython updated to 0.12
 * Preview of integration with Rockstar halo finder
 * Improvements to merger tree speed and memory usage
 * Sunrise exporter now compatible with Sunrise 4.0
 * Particle trajectory calculator now available!
 * Speed and parallel scalability improvements in projections, profiles and HOP
 * New Vorticity-related fields
 * Vast improvements to the ART frontend
 * Many improvements to the FLASH frontend, including full parameter reads,
   speedups, and support for more corner cases of FLASH 2, 2.5 and 3 data.
 * Integration of the Grid Data Format frontend, and a converter for Athena
   data to this format.
 * Improvements to command line parsing
 * Parallel import improvements on parallel filesystems
   (``from yt.pmods import *``)
 * proj_style keyword for projections, for Maximum Intensity Projections
   (``proj_style = "mip"``)
 * Fisheye rendering for planetarium rendering
 * Profiles now provide \*_std fields for standard deviation of values
 * Generalized Orientation class, providing 6DOF motion control
 * parallel_objects iteration now more robust, provides optional barrier.
   (Also now being used as underlying iteration mechanism in many internal
   routines.)
 * Dynamic load balancing in parallel_objects iteration.
 * Parallel-aware objects can now be pickled.
 * Many new colormaps included
 * Numerous improvements to the PyX-based eps_writer module
 * FixedResolutionBuffer to FITS export.
 * Generic image to FITS export.
 * Multi-level parallelism for extremely large cameras in volume rendering
 * Light cone and light ray updates to fit with current best practices for
   parallelism

Version 2.3
-----------

`(yt 2.3 docs) <http://yt-project.org/docs/2.3>`_
 * Multi-level parallelism
 * Real, extensive answer tests
 * Boolean data regions (see :ref:`boolean_data_objects`)
 * Isocontours / flux calculations (see :ref:`extracting-isocontour-information`)
 * Field reorganization
 * PHOP memory improvements
 * Bug fixes for tests
 * Parallel data loading for RAMSES, along with other speedups and improvements
   there
 * WebGL interface for isocontours and a pannable map widget added to Reason
 * Performance improvements for volume rendering
 * Adaptive HEALPix support
 * Column density calculations
 * Massive speedup for 1D profiles
 * Lots more, bug fixes etc.
 * Substantial improvements to the documentation, including
   :ref:`manual-plotting` and a revamped orientation.

Version 2.2
-----------

`(yt 2.2 docs) <http://yt-project.org/docs/2.2>`_
 * Command-line submission to the yt Hub (http://hub.yt-project.org/)
 * Initial release of the web-based GUI Reason, designed for efficient remote
   usage over SSH tunnels
 * Absorption line spectrum generator for cosmological simulations (see
   :ref:`absorption_spectrum`)
 * Interoperability with ParaView for volume rendering, slicing, and so forth
 * Support for the Nyx code
 * An order of magnitude speed improvement in the RAMSES support
 * Quad-tree projections, speeding up the process of projecting by up to an
   order of magnitude and providing better load balancing
 * "mapserver" for in-browser, Google Maps-style slice and projection
   visualization (see :ref:`mapserver`)
 * Many bug fixes and performance improvements
 * Halo loader

Version 2.1
-----------

`(yt 2.1 docs) <http://yt-project.org/docs/2.1>`_
 * HEALPix-based volume rendering for 4pi, allsky volume rendering
 * libconfig is now included
 * SQLite3 and Forthon now included by default in the install script
 * Development guide has been lengthened substantially and a development
   bootstrap script is now included.
 * Installation script now installs Python 2.7 and HDF5 1.8.6
 * iyt now tab-completes field names
 * Halos can now be stored on-disk much more easily between HaloFinding runs.
 * Halos found inline in Enzo can be loaded and merger trees calculated
 * Support for CASTRO particles has been added
 * Chombo support updated and fixed
 * New code contributions
 * Contour finder has been sped up by a factor of a few
 * Constrained two-point functions are now possible, for LOS power spectra
 * Time series analysis (:ref:`time-series-analysis`) now much easier
 * Stream Lines now a supported 1D data type
 * Stream Lines now able to be calculated and plotted (:ref:`streamlines`)
 * In situ Enzo visualization now much faster
 * "gui" source directory reorganized and cleaned up
 * Cython now a compile-time dependency, reducing the size of source tree
   updates substantially
 * ``yt-supplemental`` repository now checked out by default, containing
   cookbook, documentation, handy mercurial extensions, and advanced plotting
   examples and helper scripts.
 * Pasteboards now supported and available
 * Parallel yt efficiency improved by removal of barriers and improvement of
   collective operations

Version 2.0
-----------

 * Major reorganization of the codebase for speed, ease of modification, and maintainability
 * Re-organization of documentation and addition of Orientation Session
 * Support for FLASH code
 * Preliminary support for MAESTRO, CASTRO, ART, and RAMSES (contributions welcome!)
 * Perspective projection for volume rendering
 * Exporting to Sunrise
 * Preliminary particle rendering in volume rendering visualization
 * Drastically improved parallel volume rendering, via kD-tree decomposition
 * Simple merger tree calculation for FOF catalogs
 * New and greatly expanded documentation, with a "source" button

Version 1.7
-----------

 * Direct writing of PNGs
 * Multi-band image writing
 * Parallel halo merger tree (see :ref:`merger_tree`)
 * Parallel structure function generator (see :ref:`two_point_functions`)
 * Image pan and zoom object and display widget.
 * Parallel volume rendering (see :ref:`volume_rendering`)
 * Multivariate volume rendering, allowing for multiple forms of emission and
   absorption, including approximate scattering and Planck emissions. (see
   :ref:`volume_rendering`)
 * Added Camera interface to volume rendering (See :ref:`volume_rendering`)
 * Off-axis projection (See :ref:`volume_rendering`)
 * Stereo (toe-in) volume rendering (See :ref:`volume_rendering`)
 * DualEPS extension for better EPS construction
 * yt now uses Distribute instead of SetupTools
 * Better ``iyt`` initialization for GUI support
 * Rewritten, memory conservative and speed-improved contour finding algorithm
 * Speed improvements to volume rendering
 * Preliminary support for the Tiger code
 * Default colormap is now ``algae``
 * Lightweight projection loading with ``projload``
 * Improvements to `yt.data_objects.time_series`
 * Improvements to :class:`yt.extensions.EnzoSimulation` (See
   :ref:`analyzing-an-entire-simulation`)
 * Removed ``direct_ray_cast``
 * Fixed bug causing double data-read in projections
 * Added Cylinder support to ParticleIO
 * Fixes for 1- and 2-D Enzo datasets
 * Preliminary, largely non-functional Gadget support
 * Speed improvements to basic HOP
 * Added physical constants module
 * Beginning to standardize and enforce docstring requirements, changing to
   ``autosummary``-based API documentation.

Version 1.6.1
-------------

 * Critical fixes to ParticleIO
 * Halo mass function fixes for comoving coordinates
 * Fixes to halo finding
 * Fixes to the installation script
 * "yt instinfo" command to report current installation information as well as
   auto-update some types of installations
 * Optimizations to the volume renderer (2x-26x reported speedups)

Version 1.6
-----------

Version 1.6 is a point release, primarily notable for the new parallel halo
finder (see :ref:`halo-analysis`)

 * (New) Parallel HOP ( http://arxiv.org/abs/1001.3411 , :ref:`halo-analysis` )
 * (Beta) Software ray casting and volume rendering
   (see :ref:`volume_rendering`)
 * Rewritten, faster and better contouring engine for clump identification
 * Spectral Energy Distribution calculation for stellar populations
   (see :ref:`synthetic_spectrum`)
 * Optimized data structures such as the index
 * Star particle analysis routines
   (see :ref:`star_analysis`)
 * Halo mass function routines
 * Completely rewritten, massively faster and more memory efficient Particle IO
 * Fixes for plots, including normalized phase plots
 * Better collective communication in parallel routines
 * Consolidation of optimized C routines into ``amr_utils``
 * Many bug fixes and minor optimizations

Version 1.5
-----------

Version 1.5 features many new improvements, most prominently that of the
addition of parallel computing abilities (see :ref:`parallel-computation`) and
generalization for multiple AMR data formats, specifically both Enzo and Orion.

 * Rewritten documentation
 * Fully parallel slices, projections, cutting planes, profiles,
   quantities
 * Parallel HOP
 * Friends-of-friends halo finder
 * Object storage and serialization
 * Major performance improvements to the clump finder (factor of five)
 * Generalized domain sizes
 * Generalized field info containers
 * Dark Matter-only simulations
 * 1D and 2D simulations
 * Better IO for HDF5 sets
 * Support for the Orion AMR code
 * Spherical re-gridding
 * Halo profiler
 * Disk image stacker
 * Light cone generator
 * Callback interface improved
 * Several new callbacks
 * New data objects -- ortho and non-ortho rays, limited ray-tracing
 * Fixed resolution buffers
 * Spectral integrator for CLOUDY data
 * Substantially better interactive interface
 * Performance improvements *everywhere*
 * Command-line interface to *many* common tasks
 * Isolated plot handling, independent of PlotCollections

Version 1.0
-----------

 * Initial release!
