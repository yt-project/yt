.. _filtering-data:

Filtering your Dataset
======================

Large datasets are oftentimes too overwhelming to deal with in their
entirety, and it can be useful and faster
to analyze subsets of these datasets.  Furthermore, filtering the dataset
based on some field condition can reveal subtle information not easily
accessible by looking at the whole dataset.
Filters can be generated based on spatial position, say in a sphere
in the center of your dataset space, or more generally they can be
defined by the properties of any field in the simulation.

Because *mesh fields* are internally different from *particle fields*,
there are different ways of filtering each type as indicated below;
however, filtering fields by spatial location (i.e. geometric
objects) will apply to both types equally.

.. _filtering-mesh:

Filtering Mesh Fields
----------------------

Mesh fields can be filtered by two methods: cut region objects
(:class:`~yt.data_objects.selection_data_containers.YTCutRegion`)
and NumPy boolean masks.  Boolean masks are simpler, but they only work
for examining datasets, whereas cut regions objects create wholly new
data objects suitable for full analysis (data examination, image generation,
etc.)

Boolean Masks
^^^^^^^^^^^^^

NumPy boolean masks can be used with any NumPy array simply by passing the
array a conditional.  As a general example of this:

.. notebook-cell::

    import numpy as np

    a = np.arange(5)
    bigger_than_two = a > 2
    print("Original Array: a = \n%s" % a)
    print("Boolean Mask: bigger_than_two = \n%s" % bigger_than_two)
    print("Masked Array: a[bigger_than_two] = \n%s" % a[bigger_than_two])

Similarly, if you've created a yt data object (e.g. a region, a sphere), you
can examine its field values as a NumPy array by simply indexing it with the
field name.  Thus, it too can be masked using a NumPy boolean mask.  Let's
set a simple mask based on the contents of one of our fields.

.. notebook-cell::

    import yt

    ds = yt.load("Enzo_64/DD0042/data0042")
    ad = ds.all_data()
    hot = ad["gas", "temperature"].in_units("K") > 1e6
    print(
        'Temperature of all data: ad["gas", "temperature"] = \n%s'
        % ad["gas", "temperature"]
    )
    print("Boolean Mask: hot = \n%s" % hot)
    print(
        'Temperature of "hot" data: ad["gas", "temperature"][hot] = \n%s'
        % ad["gas", "temperature"][hot]
    )

This was a simple example, but one can make the conditionals that define
a boolean mask have multiple parts, and one can stack masks together to
make very complex cuts on one's data.  Once the data is filtered, it can be
used if you simply need to access the NumPy arrays:

.. notebook-cell::

    import yt

    ds = yt.load("Enzo_64/DD0042/data0042")
    ad = ds.all_data()
    overpressure_and_fast = (
        (ad["gas", "pressure"] > 1e-14) &
        (ad["gas", "velocity_magnitude"].in_units("km/s") > 1e2)
    )
    density = ad["gas", "density"]
    print('Density of all data: ad["gas", "density"] = \n%s' % density)
    print(
        'Density of "overpressure and fast" data: overpressure_and_fast["gas", "density"] = \n%s'
        % density[overpressure_and_fast]
    )

.. _cut-regions:

Cut Regions
^^^^^^^^^^^

Cut regions are a more general solution to filtering mesh fields.  The output
of a cut region is an entirely new data object, which can be treated like any
other data object to generate images, examine its values, etc.

.. notebook:: mesh_filter.ipynb

In addition to inputting string parameters into cut_region to specify filters,
wrapper functions exist that allow the user to use a simplified syntax for
filtering out unwanted regions. Such wrapper functions are methods of
:func: ``YTSelectionContainer3D``.

.. notebook-cell::

   import yt

   ds = yt.load("Enzo_64/DD0042/data0042")
   ad = ds.all_data()
   overpressure_and_fast = ad.include_above(("gas", "pressure"), 1e-14)
   # You can chain include_xx and exclude_xx to produce the intersection of cut regions
   overpressure_and_fast = overpressure_and_fast.include_above(
       ("gas", "velocity_magnitude"), 1e2, "km/s"
   )

   print('Density of all data: ad["gas", "density"] = \n%s' % ad["gas", "density"])
   print(
       'Density of "overpressure and fast" data: overpressure_and_fast["gas", "density"] = \n%s'
       % overpressure_and_fast["gas", "density"]
   )

The following exclude and include functions are supported:
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.include_equal` - Only include values equal to given value
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_equal`- Exclude values equal to given value
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.include_inside` - Only include values inside closed interval
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_inside` - Exclude values inside closed interval
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.include_outside` - Only include values outside closed interval
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_outside` - Exclude values outside closed interval
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_nan` - Exclude NaN values
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.include_above` - Only include values above given value
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_above` - Exclude values above given value
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.include_below` - Only include values below given balue
   - :func:`~yt.data_objects.data_containers.YTSelectionContainer3D.exclude_below` - Exclude values below given value


.. warning::

    Cut regions are unstable when used on particle fields. Though you can create
    a cut region using a mesh field or fields as a filter and then obtain a
    particle field within that region, you cannot create a cut region using
    particle fields in the filter, as yt will currently raise an error. If
    you want to filter particle fields, see the next section
    :ref:`filtering-particles` instead.

.. _filtering-particles:

Filtering Particle Fields
-------------------------

Particle filters create new particle fields based on the manipulation and
cuts on existing particle fields.  You can apply cuts to them to effectively
mask out everything except the particles with which you are concerned.

Creating a particle filter takes a few steps.  You must first define a
function which accepts a data object (e.g. all_data, sphere, etc.)
as its argument.  It uses the fields and information in this geometric
object in order to produce some sort of conditional mask that is then returned
to create a new particle type.

Here is a particle filter to create a new ``star`` particle type.  For Enzo
simulations, stars have ``particle_type`` set to 2, so our filter will select
only the particles with ``particle_type`` (i.e.  field = ``('all',
'particle_type')`` equal to 2.

.. code-block:: python

    @yt.particle_filter(requires=["particle_type"], filtered_type="all")
    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 2
        return filter

The :func:`~yt.data_objects.particle_filters.particle_filter` decorator takes a
few options.  You must specify the names of the particle fields that are
required in order to define the filter --- in this case the ``particle_type``
field.  Additionally, you must specify the particle type to be filtered --- in
this case we filter all the particle in dataset by specifying the ``all``
particle type.

In addition, you may specify a name for the newly defined particle type.  If no
name is specified, the name for the particle type will be inferred from the name
of the filter definition --- in this case the inferred name will be ``stars``.

As an alternative syntax, you can also define a new particle filter via the
:func:`~yt.data_objects.particle_filter.add_particle_filter` function.

.. code-block:: python

    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 2
        return filter


    yt.add_particle_filter(
        "stars", function=stars, filtered_type="all", requires=["particle_type"]
    )

This is equivalent to our use of the ``particle_filter`` decorator above.  The
choice to use either the ``particle_filter`` decorator or the
``add_particle_filter`` function is a purely stylistic choice.

Lastly, the filter must be applied to our dataset of choice.  Note that this
filter can be added to as many datasets as we wish.  It will only actually
create new filtered fields if the dataset has the required fields, though.

.. code-block:: python

    import yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    ds.add_particle_filter("stars")

And that's it!  We can now access all of the ('stars', field) fields from
our dataset ``ds`` and treat them as any other particle field.  In addition,
it created some ``deposit`` fields, where the particles were deposited on to
the grid as mesh fields.

We can create additional filters building on top of the filters we have.
For example, we can identify the young stars based on their age, which is
the difference between current time and their creation_time.

.. code-block:: python

    def young_stars(pfilter, data):
        age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
        filter = np.logical_and(age.in_units("Myr") <= 5, age >= 0)
        return filter


    yt.add_particle_filter(
        "young_stars",
        function=young_stars,
        filtered_type="stars",
        requires=["creation_time"],
    )

If we properly define all the filters using the decorator ``yt.particle_filter``
or the function ``yt.add_particle_filter`` in advance. We can add the filter
we need to the dataset. If the ``filtered_type`` is already defined but not
added to the dataset, it will automatically add the filter first. For example,
if we add the ``young_stars`` filter, which is filtered from ``stars``,
to the dataset, it will also add ``stars`` filter to the dataset.

.. code-block:: python

    import yt

    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    ds.add_particle_filter("young_stars")


.. notebook:: particle_filter.ipynb

.. _particle-unions:

Particle Unions
---------------

Multiple types of particles can be combined into a single, conceptual type.  As
an example, the NMSU-ART code has multiple "species" of dark matter, which we
union into a single ``darkmatter`` field.  The ``all`` particle type is a
special case of this.

To create a particle union, you need to import the ``ParticleUnion`` class from
``yt.data_objects.particle_unions``, which you then create and pass into
``add_particle_union`` on a dataset object.

Here is an example, where we union the ``halo`` and ``disk`` particle types
into a single type, ``star``.  yt will then determine which fields are
accessible to this new particle type and it will add them.

.. code-block:: python

   from yt.data_objects.particle_unions import ParticleUnion

   u = ParticleUnion("star", ["halo", "disk"])
   ds.add_particle_union(u)

.. _filtering-by-location:

Filtering Fields by Spatial Location: Geometric Objects
-------------------------------------------------------

Creating geometric objects for a dataset provides a means for filtering
a field based on spatial location.  The most commonly used of these are
spheres, regions (3D prisms), ellipsoids, disks, and rays.  The ``all_data``
object which gets used throughout this documentation section is an example of
a geometric object, but it defaults to including all the data in the dataset
volume.  To see all of the geometric objects available, see
:ref:`available-objects`.

Consult the object documentation section for all of the different objects
one can use, but here is a simple example using a sphere object to filter
a dataset.  Let's filter out everything not within 10 Mpc of some random
location, say [0.2, 0.5, 0.1], in the simulation volume.  The resulting object
will only contain grid cells with centers falling inside of our defined sphere,
which may look offset based on the presence of different resolution elements
distributed throughout the dataset.

.. notebook-cell::

    import yt

    ds = yt.load("Enzo_64/DD0042/data0042")
    center = [0.20, 0.50, 0.10]

    sp = ds.sphere(center, (10, "Mpc"))
    prj = yt.ProjectionPlot(
        ds, "x", ("gas", "density"), center=center, width=(50, "Mpc"), data_source=sp
    )

    # Mark the center with a big X
    prj.annotate_marker(center, "x", plot_args={"s": 100})

    prj.show()

    slc = yt.SlicePlot(
        ds, "x", ("gas", "density"), center=center, width=(50, "Mpc"), data_source=sp
    )

    slc.show()
