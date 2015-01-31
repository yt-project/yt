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

Because `mesh fields` are internally different from `particle fields`,
there are different ways of filtering each type as indicated below;
however, filtering fields by spatial location (i.e. geometric
objects) will apply to both types equally.

.. _filtering-mesh:

Filtering Mesh Fields
----------------------

Mesh fields can be filtered by two methods: cut region objects 
(:class:`~yt.data_objects.selection_data_containers.YTCutRegionBase`) 
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
    bigger_than_two = (a > 2)
    print "Original Array: a = \n%s" % a
    print "Boolean Mask: bigger_than_two = \n%s" % bigger_than_two
    print "Masked Array: a[bigger_than_two] = \n%s" % a[bigger_than_two]

Similarly, if you've created a yt data object (e.g. a region, a sphere), you 
can examine its field values as a NumPy array by simply indexing it with the 
field name.  Thus, it too can be masked using a NumPy boolean mask.  Let's
set a simple mask based on the contents of one of our fields.

.. notebook-cell::

    import yt
    ds = yt.load('Enzo_64/DD0042/data0042')
    ad = ds.all_data()
    hot = ad["temperature"].in_units('K') > 1e6
    print 'Temperature of all data: ad["temperature"] = \n%s' % ad["temperature"]
    print "Boolean Mask: hot = \n%s" % hot
    print 'Temperature of "hot" data: ad["temperature"][hot] = \n%s' % \
          ad['temperature'][hot]

This was a simple example, but one can make the conditionals that define
a boolean mask have multiple parts, and one can stack masks together to
make very complex cuts on one's data.  Once the data is filtered, it can be
used if you simply need to access the NumPy arrays:

.. notebook-cell::

    import yt
    ds = yt.load('Enzo_64/DD0042/data0042')
    ad = ds.all_data()
    overpressure_and_fast = (ad["pressure"] > 1e-14) & (ad["velocity_magnitude"].in_units('km/s') > 1e2)
    print 'Density of all data: ad["density"] = \n%s' % ad['density']
    print 'Density of "overpressure and fast" data: ad["density"][overpressure_and_fast] = \n%s' % \
           ad['density'][overpressure_and_fast]

.. _cut-regions:

Cut Regions
^^^^^^^^^^^

Cut regions are a more general solution to filtering mesh fields.  The output
of a cut region is an entirely new data object, which can be treated like any
other data object to generate images, examine its values, etc.

.. notebook:: mesh_filter.ipynb

Cut regions can also operator on particle fields, but a single cut region object
cannot operate on both particle fields and mesh fields at the same time.

.. _filtering-particles:

Filtering Particle Fields
-------------------------

Particle filters create new particle fields based on the manipulation and 
cuts on existing particle fields.  You can apply cuts to them to effectively
mask out everything except the particles with which you are concerned.

Creating a particle filter takes a few steps.  You must first define a 
function which accepts a geometric object (e.g. all_data, sphere, etc.)
as its argument.  It uses the fields and information in this geometric
object in order to produce some sort of conditional mask that is then returned.
Here is the function to filter only the particles with `particle_type` (i.e. 
field = `('all', 'particle_type')` equal to 2. (This is the case with
Enzo star particles.)

.. code-block:: python

    def Stars(pfilter, data):
        filter = data[("all", "particle_type")] == 2
        return filter

The particle_filter must now be defined to incorporate this function.  It takes
a few arguments: a name for the filter, our filter function, and the fields
that it requires in a dataset in order to work (in this case, it requires
the ('all', 'particle_type') field.

.. code-block:: python

    from yt.data_objects.particle_filters import add_particle_filter
    add_particle_filter("stars", function=Stars, filtered_type='all', requires=["particle_type"])

And lastly, the filter must be applied to our dataset of choice.  Note that this 
filter can be added to as many datasets as we wish.  It will only actually
create new filtered fields if the dataset has the required fields, though.

.. code-block:: python

    import yt
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    ds.add_particle_filter('stars')

And that's it!  We can now access all of the ('stars', field) fields from 
our dataset `ds` and treat them as any other particle field.  In addition,
it created some `deposit` fields, where the particles were deposited on to
the grid as mesh fields.

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

   from yt.data_objects.particle_unions import \
       ParticleUnion

   u = ParticleUnion("star", ["halo", "disk"])
   ds.add_particle_union(u)

.. _filtering-by-location:

Filtering Fields by Spatial Location: Geometric Objects
-------------------------------------------------------

Creating geometric objects for a dataset provides a means for filtering
a field based on spatial location.  The most commonly used of these are
spheres, regions (3D prisms), ellipsoids, disks, and rays.  The `all_data`
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
    ds = yt.load('Enzo_64/DD0042/data0042')
    center = [0.20, 0.50, 0.10]

    sp = ds.sphere(center, (10, 'Mpc'))
    prj = yt.ProjectionPlot(ds, "x", "density", center=center, width=(50, "Mpc"),
                            data_source=sp)

    # Mark the center with a big X
    prj.annotate_marker(center, 'x', plot_args={'s':100})

    prj.show()

    slc = yt.SlicePlot(ds, "x", "density", center=center, width=(50, "Mpc"),
                       data_source=sp)

    slc.show()
