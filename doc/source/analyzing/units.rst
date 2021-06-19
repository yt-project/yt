.. _units:

Symbolic Units
==============

This section describes yt's symbolic unit capabilities. This is provided as
quick introduction for those who are already familiar with yt but want to learn
more about the unit system.  Please see :ref:`analyzing` and :ref:`visualizing`
for more detail about querying, analyzing, and visualizing data in yt.

Originally the unit system was a part of yt proper but since the yt 4.0 release,
the unit system has been split off into `its own library
<https://github.com/yt-project/unyt>`_, ``unyt``.

For a detailed discussion of how to use ``unyt``, we suggest taking a look at
the unyt documentation available at https://unyt.readthedocs.io/, however yt
adds additional capabilities above and beyond what is provided by ``unyt``
alone, we describe those capabilities below.

Selecting data from a data object
---------------------------------

The data returned by yt will have units attached to it. For example, let's query
a data object for the ``('gas', 'density')`` field:

    >>> import yt
    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> dd = ds.all_data()
    >>> dd['gas', 'density']
    unyt_array([4.92775113e-31, 4.94005233e-31, 4.93824694e-31, ...,
                1.12879234e-25, 1.59561490e-25, 1.09824903e-24], 'g/cm**3')

We can see how we get back a ``unyt_array`` instance. A ``unyt_array`` is a
subclass of NumPy's NDarray type that has units attached to it:

    >>> dd['gas', 'density'].units
    g/cm**3

It is straightforward to convert data to different units:

    >>> dd['gas', 'density'].to('Msun/kpc**3')
    unyt_array([7.28103608e+00, 7.29921182e+00, 7.29654424e+00, ...,
               1.66785569e+06, 2.35761291e+06, 1.62272618e+07], 'Msun/kpc**3')

For more details about working with ``unyt_array``, see the `the documentation
<https://unyt.readthedocs.io>`__ for ``unyt``.

Applying Units to Data
----------------------

A ``unyt_array`` can be created from a list, tuple, or NumPy array using
multiplication with a ``Unit`` object. For convenience, each yt dataset has a
``units`` attribute one can use to obtain unit objects for this purpose:

   >>> data = np.random.random((100, 100))
   >>> data_with_units = data * ds.units.gram

All units known to the dataset will be available via ``ds.units``, including
code units and comoving units.

Derived Field Units
-------------------

Special care often needs to be taken to ensure the result of a derived field
will come out in the correct units. The yt unit system will double-check for you
to make sure you are not accidentally making a unit conversion mistake. To see
what that means in practice, let's define a derived field corresponding to the
square root of the gas density:

    >>> import yt
    >>> import numpy as np

    >>> def root_density(field, data):
    ...     return np.sqrt(data['gas', 'density'])

    >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

    >>> ds.add_field(("gas", "root_density"), units="(g/cm**3)**(1/2)",
    ...              function=root_density, sampling_type='cell')

    >>> ad = ds.all_data()
    >>> ad['gas', 'root_density']
    unyt_array([7.01979425e-16, 7.02855059e-16, 7.02726614e-16, ...,
                3.35975050e-13, 3.99451486e-13, 1.04797377e-12], 'sqrt(g)/cm**(3/2)')

No special unit logic needs to happen inside of the function: the result of
``np.sqrt`` will have the correct units:

    >>> np.sqrt(ad['gas', 'density'])
    unyt_array([7.01979425e-16, 7.02855059e-16, 7.02726614e-16, ...,
                3.35975050e-13, 3.99451486e-13, 1.04797377e-12], 'sqrt(g)/cm**(3/2)')

One could also specify any other units that have dimensions of square root of
density and yt would automatically convert the return value of the field
function to the specified units. An error would be raised if the units are not
dimensionally equivalent to the return value of the field function.

Code Units
----------

All yt datasets are associated with a "code" unit system that corresponds to
whatever unit system the data is represented in on-disk. Let's take a look at
the data in an Enzo simulation, specifically the ``("enzo", "Density")`` field:

    >>> import yt
    >>> ds = yt.load('Enzo_64/DD0043/data0043')
    >>> ad = ds.all_data()
    >>> ad["enzo", "Density"]
    unyt_array([6.74992726e-02, 6.12111635e-02, 8.92988636e-02, ...,
                9.09875931e+01, 5.66932465e+01, 4.27780263e+01], 'code_mass/code_length**3')

we see we get back data from yt in units of ``code_mass/code_length**3``. This
is the density unit formed out of the base units of mass and length in the
internal unit system in the simulation. We can see the values of these units by
looking at the ``length_unit`` and ``mass_unit`` attributes of the dataset
object:

    >>> ds.length_unit
    unyt_quantity(128, 'Mpccm/h')
    >>> ds.mass_unit
    unyt_quantity(4.89045159e+50, 'g')

And we can see that both of these have values of 1 in the code unit system.

    >>> ds.length_unit.to('code_length')
    unyt_quantity(1., 'code_length')
    >>> ds.mass_unit.to('code_mass')
    unyt_quantity(1., 'code_mass')

In addition to ``length_unit`` and ``mass_unit``, there are also ``time_unit``,
``velocity_unit``, and ``magnetic_unit`` attributes for this dataset. Some
frontends also define a ``density_unit``, ``pressure_unit``,
``temperature_unit``, and ``specific_energy`` attribute. If these are not defined
then the corresponding unit is calculated from the base length, mass, and time unit.
Each of these attributes corresponds to a unit in the code unit system:

    >>> [un for un in dir(ds.units) if un.startswith('code')]
    ['code_density',
     'code_length',
     'code_magnetic',
     'code_mass',
     'code_metallicity',
     'code_pressure',
     'code_specific_energy',
     'code_temperature',
     'code_time',
     'code_velocity']

You can use these unit names to convert arbitrary data into a dataset's code
unit system:

    >>> u = ds.units
    >>> data = 10**-30 * u.g / u.cm**3
    >>> data.to('code_density')
    unyt_quantity(0.36217187, 'code_density')

Note how in this example we used ``ds.units`` instead of the top-level ``unyt``
namespace or ``yt.units``. This is because the units from ``ds.units`` know
about the dataset's code unit system and can convert data into it. Unit objects
from ``unyt`` or ``yt.units`` will not know about any particular dataset's unit
system.


.. _cosmological-units:

Comoving units for Cosmological Simulations
-------------------------------------------

The length unit of the dataset I used above uses a cosmological unit:

    >>> print(ds.length_unit)
    128 Mpccm/h

In English, this says that the length unit is 128 megaparsecs in the comoving
frame, scaled as if the hubble constant were 100 km/s/Mpc. Although :math:`h`
isn't really a unit, yt treats it as one for the purposes of the unit system.

As an aside, `Darren Croton's research note <https://arxiv.org/abs/1308.4150>`_
on the history, use, and interpretation of :math:`h` as it appears in the
astronomical literature is pretty much required reading for anyone who has to
deal with factors of :math:`h` every now and then.

In yt, comoving length unit symbols are named following the pattern ``< length
unit >cm``, i.e. ``pccm`` for comoving parsec or ``mcm`` for a comoving
meter. A comoving length unit is different from the normal length unit by a
factor of :math:`(1+z)`:

    >>> u = ds.units
    >>> print((1*u.Mpccm)/(1*u.Mpc))
    0.9986088499304777 dimensionless
    >>> 1 / (1 + ds.current_redshift)
    0.9986088499304776

As we saw before, h is treated like any other unit symbol. It has dimensionless
units, just like a scalar:

    >>> (1*u.Mpc)/(1*u.Mpc/u.h)
    unyt_quantity(0.71, '(dimensionless)')
    >>> ds.hubble_constant
    0.71

Using parsec as an example,

  * ``pc``
    Proper parsecs, :math:`\rm{pc}`.

  * ``pccm``
    Comoving parsecs, :math:`\rm{pc}/(1+z)`.

  * ``pccm/h``
    Comoving parsecs normalized by the scaled hubble constant, :math:`\rm{pc}/h/(1+z)`.

  * ``pc/h``
    Proper parsecs, normalized by the scaled hubble constant, :math:`\rm{pc}/h`.

Overriding Code Unit Definitions
--------------------------------

On occasion, you might have a dataset for a supported frontend that does not
have the conversions to code units accessible or you may want to change them
outright. ``yt`` provides a mechanism so that one may provide their own code
unit definitions to ``yt.load``, which override the default rules for a given
frontend for defining code units.

This is provided through the ``units_override`` argument to ``yt.load``. We'll
use an example of an Athena dataset. First, a call to ``yt.load`` without
``units_override``:

    >>> ds = yt.load("MHDSloshing/virgo_low_res.0054.vtk")
    >>> ds.length_unit
    unyt_quantity(1., 'cm')
    >>> ds.mass_unit
    unyt_quantity(1., 'g')
    >>> ds.time_unit
    unyt_quantity(1., 's')
    >>> sp1 = ds1.sphere("c", (0.1, "unitary"))
    >>> print(sp1["gas", "density"])
    [0.05134981 0.05134912 0.05109047 ... 0.14608461 0.14489453 0.14385277] g/cm**3

This particular simulation is of a galaxy cluster merger so these density values
are way, way too high. This is happening because Athena does not encode any
information about the unit system used in the simulation or the output data, so
yt cannot infer that information and must make an educated guess. In this case
it incorrectly assumes the data are in CGS units.

However, we know *a priori* what the unit system *should* be, and we can supply
a ``units_override`` dictionary to ``yt.load`` to override the incorrect
assumptions yt is making about this dataset. Let's define:

    >>> units_override = {"length_unit": (1.0, "Mpc"),
    ...                   "time_unit": (1.0, "Myr"),
    ...                   "mass_unit": (1.0e14, "Msun")}

The ``units_override`` dictionary can take the following keys:

    * ``length_unit``
    * ``time_unit``
    * ``mass_unit``
    * ``magnetic_unit``
    * ``temperature_unit``

and the associated values can be ``(value, "unit")`` tuples, ``unyt_quantity``
instances, or floats (in the latter case they are assumed to have the
corresponding cgs unit). Now let's reload the dataset using our
``units_override`` dict:

    >>> ds = yt.load("MHDSloshing/virgo_low_res.0054.vtk",
    ...              units_override=units_override)
    >>> sp = ds.sphere("c",(0.1,"unitary"))
    >>> print(sp["gas", "density"])
    [3.47531683e-28 3.47527018e-28 3.45776515e-28 ... 9.88689766e-28
     9.80635384e-28 9.73584863e-28] g/cm**3

and we see how the data now have much more sensible values for a galaxy cluster
merge simulation.

Comparing Units From Different Simulations
------------------------------------------

The code units from different simulations will have different conversions to
physical coordinates. This can get confusing when working with data from more
than one simulation or from a single simulation where the units change with
time.

As an example, let's load up two enzo datasets from different redshifts in the
same cosmology simulation, one from high redshift:

    >>> ds1 = yt.load('Enzo_64/DD0002/data0002')
    >>> ds1.current_redshift
    7.8843748886903
    >>> ds1.length_unit
    unyt_quantity(128, 'Mpccm/h')
    >>> ds1.length_unit.in_cgs()
    unyt_quantity(6.26145538e+25, 'cm')

And another from low redshift:

    >>> ds2 = yt.load('Enzo_64/DD0043/data0043')
    >>> ds2.current_redshift
    0.0013930880640796
    >>> ds2.length_unit
    unyt_quantity(128, 'Mpccm/h')
    >>> ds2.length_unit.in_cgs()
    unyt_quantity(5.55517285e+26, 'cm')

Now despite the fact that ``'Mpccm/h'`` means different things for the two
datasets, it's still a well-defined operation to take the ratio of the two
length units:

    >>> ds2.length_unit / ds1.length_unit
    unyt_quantity(8.87201539, '(dimensionless)')

Because code units and comoving units are defined relative to a physical unit
system, ``unyt`` is able to give the correct answer here. So long as the result
comes out dimensionless or in a physical unit then the answer will be
well-defined. However, if we want the answer to come out in the internal units
of one particular dataset, additional care must be taken. For an example where
this might be an issue, let's try to compute the sum of two comoving distances
from each simulation:

    >>> d1 = 12 * ds1.units.Mpccm
    >>> d2 = 12 * ds2.units.Mpccm
    >>> d1 + d2
    unyt_quantity(118.46418468, 'Mpccm')
    >>> d2 + d1
    unyt_quantity(13.35256754, 'Mpccm')

So this is definitely weird - addition appears to not be associative anymore!
However, both answers are correct, the confusion is arising because ``"Mpccm"``
is ambiguous in these expressions. In situations like this, ``unyt`` will use
the definition for units from the leftmost term in an expression, so the first
example is returning data in high-redshift comoving megaparsecs, while the
second example returns data in low-redshift comoving megaparsecs.

Wherever possible it's best to do calculations in physical units when working
with more than one dataset. If you need to use comoving units or code units then
extra care must be taken in your code to avoid ambiguity.
