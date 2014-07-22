.. _creating-derived-fields:

Creating Derived Fields
=======================

One of the more powerful means of extending ``yt`` is through the usage of derived
fields.  These are fields that describe a value at each cell in a simulation.

Defining a New Field
--------------------

So once a new field has been conceived of, the best way to create it is to
construct a function that performs an array operation -- operating on a 
collection of data, neutral to its size, shape, and type. (All fields should
be provided as 64-bit floats.)

A simple example of this is the pressure field, which demonstrates the ease of
this approach.

.. code-block:: python

   import yt

   def _pressure(field, data):
       return (data.ds.gamma - 1.0) * \
              data["density"] * data["thermal_energy"]

Note that we do a couple different things here.  We access the "gamma"
parameter from the dataset, we access the "density" field and we access
the "thermal_energy" field.  "thermal_energy" is, in fact, another derived field!
("thermal_energy" deals with the distinction in storage of energy between dual
energy formalism and non-DEF.)  We don't do any loops, we don't do any
type-checking, we can simply multiply the three items together.

Once we've defined our function, we need to notify ``yt`` that the field is
available.  The :func:`add_field` function is the means of doing this; it has a
number of fairly specific parameters that can be passed in, but here we'll only
look at the most basic ones needed for a simple scalar baryon field.

.. code-block:: python

   yt.add_field("pressure", function=_pressure, units="dyne/cm**2")

We feed it the name of the field, the name of the function, and the
units.  Note that the units parameter is a "raw" string, in the format that ``yt`` uses
in its `symbolic units implementation <units>`_ (e.g., employing only unit names, numbers,
and mathematical operators in the string, and using ``"**"`` for exponentiation). We suggest
that you name the function that creates a derived field with the intended field name prefixed
by a single underscore, as in the ``_pressure`` example above.

:func:`add_field` can be invoked in two other ways. The first is by the function
decorator :func:`derived_field`. The following code is equivalent to the previous
example:

.. code-block:: python

   from yt import derived_field

   @derived_field(name="pressure", units="dyne/cm**2")
   def _pressure(field, data):
       return (data.ds.gamma - 1.0) * \
              data["density"] * data["thermal_energy"]

The :func:`derived_field` decorator takes the same arguments as :func:`add_field`,
and is often a more convenient shorthand in cases where you want to quickly set up
a new field.

Defining derived fields in the above fashion must be done before a dataset is loaded,
in order for the dataset to recognize it. If you want to set up a derived field after you
have loaded a dataset, or if you only want to set up a derived field for a particular
dataset, there is an :meth:`add_field` method that hangs off dataset objects. The calling
syntax is the same:

.. code-block:: python

   ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   ds.add_field("pressure", function=_pressure, units="dyne/cm**2")

If you find yourself using the same custom-defined fields over and over, you
should put them in your plugins file as described in :ref:`plugin-file`.

A More Complicated Example
--------------------------

But what if we want to do something a bit more fancy?  Here's an example of getting
parameters from the data object and using those to define the field;
specifically, here we obtain the ``center`` and ``bulk_velocity`` parameters
and use those to define a field for radial velocity (there is already a ``"radial_velocity"``
field in ``yt``, but we create this one here just as a transparent and simple example).

.. code-block:: python

   from yt.fields.api import ValidateParameter
   import numpy as np

   def _my_radial_velocity(field, data):
       if data.has_field_parameter("bulk_velocity"):
           bv = data.get_field_parameter("bulk_velocity").in_units("cm/s")
       else:
           bv = data.ds.arr(np.zeros(3), "cm/s")
       xv = data["gas","velocity_x"] - bv[0]
       yv = data["gas","velocity_y"] - bv[1]
       zv = data["gas","velocity_z"] - bv[2]
       center = data.get_field_parameter('center')
       x_hat = data["x"] - center[0]
       y_hat = data["y"] - center[1]
       z_hat = data["z"] - center[2]
       r = np.sqrt(x_hat*x_hat+y_hat*y_hat+z_hat*z_hat)
       x_hat /= r
       y_hat /= r
       z_hat /= r
       return xv*x_hat + yv*y_hat + zv*z_hat
   yt.add_field("my_radial_velocity",
                function=_my_radial_velocity,
                units="cm/s",
                take_log=False,
                validators=[ValidateParameter('center'),
                            ValidateParameter('bulk_velocity')])

Note that we have added a few parameters below the main function; we specify
that we do not wish to display this field as logged, that we require both
``bulk_velocity`` and ``center`` to be present in a given data object we wish
to calculate this for, and we say that it should not be displayed in a
drop-down box of fields to display. This is done through the parameter
*validators*, which accepts a list of :class:`FieldValidator` objects. These
objects define the way in which the field is generated, and when it is able to
be created. In this case, we mandate that parameters *center* and
*bulk_velocity* are set before creating the field. These are set via
:meth:`~yt.data_objects.data_containers.set_field_parameter`, which can 
be called on any object that has fields:

.. code-block:: python

   ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   sp = ds.sphere("max", (200.,"kpc"))
   sp.set_field_parameter("bulk_velocity", yt.YTArray([-100.,200.,300.], "km/s"))

In this case, we already know what the *center* of the sphere is, so we do not set it. Also,
note that *center* and *bulk_velocity* need to be :class:`YTArray` objects with units.

Field Options
-------------

The arguments to :func:`add_field` are passed on to the constructor of :class:`DerivedField`.
There are a number of options available, but the only mandatory ones are ``name``,
``units``, and ``function``.

   ``name``
     This is the name of the field -- how you refer to it.  For instance,
     ``pressure`` or ``magnetic_field_strength``.
   ``function``
     This is a function handle that defines the field
   ``units``
     This is a string that describes the units. Powers must be in
     python syntax (** instead of ^).
   ``display_name``
     This is a name used in the plots, for instance ``"Divergence of
     Velocity"``.  If not supplied, the ``name`` value is used.
   ``take_log``
     This is *True* or *False* and describes whether the field should be logged
     when plotted.
   ``particle_type``
     Is this field a *particle* field?
   ``validators``
     (*Advanced*) This is a list of :class:`FieldValidator` objects, for instance to mandate
     spatial data.
   ``display_field``
     (*Advanced*) Should this field appear in the dropdown box in Reason?
   ``not_in_all``
     (*Advanced*) If this is *True*, the field may not be in all the grids.
   ``output_units``
     (*Advanced*) For fields that exist on disk, which we may want to convert to other
     fields or that get aliased to themselves, we can specify a different
     desired output unit than the unit found on disk.

Units for Cosmological Datasets
-------------------------------

``yt`` has additional capabilities to handle the comoving coordinate system used
internally in cosmological simulations. Simulations that use comoving
coordinates, all length units have three other counterparts correspoding to
comoving units, scaled comoving units, and scaled proper units. In all cases
'scaled' units refer to scaling by the reduced Hubble constant - i.e. the length
unit is what it would be in a universe where Hubble's constant is 100 km/s/Mpc.  

To access these different units, yt has a common naming system. Scaled units
are denoted by appending ``h`` to the end of the unit name. Comoving units are
denoted by appending ``cm`` to the end of the unit name. If both are used, the
strings should be appended in that order: 'Mpchcm', *but not* 'Mpccmh'.

Using the parsec as an example,

``pc``
    Proper parsecs, :math:`\rm{pc}`.

``pccm``
    Comoving parsecs, :math:`\rm{pc}/(1+z)`.

``pchcm``
    Comoving parsecs normalized by the scaled hubble constant, :math:`\rm{pc}/h/(1+z)`.

``pch``
    Proper parsecs, normalized by the scaled hubble constant, :math:`\rm{pc}/h`.
