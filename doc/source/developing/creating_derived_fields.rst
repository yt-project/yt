.. _creating-derived-fields:

Creating Derived Fields
=======================

One of the more powerful means of extending yt is through the usage of derived
fields.  These are fields that describe a value at each cell in a simulation.

Defining a New Field
--------------------

Once a new field has been conceived of, the best way to create it is to
construct a function that performs an array operation -- operating on a
collection of data, neutral to its size, shape, and type.

A simple example of this is the pressure field, which demonstrates the ease of
this approach.

.. code-block:: python

   import yt


   def _pressure(field, data):
       return (
           (data.ds.gamma - 1.0)
           * data["gas", "density"]
           * data["gas", "specific_thermal_energy"]
       )

Note that we do a couple different things here.  We access the ``gamma``
parameter from the dataset, we access the ``density`` field and we access
the ``specific_thermal_energy`` field.  ``specific_thermal_energy`` is, in
fact, another derived field!  We don't do any loops, we don't do any
type-checking, we can simply multiply the three items together.

In this example, the ``density`` field will return data with units of
``g/cm**3`` and the ``specific_thermal_energy`` field will return data units of
``erg/g``, so the result will automatically have units of pressure,
``erg/cm**3``. This assumes the unit system is set to the default, which is
CGS: if a different unit system is selected, the result will be in the same
dimensions of pressure but different units. See :ref:`units` for more
information.

Once we've defined our function, we need to notify yt that the field is
available.  The :func:`add_field` function is the means of doing this; it has a
number of fairly specific parameters that can be passed in, but here we'll only
look at the most basic ones needed for a simple scalar baryon field.

.. note::

    There are two different :func:`add_field` functions.  For the differences,
    see :ref:`faq-add-field-diffs`.

.. code-block:: python

    yt.add_field(
        name=("gas", "pressure"),
        function=_pressure,
        sampling_type="local",
        units="dyne/cm**2",
    )

We feed it the name of the field, the name of the function, the sampling type,
and the units. The ``sampling_type`` keyword determines which elements are
used to make the field (i.e., grid cell or particles) and controls how volume
is calculated. It can be set to "cell" for grid/mesh fields, "particle" for
particle and SPH fields, or "local" to use the primary format of the loaded
dataset. In most cases, "local" is sufficient, but "cell" and "particle"
can be used to specify the source for datasets that have both grids and
particles. In a dataset with both grids and particles, using "cell" will
ensure a field is created with a value for every grid cell, while using
"particle" will result in a field with a value for every particle.

The units parameter is a "raw" string, in the format that yt
uses in its :ref:`symbolic units implementation <units>` (e.g., employing only
unit names, numbers, and mathematical operators in the string, and using
``"**"`` for exponentiation). For cosmological datasets and fields, see
:ref:`cosmological-units <cosmological-units>`.  We suggest that you name the function that creates
a derived field with the intended field name prefixed by a single underscore,
as in the ``_pressure`` example above.

Field definitions return array data with units. If the field function returns
data in a dimensionally equivalent unit (e.g. a ``"dyne"`` versus a ``"N"``), the
field data will be converted to the units specified in ``add_field`` before
being returned in a data object selection. If the field function returns data
with dimensions that are incompatible with units specified in ``add_field``,
you will see an error. To clear this error, you must ensure that your field
function returns data in the correct units. Often, this means applying units to
a dimensionless float or array.

If your field definition includes physical constants rather than defining a
constant as a float, you can import it from ``yt.units``
to get a predefined version of the constant with the correct units. If you know
the units your data is supposed to have ahead of time, you can also import unit
symbols like ``g`` or ``cm`` from the ``yt.units`` namespace and multiply the
return value of your field function by the appropriate combination of unit
symbols for your field's units. You can also convert floats or NumPy arrays into
:class:`~yt.units.yt_array.YTArray` or :class:`~yt.units.yt_array.YTQuantity`
instances by making use of the
:func:`~yt.data_objects.static_output.Dataset.arr` and
:func:`~yt.data_objects.static_output.Dataset.quan` convenience functions.

Lastly, if you do not know the units of your field ahead of time, you can
specify ``units='auto'`` in the call to ``add_field`` for your field.  This will
automatically determine the appropriate units based on the units of the data
returned by the field function. This is also a good way to let your derived
fields be automatically converted to the units of the unit system in your
dataset.

If ``units='auto'`` is set, it is also required to set the ``dimensions`` keyword
argument so that error-checking can be done on the derived field to make sure that
the dimensionality of the returned array and the field are the same:

.. code-block:: python

    import yt
    from yt.units import dimensions


    def _pressure(field, data):
        return (
            (data.ds.gamma - 1.0)
            * data["gas", "density"]
            * data["gas", "specific_thermal_energy"]
        )


    yt.add_field(
        ("gas", "pressure"),
        function=_pressure,
        sampling_type="local",
        units="auto",
        dimensions=dimensions.pressure,
    )

If ``dimensions`` is not set, an error will be thrown. The ``dimensions`` keyword
can be a SymPy ``symbol`` object imported from ``yt.units.dimensions``, a compound
dimension of these, or a string corresponding to one of these objects.

:func:`add_field` can be invoked in two other ways. The first is by the
function decorator :func:`derived_field`. The following code is equivalent to
the previous example:

.. code-block:: python

   from yt import derived_field


   @derived_field(name="pressure", sampling_type="cell", units="dyne/cm**2")
   def _pressure(field, data):
       return (
           (data.ds.gamma - 1.0)
           * data["gas", "density"]
           * data["gas", "specific_thermal_energy"]
       )

The :func:`derived_field` decorator takes the same arguments as
:func:`add_field`, and is often a more convenient shorthand in cases where
you want to quickly set up a new field.

Defining derived fields in the above fashion must be done before a dataset is
loaded, in order for the dataset to recognize it. If you want to set up a
derived field after you have loaded a dataset, or if you only want to set up
a derived field for a particular dataset, there is an
:func:`~yt.data_objects.static_output.Dataset.add_field` method that hangs off
dataset objects. The calling syntax is the same:

.. code-block:: python

   ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   ds.add_field(
       ("gas", "pressure"),
       function=_pressure,
       sampling_type="cell",
       units="dyne/cm**2",
   )

If you specify fields in this way, you can take advantage of the dataset's unit
system to define the units for you, so that the units will be returned in the
units of that system:

.. code-block:: python

    ds.add_field(
        ("gas", "pressure"),
        function=_pressure,
        sampling_type="cell",
        units=ds.unit_system["pressure"],
    )

Since the :class:`yt.units.unit_systems.UnitSystem` object returns a :class:`yt.units.unit_object.Unit` object when
queried, you're not limited to specifying units in terms of those already available. You can specify units for fields
using basic arithmetic if necessary:

.. code-block:: python

    ds.add_field(
        ("gas", "my_acceleration"),
        function=_my_acceleration,
        sampling_type="cell",
        units=ds.unit_system["length"] / ds.unit_system["time"] ** 2,
    )

If you find yourself using the same custom-defined fields over and over, you should put them in your plugins file as
described in :ref:`plugin-file`.

A More Complicated Example
--------------------------

But what if we want to do something a bit more fancy?  Here's an example of getting
parameters from the data object and using those to define the field;
specifically, here we obtain the ``center`` and ``bulk_velocity`` parameters
and use those to define a field for radial velocity (there is already
a ``radial_velocity`` field in yt, but we create this one here just as a
transparent and simple example).

.. code-block:: python

   import numpy as np

   from yt.fields.api import ValidateParameter


   def _my_radial_velocity(field, data):
       if data.has_field_parameter("bulk_velocity"):
           bv = data.get_field_parameter("bulk_velocity").in_units("cm/s")
       else:
           bv = data.ds.arr(np.zeros(3), "cm/s")
       xv = data["gas", "velocity_x"] - bv[0]
       yv = data["gas", "velocity_y"] - bv[1]
       zv = data["gas", "velocity_z"] - bv[2]
       center = data.get_field_parameter("center")
       x_hat = data["gas", "x"] - center[0]
       y_hat = data["gas", "y"] - center[1]
       z_hat = data["gas", "z"] - center[2]
       r = np.sqrt(x_hat * x_hat + y_hat * y_hat + z_hat * z_hat)
       x_hat /= r
       y_hat /= r
       z_hat /= r
       return xv * x_hat + yv * y_hat + zv * z_hat


   yt.add_field(
       ("gas", "my_radial_velocity"),
       function=_my_radial_velocity,
       sampling_type="cell",
       units="cm/s",
       take_log=False,
       validators=[ValidateParameter(["center", "bulk_velocity"])],
   )

Note that we have added a few optional arguments to ``yt.add_field``; we specify
that we do not wish to display this field as logged, that we require both the
``bulk_velocity`` and ``center`` field parameters to be present in a given data
object we wish to calculate this for, and we say that it should not be displayed
in a drop-down box of fields to display. This is done through the parameter
*validators*, which accepts a list of
:class:`~yt.fields.derived_field.FieldValidator` objects. These objects define
the way in which the field is generated, and when it is able to be created. In
this case, we mandate that parameters ``center`` and ``bulk_velocity`` are set
before creating the field. These are set via
:meth:`~yt.data_objects.data_containers.set_field_parameter`, which can be
called on any object that has fields:

.. code-block:: python

   ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   sp = ds.sphere("max", (200.0, "kpc"))
   sp.set_field_parameter("bulk_velocity", yt.YTArray([-100.0, 200.0, 300.0], "km/s"))

In this case, we already know what the ``center`` of the sphere is, so we do
not set it. Also, note that ``center`` and ``bulk_velocity`` need to be
:class:`~yt.units.yt_array.YTArray` objects with units.

If you are writing a derived field that uses a field parameter that changes the
behavior of the field depending on the value of the field parameter, you can
make yt test to make sure the field handles all possible values for the field
parameter using a special form of the ``ValidateParameter`` field validator. In
particular, ``ValidateParameter`` supports an optional second argument, which
takes a dictionary mapping from parameter names to parameter values that
you would like yt to test. This is useful when a field will select different
fields to access based on the value of a field parameter. This option allows you
to force yt to select *all* needed dependent fields for your derived field
definition at field detection time. This can avoid errors related to missing fields.

For example, let's write a field that depends on a field parameter named ``'axis'``:

.. code-block:: python

   def my_axis_field(field, data):
       axis = data.get_field_parameter("axis")
       if axis == 0:
           return data["gas", "velocity_x"]
       elif axis == 1:
           return data["gas", "velocity_y"]
       elif axis == 2:
           return data["gas", "velocity_z"]
       else:
           raise ValueError


   ds.add_field(
       "my_axis_field",
       function=my_axis_field,
       units="cm/s",
       validators=[ValidateParameter("axis", {"axis": [0, 1, 2]})],
   )

In this example, we've told yt's field system that the data object we are
querying ``my_axis_field`` must have the ``axis`` field parameter set. In
addition, it forces yt to recognize that this field might depend on any one of
``x-velocity``, ``y-velocity``, or ``z-velocity``. By specifying that ``axis``
might be 0, 1, or 2 in the ``ValidataParameter`` call, this ensures that this
field will only be valid and available for datasets that have all three fields
available.

Other examples for creating derived fields can be found in the cookbook recipe
:ref:`cookbook-simple-derived-fields`.

.. _derived-field-options:

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
     This is a string that describes the units, or a query to a UnitSystem
     object, e.g. ``ds.unit_system["energy"]``. Powers must be in Python syntax (``**``
     instead of ``^``). Alternatively, it may be set to ``"auto"`` to have the units
     determined automatically. In this case, the ``dimensions`` keyword must be set to the
     correct dimensions of the field.
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
``force_override``
     (*Advanced*) Overrides the definition of an old field if a field with the
     same name has already been defined.
``dimensions``
     Set this if ``units="auto"``. Can be either a string or a dimension object from
     ``yt.units.dimensions``.

Debugging a Derived Field
-------------------------

If your derived field is not behaving as you would like, you can insert a call
to ``data._debug()`` to spawn an interactive interpreter whenever that line is
reached.  Note that this is slightly different from calling
``pdb.set_trace()``, as it will *only* trigger when the derived field is being
called on an actual data object, rather than during the field detection phase.
The starting position will be one function lower in the stack than you are
likely interested in, but you can either step through back to the derived field
function, or simply type ``u`` to go up a level in the stack.

For instance, if you had defined this derived field:

.. code-block:: python

   @yt.derived_field(name=("gas", "funthings"))
   def funthings(field, data):
       return data["sillythings"] + data["humorousthings"] ** 2.0

And you wanted to debug it, you could do:

.. code-block:: python

   @yt.derived_field(name=("gas", "funthings"))
   def funthings(field, data):
       data._debug()
       return data["sillythings"] + data["humorousthings"] ** 2.0

And now, when that derived field is actually used, you will be placed into a
debugger.
