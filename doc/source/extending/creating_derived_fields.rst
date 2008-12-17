Creating Derived Fields
=======================

One of the more powerful means of extending yt is through the usage of derived
fields.  These are fields that describe a value at each cell in a simulation.

Defining a New Field
--------------------

So once a new field has been conceived of, the best way to create it is to
construct a function that performs an array operation -- operating on a 
collection of data, neutral to its size, shape, and type.  (All fields should
be provided as 64-bit floats.)

A simple example of this is the pressure field, which demonstrates the ease of
this approach.

.. code-block:: python

   def _Pressure(field, data):
       return (data.pf["Gamma"] - 1.0) * \
              data["Density"] * data["ThermalEnergy"]

Note that we do a couple different things here.  We access the "Gamma"
parameter from the parameter file, we access the "Density" field and we access
the "ThermalEnergy" field.  "ThermalEnergy" is, in fact, another derived field!
("ThermalEnergy" deals with the distinction in storage of energy between dual
energy formalism and non-DEF.)  We don't do any loops, we don't do any
type-checking, we can simply multiply the three items together.

Once we've defined our function, we need to notify yt that the field is
available.  The :func:`add_field` function is the means of doing this; it has a
number of fairly specific parameters that can be passed in, but here we'll only
look at the most basic ones needed for a simple scalar baryon field.

.. code-block:: python

   add_field("Pressure", units=r"\rm{dyne}/\rm{cm}^{2}")

We feed it the name of the field, and the units.  Note that the units parameter
is a "raw" string, with some LaTeX-style formatting -- Matplotlib actually has
a MathText rendering engine, so if you include LaTeX it will be rendered
appropriately.

One very important thing to note about the call to ``add_field`` is that it
**does not** need to specify the function name **if** the function is the name
of the field prefixed with an underscore.  If it is not -- and it won't be, for
fields in different units (such as "CellMassMsun", for instance) -- then you
need to specify it with the argument ``function``.

Note one last thing about this definition; we do not unit conversion.  All of
the fields fed into the field are pre-supposed to be in CGS.  If we need no
constants applied after that, we need not include them.  (If we do, it is
better if you define a second function and use the argument
``convert_function`` to ``add_field`` to point to it.)

Some More Complicated Examples
------------------------------

But what if we want to do some more fancy stuff?  Here's an example of getting
parameters from the data object and using those to define the field;
specifically, here we obtain the ``center`` and ``height_vector`` parameters
and use those to define an angle of declination of a point with respect to a
disk.

.. code-block:: python

   def _DiskAngle(field, data):
       # We make both r_vec and h_vec into unit vectors
       center = data.get_field_parameter("center")
       r_vec = na.array([data["x"] - center[0],
                         data["y"] - center[1],
                         data["z"] - center[2]])
       r_vec = r_vec/na.sqrt((r_vec**2.0).sum(axis=0))
       h_vec = na.array(data.get_field_parameter("height_vector"))
       dp = r_vec[0,:] * h_vec[0] \
          + r_vec[1,:] * h_vec[1] \
          + r_vec[2,:] * h_vec[2]
       return na.arccos(dp)
   add_field("DiskAngle", take_log=False,
             validators=[ValidateParameter("height_vector"),
                         ValidateParameter("center")],
             display_field=False)

Note that we have added a few parameters below the main function; we specify
that we do not wish to display this field as logged, that we require both
``height_vector`` and ``center`` to be present in a given data object we wish
to calculate this for, and we say that it should not be displayed in a
drop-down box of fields to display.  This is done through the parameter
*validators*, which accepts a list of :class:`FieldValidator` objects.  These
objects define the way in which the field is generated, and when it is able to
be created.  In this case, we mandate that parameters *center* and
*height_vector* are set before creating the field.  These are set via 
:meth:`~yt.lagos.EnzoData.set_field_parameter`, which can be called on any
object that has fields.

We can also define vector fields.

.. code-block:: python

   def _SpecificAngularMomentum(field, data):
       if data.has_field_parameter("bulk_velocity"):
           bv = data.get_field_parameter("bulk_velocity")
       else: bv = na.zeros(3, dtype='float64')
       xv = data["x-velocity"] - bv[0]
       yv = data["y-velocity"] - bv[1]
       zv = data["z-velocity"] - bv[2]
       center = data.get_field_parameter('center')
       coords = na.array([data['x'],data['y'],data['z']], dtype='float64')
       new_shape = tuple([3] + [1]*(len(coords.shape)-1))
       r_vec = coords - na.reshape(center,new_shape)
       v_vec = na.array([xv,yv,zv], dtype='float64')
       return na.cross(r_vec, v_vec, axis=0)
   def _convertSpecificAngularMomentum(data):
       return data.convert("cm")
   add_field("SpecificAngularMomentum",
             convert_function=_convertSpecificAngularMomentum, vector_field=True,
             units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])

Here we define the SpecificAngularMomentum field, optionally taking a
``bulk_velocity``, and returning a vector field that needs conversion by the
function ``_convertSpecificAngularMomentum``.

Field Options
-------------

The arguments to :func:`add_field` are passed on to the constructor of
:class:`DerivedField`.  :func:`add_field` takes care of finding the arguments
`function` and `convert_function` if it can, however.

.. autoclass:: yt.lagos.DerivedField

How Do Units Work?
------------------

Everything is done under the assumption that all of the native Enzo fields that
yt knows about are converted to cgs before being handed to any processing
routines.

Which Enzo Fields Does yt Know About?
-------------------------------------

* Density
* Temperature
* Gas Energy
* Total Energy
* [xyz]-velocity
* Species fields: HI, HII, Electron, HeI, HeII, HeIII, HM, H2I, H2II, DI, DII, HDI
* Particle mass, velocity, 

