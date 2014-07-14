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
       return (data.pf.gamma - 1.0) * \
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
and mathematical operators in the string, and using ``"**"`` for exponentiation).

.. One very important thing to note about the call to ``add_field`` is
.. that it **does not** need to specify the function name **if** the
.. function is the name of the field prefixed with an underscore.  If it
.. is not -- and it won't be for fields in different units (such as
.. "cell_mass") -- then you need to specify it with the argument
.. ``function``.

We suggest that you name the function that creates a derived field
with the intended field name prefixed by a single underscore, as in
the ``_pressure`` example above.

:func:`add_field` can be invoked in two other ways. The first is by the function
decorator :func:`derived_field`. The following code is equivalent to the previous
example:

.. code-block:: python

   from yt import derived_field

   @derived_field(name="pressure", units="dyne/cm**2")
   def _pressure(field, data):
       return (data.pf.gamma - 1.0) * \
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

Some More Complicated Examples
------------------------------

But what if we want to do some more fancy stuff?  Here's an example of getting
parameters from the data object and using those to define the field;
specifically, here we obtain the ``center`` and ``height_vector`` parameters
and use those to define an angle of declination of a point with respect to a
disk.

.. code-block:: python

   def _disk_angle(field, data):
       # We make both r_vec and h_vec into unit vectors
       center = data.get_field_parameter("center")
       r_vec = np.array([data["x"] - center[0],
                         data["y"] - center[1],
                         data["z"] - center[2]])
       r_vec = r_vec/np.sqrt((r_vec**2.0).sum(axis=0))
       h_vec = np.array(data.get_field_parameter("height_vector"))
       dp = r_vec[0,:] * h_vec[0] \
          + r_vec[1,:] * h_vec[1] \
          + r_vec[2,:] * h_vec[2]
       return np.arccos(dp)
   yt.add_field("disk_angle", take_log=False,
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
:meth:`~yt.data_objects.data_containers.set_field_parameter`, which can 
be called on any object that has fields.

We can also define vector fields.

.. code-block:: python

   def _specific_angular_momentum(field, data):
       if data.has_field_parameter("bulk_velocity"):
           bv = data.get_field_parameter("bulk_velocity")
       else:
           bv = np.zeros(3, dtype='float64')
       xv = data["velocity_x"] - bv[0]
       yv = data["velocity_y"] - bv[1]
       zv = data["velocity_z"] - bv[2]
       center = data.get_field_parameter('center')
       coords = np.array([data['x'],data['y'],data['z']], dtype='float64')
       new_shape = tuple([3] + [1]*(len(coords.shape)-1))
       r_vec = coords - np.reshape(center,new_shape)
       v_vec = np.array([xv,yv,zv], dtype='float64')
       return np.cross(r_vec, v_vec, axis=0)
   add_field("specific_angular_momentum",
             vector_field=True, units="cm**2/s",
             validators=[ValidateParameter('center')])

Here we define the ``specific_angular_momentum`` field, optionally taking a
``bulk_velocity``, and returning a vector field.

It is also possible to define fields that depend on spatial derivatives of 
other fields.  Calculating the derivative for a single grid cell requires 
information about neighboring grid cells.  Therefore, properly calculating 
a derivative for a cell on the edge of the grid will require cell values from 
neighboring grids.  Below is an example of a field that is the divergence of the 
velocity.

.. code-block:: python

    def _DivV(field, data):
        # We need to set up stencils
        if data.pf["HydroMethod"] == 2:
            sl_left = slice(None,-2,None)
            sl_right = slice(1,-1,None)
            div_fac = 1.0
        else:
            sl_left = slice(None,-2,None)
            sl_right = slice(2,None,None)
            div_fac = 2.0
        ds = div_fac * data['dx'].flat[0]
        f  = data["velocity_x"][sl_right,1:-1,1:-1]/ds
        f -= data["velocity_x"][sl_left ,1:-1,1:-1]/ds
        if data.pf.dimensionality > 1:
            ds = div_fac * data['dy'].flat[0]
            f += data["velocity_y"][1:-1,sl_right,1:-1]/ds
            f -= data["velocity_y"][1:-1,sl_left ,1:-1]/ds
        if data.pf.dimensionality > 2:
            ds = div_fac * data['dz'].flat[0]
            f += data["velocity_z"][1:-1,1:-1,sl_right]/ds
            f -= data["velocity_z"][1:-1,1:-1,sl_left ]/ds
        new_field = np.zeros(data["velocity_x"].shape, dtype='float64')
        new_field[1:-1,1:-1,1:-1] = f
        return new_field
    def _convertDivV(data):
        return data.convert("cm")**-1.0
    add_field("DivV", function=_DivV,
               validators=[ValidateSpatial(ghost_zones=1,
	                   fields=["velocity_x","velocity_y","velocity_z"])],
              units=r"\rm{s}^{-1}", take_log=False,
              convert_function=_convertDivV)

Note that *slice* is simply a native Python object used for taking slices of 
arrays or lists.  Another :class:`FieldValidator` object, ``ValidateSpatial`` 
is given in the list of *validators* in the call to ``add_field`` with 
*ghost_zones* = 1, specifying that the original grid be padded with one additional 
cell from the neighboring grids on all sides.  The *fields* keyword simply 
mandates that the listed fields be present.  With one ghost zone added to all sides 
of the grid, the data fields (data["velocity_x"], data["velocity_y"], and 
data["velocity_z"]) will have a shape of (NX+2, NY+2, NZ+2) inside of this function, 
where the original grid has dimension (NX, NY, NZ).  However, when the final field 
data is returned, the ghost zones will be removed and the shape will again be 
(NX, NY, NZ).

.. _derived-field-options:

Saving Derived Fields
---------------------

Complex fields can be time-consuming to generate, especially on large datasets.
To mitigate this, ``yt`` provides a mechanism for saving fields to a backup file
using the Grid Data Format. The next time you start yt, it will check this file
and your field will be treated as native if present. 

The code below creates a new derived field called "dinosaurs" and saves it to disk:

.. code-block:: python

    import yt
    from yt.utilities.grid_data_format import writer
    import numpy as np

    def _dinosaurs(field, data) :
        return data["temperature"]*np.sqrt(data["density"])
    yt.add_field("dinosaurs", units="K*sqrt(g)/sqrt(cm**3)")

    ds = yt.load('GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100')
    writer.save_field(ds, "dinosaurs")

This creates a "_backup.gdf" file next to your datadump. If you load up the dataset again:

.. code-block:: python

    import yt

    ds = yt.load('GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100')
    dd = ds.all_data()
    print dd["entropy"]

you can work with the field exactly as before, without having to recompute it.

Field Options
-------------

The arguments to :func:`add_field` are passed on to the constructor of
:class:`DerivedField`.  :func:`add_field` takes care of finding the arguments
`function` and `convert_function` if it can, however.  There are a number of
options available, but the only mandatory ones are ``name``, ``units``, and possibly
``function``.

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
comoving units, scaled comoving units, and scaled proper units.  In all cases
'scaled' units refer to scaling by the reduced Hubble constant - i.e. the length
unit is what it would be in a universe where Hubble's constant is 100 km/s/Mpc.  

To access these different units, yt has a common naming system.  Scaled units
are denoted by appending ``h`` to the end of the unit name.  Comoving units are
denoted by appending ``cm`` to the end of the unit name.  If both are used, the
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
