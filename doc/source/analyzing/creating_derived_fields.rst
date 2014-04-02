.. _creating-derived-fields:

Creating Derived Fields
=======================

One of the more powerful means of extending ``yt`` is through the usage of derived
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
              data["density"] * data["thermal_energy"]

Note that we do a couple different things here.  We access the "Gamma"
parameter from the parameter file, we access the "density" field and we access
the "thermal_energy" field.  "thermal_energy" is, in fact, another derived field!
("thermal_energy" deals with the distinction in storage of energy between dual
energy formalism and non-DEF.)  We don't do any loops, we don't do any
type-checking, we can simply multiply the three items together.

Once we've defined our function, we need to notify ``yt`` that the field is
available.  The :func:`add_field` function is the means of doing this; it has a
number of fairly specific parameters that can be passed in, but here we'll only
look at the most basic ones needed for a simple scalar baryon field.

.. code-block:: python

   add_field("pressure", function=_Pressure, units=r"\rm{dyne}/\rm{cm}^{2}")

We feed it the name of the field, the name of the function, and the
units.  Note that the units parameter is a "raw" string, with some
LaTeX-style formatting -- Matplotlib actually has a MathText rendering
engine, so if you include LaTeX it will be rendered appropriately.

.. One very important thing to note about the call to ``add_field`` is
.. that it **does not** need to specify the function name **if** the
.. function is the name of the field prefixed with an underscore.  If it
.. is not -- and it won't be for fields in different units (such as
.. "cell_mass") -- then you need to specify it with the argument
.. ``function``.

We suggest that you name the function that creates a derived field
with the intended field name prefixed by a single underscore, as in
the ``_Pressure`` example above.

If you find yourself using the same custom-defined fields over and over, you
should put them in your plugins file as described in :ref:`plugin-file`.

.. _conversion-factors:

Conversion Factors
~~~~~~~~~~~~~~~~~~

When creating a derived field, ``yt`` does not by default do unit
conversion.  All of the fields fed into the field are pre-supposed to
be in CGS.  If the field does not need any constants applied after
that, you are done. If it does, you should define a second function
that applies the proper multiple in order to return the desired units
and use the argument ``convert_function`` to ``add_field`` to point to
it.  

The argument that you pass to ``convert_function`` will be dependent on 
what fields are input into your derived field, and in what form they
are passed from their native format.  For enzo fields, nearly all the
native on-disk fields are in CGS units already (except for ``dx``, ``dy``,
and ``dz`` fields), so you typically only need to convert for 
off-standard fields taking into account where those fields are 
used in the final output derived field.  For other codes, it can vary.

You can check to see the units associated with any field in a dataset
from any code by using the ``_units`` attribute.  Here is an example 
with one of our sample FLASH datasets available publicly at 
http://yt-project.org/data :

.. code-block:: python

   >>> from yt.mods import *
   >>> pf = load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
   >>> pf.field_list
   ['dens', 'temp', 'pres', 'gpot', 'divb', 'velx', 'vely', 'velz', 'magx', 'magy', 'magz', 'magp']
   >>> pf.field_info['dens']._units
   '\\rm{g}/\\rm{cm}^{3}'
   >>> pf.field_info['temp']._units
   '\\rm{K}'
   >>> pf.field_info['velx']._units
   '\\rm{cm}/\\rm{s}'

Thus if you were using any of these fields as input to your derived field, you 
wouldn't have to worry about unit conversion because they're already in CGS.

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
       r_vec = np.array([data["x"] - center[0],
                         data["y"] - center[1],
                         data["z"] - center[2]])
       r_vec = r_vec/np.sqrt((r_vec**2.0).sum(axis=0))
       h_vec = np.array(data.get_field_parameter("height_vector"))
       dp = r_vec[0,:] * h_vec[0] \
          + r_vec[1,:] * h_vec[1] \
          + r_vec[2,:] * h_vec[2]
       return np.arccos(dp)
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
:meth:`~yt.data_objects.data_containers.set_field_parameter`, which can 
be called on any object that has fields.

We can also define vector fields.

.. code-block:: python

   def _SpecificAngularMomentum(field, data):
       if data.has_field_parameter("bulk_velocity"):
           bv = data.get_field_parameter("bulk_velocity")
       else: bv = np.zeros(3, dtype='float64')
       xv = data["velocity_x"] - bv[0]
       yv = data["velocity_y"] - bv[1]
       zv = data["velocity_z"] - bv[2]
       center = data.get_field_parameter('center')
       coords = np.array([data['x'],data['y'],data['z']], dtype='float64')
       new_shape = tuple([3] + [1]*(len(coords.shape)-1))
       r_vec = coords - np.reshape(center,new_shape)
       v_vec = np.array([xv,yv,zv], dtype='float64')
       return np.cross(r_vec, v_vec, axis=0)
   def _convertSpecificAngularMomentum(data):
       return data.convert("cm")
   add_field("SpecificAngularMomentum",
             convert_function=_convertSpecificAngularMomentum, vector_field=True,
             units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])

Here we define the SpecificAngularMomentum field, optionally taking a
``bulk_velocity``, and returning a vector field that needs conversion by the
function ``_convertSpecificAngularMomentum``.

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

The code below creates a new derived field called "Entr" and saves it to disk:

.. code-block:: python

    from yt.mods import *
    from yt.utilities.grid_data_format import writer

    def _Entropy(field, data) :
        return data["temperature"]*data["density"]**(-2./3.)
    add_field("Entr", function=_Entropy)

    pf = load('GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100')
    writer.save_field(pf, "Entr")

This creates a "_backup.gdf" file next to your datadump. If you load up the dataset again:

.. code-block:: python

    from yt.mods import *

    pf = load('GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100')
    data = pf.h.all_data()
    print data["Entr"]

you can work with the field exactly as before, without having to recompute it.

Field Options
-------------

The arguments to :func:`add_field` are passed on to the constructor of
:class:`DerivedField`.  :func:`add_field` takes care of finding the arguments
`function` and `convert_function` if it can, however.  There are a number of
options available, but the only mandatory ones are ``name`` and possibly
``function``.

   ``name``
     This is the name of the field -- how you refer to it.  For instance,
     ``Pressure`` or ``H2I_Fraction``.
   ``function``
     This is a function handle that defines the field
   ``convert_function``
     This is the function that converts the field to CGS.  All inputs to this
     function are mandated to already *be* in CGS.
   ``units``
     This is a mathtext (LaTeX-like) string that describes the units.
   ``projected_units``
     This is a mathtext (LaTeX-like) string that describes the units if the
     field has been projected without a weighting.
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
   ``vector_field``
     (*Advanced*) Is this field more than one value per cell?
   ``display_field``
     (*Advanced*) Should this field appear in the dropdown box in Reason?
   ``not_in_all``
     (*Advanced*) If this is *True*, the field may not be in all the grids.

How Do Units Work?
------------------

The best way to understand yt's unit system is to keep in mind that ``yt`` is really
handling *two* unit systems: the internal unit system of the dataset and the
physical (usually CGS) unit system.  For simulation codes like FLASH and ORION
that do all computations in CGS units internally, these two unit systems are the
same.  Most other codes do their calculations in a non-dimensionalized unit
system chosen so that most primitive variables are as close to unity as
possible.  ``yt`` allows data access both in code units and physical units by
providing a set of standard yt fields defined by all frontends.

When a dataset is loaded, ``yt`` reads the conversion factors necessary convert the
data to CGS units from the datafile itself or from a dictionary passed to the
``load`` command.  Raw on-disk fields are presented to the user via the string
names used in the dataset.  For a full enumeration of the known field names for
each of the different frontends, see the :ref:`field-list`. In general, no
conversion factors are applied to on-disk fields.

To access data in physical CGS units, yt recognizes a number of 'universal'
field names.  All primitive fields (density, pressure, magnetic field strength,
etc.) are mapped to Enzo field names, listed in the :ref:`enzo-field-names`.
The reason Enzo field names are used here is because ``yt`` was originally written
to only read Enzo data.  In the future we will switch to a new system of
universal field names - this will also make it much easier to access raw on-disk
Enzo data!

In addition to primitive fields, yt provides an extensive list of "universal"
derived fields that are accessible from any of the frontends.  For a full
listing of the universal derived fields, see :ref:`universal-field-list`.

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

Which Enzo Field names Does ``yt`` Know About?
----------------------------------------------

These are the names of primitive fields in the Enzo AMR code.  ``yt`` was originally
written to analyze Enzo data so the default field names used by the various
frontends are the same as Enzo fields.

.. note::

   Enzo field names are *universal* yt fields.  All frontends define conversions
   to Enzo fields.  Enzo fields are always in CGS.

* Density
* Temperature
* Gas Energy
* Total Energy
* [xyz]-velocity
* Species fields: HI, HII, Electron, HeI, HeII, HeIII, HM, H2I, H2II, DI, DII, HDI
* Particle mass, velocity, 
