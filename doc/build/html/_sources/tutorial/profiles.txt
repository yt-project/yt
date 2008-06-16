Profiles
========

Profiles are one of the more interesting data types.  They are a binned average
of one or more variables, along one, two or three different dimensions.  For example,
you could get an average temperature in Density and H2 Fraction space.  Or you could
Look at the distribution of Volume through Density and Temperature space.  Or you could
take the spherically-averaged Number Density and see how it changes with radius.

One Dimensional Profiles
------------------------

In order to make a profile, we first need a data object we want to get the profile
*from*.  This can be a sphere, a region, or even an extracted region.  We'll
first take a look at doing a spherically-average profile of Temperature with Radius.  ::

   >>> sphere = a.h.sphere([0.5,0.5,0.5], 0.1/a["kpc"])
   >>> prof = lagos.BinnedProfile1D(sphere, 16, "Radiuskpc",
                                 0.0001, 0.1, True, True)
   >>> prof.add_fields(["Density","Temperature"], weight="CellMassMsun")


Alright, again, a bit of explanation.  First, we get our sphere -- and in doing
so, we tell it where the center is, and how big (in *code units*) we want it to be.
(Note that in order to convert a number into some unit, we *divide* by the conversion
factor for that unit.

Once we have that data object, we can do anything with it.  (See above!)
However, what we *want* to do is generate a profile.  So we create
a profile object, telling it where to get the data, which field we want to bin by,
the lower and upper bounds (in the units of that field), and whether we want those
bins to be log spaced.  But, there's one more option -- the last one is the "lazy_reader"
option.  This controls whether we load all the data at once (faster) or whether
we grab it from each grid individually (slower, but more memory-conservative.) 

Once we have instantiated the profile object, we add fields to it.
These fields are weighted by CellMassMsun, and we get back an average value for
each bin in 'Radiuskpc'.  Now we can plot them, or inspect them, with simple
data-access methods.  Additionally, you can feed them directly into pylab or matplotlib
(neither of which we will cover in any length here) and plot them: ::

   >>> print prof["Density"]
   >>> import pylab
   >>> pylab.loglog(prof["Radiuskpc"], prof["Density"])
   >>> pylab.xlabel(r"$\rm{Radius} (\rm{kpc})$")
   >>> pylab.ylabel(r"$\rm{Density} (\rm{g}/\rm{cm}^3)$")
   >>> pylab.savefig("my_profile.png")

Note that we manually plot it, but :mod:`raven` also includes an automated
plotting method for obtaining and displaying radial profiles.

Two-Dimensional Profiles
------------------------

Two-dimensional profiles are along the same lines.  ::

   >>> prof2d = lagos.BinnedProfile2D(sphere,
   ...                 16, "Density", 1e-32, 1e-24, True,
   ...                 16, "Temperature", 1e2, 1e5, True, True)
   >>> prof2d.add_fields("CellVolume", weight=None)
   >>> pylab.pcolormesh(prof2d["Density"], prof2d["Temperature"],
   ...                  prof2d["CellVolume"].transpose())

Again, we have manually plotted this, but an automated plotter is available.::

   >>> pc = raven.PlotCollection(a)
   >>> pc.add_phase_object(sphere, ["Density", "Temperature", "x-velocity"],
   ...                 16, True, (1e-32, 1e-24), 
   ...                 16, True, (1e2, 1e5), lazy_reader=True)

The API is slightly different here, but it will accept any object (a sphere in
our case) and generate the profile based on the arguments.  If you simply want
a sphere, you can ask for that and the sphere will be generated as well -- and
the function will automatically generate the bounds.  (Although, keep in mind,
if you do not supply bounds it will have to read the data into memory, which
means you are no longer able to use the :keyword:`lazy_reader` functionality.
This will change in future versions.) ::


   >>> pc.add_phase_sphere(0.1, 'kpc', ["Density", "Temperature", "x-velocity"])

Three-Dimensional Profiles
--------------------------

.. image:: ../_images/example_3dphase.png
   :align: left
   :width: 250

If you have the `S2PLOT <http://astronomy.swin.edu.au/s2plot/index.php?title=S2PLOT>`_
bindings installed, you can visualize three-dimensional profiles with yt.  The
API is nearly identical to that used to generate the lower-dimensionality
profiles. ::

   >>> prof = lagos.BinnedProfile3D(sphere,
   ...                 128, "Density",     extrema[0][0], extrema[0][1], True,
   ...                 128, "Temperature", extrema[1][0], extrema[1][1], True,
   ...                 128, "z-velocity",  extrema[2][0], extrema[2][1], True,
   ...                 lazy_reader=True)
   >>> prof.add_fields("CellMassMsun", None)

I've wrapped the bounds into the variable *extrema*, to make it more clear.
Once you have this profile, you can plot it with S2PLOT by calling initializing
a 3D plot and calling :func:`run`. ::

   >>> vrc = raven.VolumeRendering3DProfile(k, "CellMassMsun", amax=1.0)
   >>> vrc.run()

They can be generated anywhere -- and serialized -- so if you want to generate
them on a remote machine and use a local-machine for display, you only have to
call :func:`store` and then initialize, instead of the standard
:class:`BinnedProfile3D`, :class:`StoredBinnedProfile3D`.  ::

   >>> prof.store_profile('MyProfile')
   >>> new_prof = lagos.StoredBinnedProfile3D(a, 'MyProfile')

The retrieved profile is static, and new data cannot be added.  However, it
becomes useful in the likely circumstance that your data is stored on a machine
on which is inconvenient to do OpenGL-based data exploration.
