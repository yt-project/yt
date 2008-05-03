Profiles
--------

Profiles are one of the more interesting data types.  They are a binned average
of one or more variables, along one or two different dimensions.  For example,
you could get an average temperature in Density and H2 Fraction space.  Or you could
Look at the distribution of Volume through Density and Temperature space.  Or you could
take the spherically-averaged Number Density and see how it changes with radius.

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


Two-dimensional profiles are along the same lines.  ::

   >>> prof2d = lagos.BinnedProfile2D(sphere,
   ...                 16, "Density", 1e-32, 1e-24, True,
   ...                 16, "Temperature", 1e2, 1e5, True, True)
   >>> prof2d.add_fields("CellVolume", weight=None)
   >>> pylab.pcolormesh(prof2d["Density"], prof2d["Temperature"],
   ...                  prof2d["CellVolume"].transpose())
