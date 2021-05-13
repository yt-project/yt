.. _creating_derived_quantities:

Creating Derived Quantities
---------------------------

.. warning:: This section is not yet updated to work with yt 3.0.  If you
             have a question about making a custom derived quantity, please
             contact the mailing list.

The basic idea is that you need to be able to operate both on a set of data,
and a set of sets of data.  (If this is not possible, the quantity needs to be
added with the ``force_unlazy`` option.)

Two functions are necessary.  One will operate on arrays of data, either fed
from each grid individually or fed from the entire data object at once.  The
second one takes the results of the first, either as lists of arrays or as
single arrays, and returns the final values.  For an example, we look at the
``TotalMass`` function:

.. code-block:: python

   def _TotalMass(data):
       baryon_mass = data["mass"].sum()
       particle_mass = data["particle_mass"].sum()
       return baryon_mass, particle_mass


   def _combTotalMass(data, baryon_mass, particle_mass):
       return baryon_mass.sum() + particle_mass.sum()


   add_quantity("TotalMass", function=_TotalMass, combine_function=_combTotalMass, n_ret=2)

Once the two functions have been defined, we then call :func:`add_quantity` to
tell it the function that defines the data, the collator function, and the
number of values that get passed between them.  In this case we return both the
particle and the baryon mass, so we have two total values passed from the main
function into the collator.
