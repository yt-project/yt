

.. function:: Action(action, combine_action, filter=None):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._Action`.)
   This function evals the string given by the action arg and uses 
   the function thrown with the combine_action to combine the values.  
   A filter can be thrown to be evaled to short-circuit the calculation 
   if some criterion is not met.
   :param action: a string containing the desired action to be evaled.
   :param combine_action: the function used to combine the answers when done lazily.
   :param filter: a string to be evaled to serve as a data filter.



.. function:: AngularMomentumVector():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._AngularMomentumVector`.)
   This function returns the mass-weighted average angular momentum vector.



.. function:: BaryonSpinParameter():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._BaryonSpinParameter`.)
   This function returns the spin parameter for the baryons, but it uses
   the particles in calculating enclosed mass.



.. function:: BulkVelocity():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._BulkVelocity`.)
   This function returns the mass-weighted average velocity in the object.



.. function:: CenterOfMass(use_cells=True, use_particles=False):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._CenterOfMass`.)
   This function returns the location of the center
   of mass. By default, it computes of the *non-particle* data in the object. 
   
   Parameters
   ----------
   
   use_cells : bool
       If True, will include the cell mass (default: True)
   use_particles : bool
       if True, will include the particles in the object (default: False)



.. function:: Extrema(fields, non_zero=False, filter=None):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._Extrema`.)
   This function returns the extrema of a set of fields
   
   :param fields: A field name, or a list of field names
   :param filter: a string to be evaled to serve as a data filter.



.. function:: IsBound(truncate=True, include_thermal_energy=False, treecode=True, opening_angle=1.0, periodic_test=False, include_particles=True):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._IsBound`.)
   This returns whether or not the object is gravitationally bound. If this
   returns a value greater than one, it is bound, and otherwise not.
   
   Parameters
   ----------
   truncate : Bool
       Should the calculation stop once the ratio of
       gravitational:kinetic is 1.0?
   include_thermal_energy : Bool
       Should we add the energy from ThermalEnergy
       on to the kinetic energy to calculate 
       binding energy?
   treecode : Bool
       Whether or not to use the treecode.
   opening_angle : Float 
       The maximal angle a remote node may subtend in order
       for the treecode method of mass conglomeration may be
       used to calculate the potential between masses.
   periodic_test : Bool 
       Used for testing the periodic adjustment machinery
       of this derived quantity.
   include_particles : Bool
       Should we add the mass contribution of particles
       to calculate binding energy?
   
   Examples
   --------
   >>> sp.quantities["IsBound"](truncate=False,
   ... include_thermal_energy=True, treecode=False, opening_angle=2.0)
   0.32493



.. function:: MaxLocation(field):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._MaxLocation`.)
   This function returns the location of the maximum of a set
   of fields.



.. function:: MinLocation(field):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._MinLocation`.)
   This function returns the location of the minimum of a set
   of fields.



.. function:: ParticleSpinParameter():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._ParticleSpinParameter`.)
   This function returns the spin parameter for the baryons, but it uses
   the particles in calculating enclosed mass.



.. function:: StarAngularMomentumVector():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._StarAngularMomentumVector`.)
   This function returns the mass-weighted average angular momentum vector 
   for stars.



.. function:: TotalMass():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._TotalMass`.)
   This function takes no arguments and returns the sum of cell masses and
   particle masses in the object.



.. function:: TotalQuantity(fields):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._TotalQuantity`.)
   This function sums up a given field over the entire region
   
   :param fields: The fields to sum up



.. function:: WeightedAverageQuantity(field, weight):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._WeightedAverageQuantity`.)
   This function returns an averaged quantity.
   
   :param field: The field to average
   :param weight: The field to weight by

.. function:: WeightedVariance(field, weight):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities._WeightedVariance`.)
    This function returns the variance of a field.

    :param field: The target field
    :param weight: The field to weight by

    Returns the weighted variance and the weighted mean.
