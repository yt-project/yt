

.. function:: angular_momentum_vector()

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.AngularMomentumVector`.)
   This function returns the mass-weighted average angular momentum vector.


.. function:: bulk_velocity():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.BulkVelocity`.)
   This function returns the mass-weighted average velocity in the object.


.. function:: center_of_mass(use_cells=True, use_particles=False):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.CenterOfMass`.)
   This function returns the location of the center
   of mass. By default, it computes of the *non-particle* data in the object. 
   
   Parameters
   ----------
   
   use_cells : bool
       If True, will include the cell mass (default: True)
   use_particles : bool
       if True, will include the particles in the object (default: False)



.. function:: extrema(fields, non_zero=False, filter=None):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.Extrema`.)
   This function returns the extrema of a set of fields
   
   :param fields: A field name, or a list of field names
   :param filter: a string to be evaled to serve as a data filter.


.. function:: max_location(field):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.max_location`.)
   This function returns the location of the maximum of a set
   of fields.


.. function:: min_location(field):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.MinLocation`.)
   This function returns the location of the minimum of a set
   of fields.



.. function:: spin_parameter(use_gas=True, use_particles=True):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.SpinParameter`.)
   This function returns the spin parameter for the baryons, but it uses
   the particles in calculating enclosed mass.


.. function:: total_mass():

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.TotalMass`.)
   This function takes no arguments and returns the sum of cell masses and
   particle masses in the object.


.. function:: total_quantity(fields):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.TotalQuantity`.)
   This function sums up a given field over the entire region
   
   :param fields: The fields to sum up



.. function:: weighted_average_quantity(field, weight):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.WeightedAverageQuantity`.)
   This function returns an averaged quantity.
   
   :param field: The field to average
   :param weight: The field to weight by

.. function:: weighted_variance(field, weight):

   (This is a proxy for :func:`~yt.data_objects.derived_quantities.WeightedVariance`.)
    This function returns the variance of a field.

    :param field: The target field
    :param weight: The field to weight by

    Returns the weighted variance and the weighted mean.
