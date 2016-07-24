"""
Halo recipe object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.operator_registry import \
    OperatorRegistry

recipe_registry = OperatorRegistry()

def add_recipe(name, function):
    recipe_registry[name] =  HaloRecipe(function)

class HaloRecipe(object):
    r"""
    A HaloRecipe is a function that minimally takes in a Halo object
    and performs some analysis on it.  This function may attach attributes
    to the Halo object, write out data, etc, but does not return anything.
    """
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, halo_catalog):
        return self.function(halo_catalog, *self.args, **self.kwargs)

def calculate_virial_quantities(hc, fields,
                                weight_field=None, accumulation=True,
                                radius_field="virial_radius", factor=2.0,
                                overdensity_field=("gas", "overdensity"),
                                critical_overdensity=200):
    r"""
    Calculate virial quantities with the following procedure:
    1. Create a sphere data container.
    2. Create 1D radial profiles of overdensity and any requested fields.
    3. Call virial_quantities callback to interpolate profiles for value of critical overdensity.
    4. Delete profile and sphere objects from halo.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    fields: string or list of strings
        The fields for which virial values are to be calculated.
    weight_field : string
        Weight field for profiling.
        Default : "cell_mass"
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of
        all the values from bin 0 to n.  If -True, the sum is reversed so
        that the value for bin n is the cumulative sum from bin N (total bins)
        to n.  If the profile is 2D or 3D, a list of values can be given to
        control the summation in each dimension independently.
        Default: False.
    radius_field : string
        Field to be retrieved from the quantities dictionary as
        the basis of the halo radius.
        Default: "virial_radius".
    factor : float
        Factor to be multiplied by the base radius for defining
        the radius of the sphere.
        Default: 2.0.
    overdensity_field : string or tuple of strings
        The field used as the overdensity from which interpolation is done to
        calculate virial quantities.
        Default: ("gas", "overdensity")
    critical_overdensity : float
        The value of the overdensity at which to evaulate the virial quantities.
        Overdensity is with respect to the critical density.
        Default: 200

    """

    storage = "virial_quantities_profiles"
    pfields = [field for field in fields if field != "radius"]

    hc.add_callback("sphere", factor=factor)
    if pfields:
        hc.add_callback("profile", ["radius"], pfields,
                        weight_field=weight_field,
                        accumulation=accumulation,
                        storage=storage)
    hc.add_callback("profile", ["radius"], [overdensity_field],
                    weight_field="cell_volume", accumulation=True,
                    storage=storage)
    hc.add_callback("virial_quantities", fields,
                    overdensity_field=overdensity_field,
                    critical_overdensity=critical_overdensity,
                    profile_storage=storage)
    hc.add_callback("delete_attribute", storage)
    hc.add_callback("delete_attribute", "data_object")

add_recipe("calculate_virial_quantities", calculate_virial_quantities)
