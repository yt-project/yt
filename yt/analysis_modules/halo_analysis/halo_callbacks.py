"""
Halo callback object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py

from yt.data_objects.profiles import \
     Profile1D
from yt.data_objects.yt_array import \
     YTArray, YTQuantity
from yt.funcs import \
     ensure_list
from yt.utilities.exceptions import YTUnitConversionError
from yt.utilities.logger import ytLogger as mylog
     
from .operator_registry import \
    callback_registry

def add_callback(name, function):
    callback_registry[name] =  HaloCallback(function)

class HaloCallback(object):
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, halo):
        self.function(halo, *self.args, **self.kwargs)
        return True

def halo_sphere(halo, radius_field="virial_radius", factor=1.0):
    r"""
    Create a sphere data container to associate with a halo.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    radius_field : string
        Field to be retrieved from the quantities dictionary as 
        the basis of the halo radius.
        Default: "virial_radius".
    factor : float
        Factor to be multiplied by the base radius for defining 
        the radius of the sphere.
        Defautl: 1.0.
        
    """

    dpf = halo.halo_catalog.data_pf
    hpf = halo.halo_catalog.halos_pf
    center = dpf.arr([halo.quantities['particle_position_%s' % axis] \
                      for axis in "xyz"]) / dpf.length_unit
    radius = factor * halo.quantities[radius_field] / dpf.length_unit
    sphere = dpf.h.sphere(center, (radius, "code_length"))
    setattr(halo, "data_object", sphere)

add_callback("sphere", halo_sphere)

def profile(halo, x_field, y_fields, x_bins=32, x_range=None, x_log=True,
            weight_field="cell_mass"):
    r"""
    Create a 1d profile.

    Parameter
    ---------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    x_field : string
        The binning field for the profile.
    y_fields : string or list of strings
        The fields to be profiled.
    x_bins : int
        The number of bins in the profile.
        Default: 32
    x_range : (float, float)
        The range of the x_field.  If None, the extrema are used.
        Default: None
    x_log : bool
        Flag for logarithmmic binning.
        Default: True
    weight_field : string
        Weight field for profiling.
        Default : "cell_mass"

    """

    dpf = halo.halo_catalog.data_pf
    
    if x_range is None:
        x_range = halo.data_object.quantities["Extrema"](x_field)[0]
        # temporary check until derived quantities are fixed
        # right now derived quantities do not work with units, so once they do, let us know
        try:
            x_range[0]._unit_repr_check_same("cm")
            raise RuntimeError("Looks like derived quantities have been fixed.  Fix this code!")
        except YTUnitConversionError:
            # for now, Extrema return dimensionless, but assume it is code_length
            x_range = [dpf.arr(x.to_ndarray(), "code_length") for x in x_range]
            print halo.data_object.radius.in_units("Mpc")
            print x_range
        #    my_profile = Profile1D(halo.data_object
    return

add_callback("profile", profile)
