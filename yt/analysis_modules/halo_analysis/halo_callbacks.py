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
     create_profile
from yt.data_objects.yt_array import \
     YTArray, YTQuantity

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
    radius_field : str
        Field to be retrieved from the quantities dictionary as 
        the basis of the halo radius.
        Default: "virial_radius".
    factor : float
        Factor to be multiplied by the base radius for defining 
        the radius of the sphere.
        Defautl: 1.0.
        
    """

    center = YTArray([halo.quantities['particle_position_%s' % axis] \
                      for axis in "xyz"]) / \
      halo.halo_catalog.data_pf.length_unit
    radius = factor * halo.quantities[radius_field] / \
      halo.halo_catalog.data_pf.length_unit
    sphere = halo.halo_catalog.data_pf.h.sphere(center, (radius, "code_length"))
    setattr(halo, "object", sphere)

add_callback("sphere", halo_sphere)
