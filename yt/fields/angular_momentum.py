"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import types
import numpy as np
import inspect
import copy

from .derived_field import \
    ValidateParameter

from .field_plugin_registry import \
    register_field_plugin

def obtain_velocities(data):
    return obtain_rv_vec(data)

@register_field_plugin
def setup_angular_momentum(registry, ftype = "gas", slice_info = None):
    def _specific_angular_momentum_x(field, data):
        xv, yv, zv = obtain_velocities(data)
        center = data.get_field_parameter('center')
        v_vec = obtain_rvec(data)
        v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
        rv = v_vec - center
        return yv * rv[...,2] - zv * rv[...,1]

    def _specific_angular_momentum_y(field, data):
        xv, yv, zv = obtain_velocities(data)
        center = data.get_field_parameter('center')
        v_vec = obtain_rvec(data)
        v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
        rv = v_vec - center
        return - (xv * rv[...,2] - zv * rv[...,0])

    def _specific_angular_momentum_z(field, data):
        xv, yv, zv = obtain_velocities(data)
        center = data.get_field_parameter('center')
        v_vec = obtain_rvec(data)
        v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
        rv = v_vec - center
        return xv * rv[...,1] - yv * rv[...,0]

    registry.add_field((ftype, "specific_angular_momentum_x"),
                        function=_specific_angular_momentum_x,
                        units="cm**2/s",
                        validators=[ValidateParameter("center")])
    registry.add_field((ftype, "specific_angular_momentum_y"),
                        function=_specific_angular_momentum_y,
                        units="cm**2/s",
                        validators=[ValidateParameter("center")])
    registry.add_field((ftype, "specific_angular_momentum_z"),
                        function=_specific_angular_momentum_z,
                        units="cm**2/s",
                        validators=[ValidateParameter("center")])

    def _angular_momentum_x(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_x"]
    registry.add_field((ftype, "angular_momentum_x"),
                       function=_angular_momentum_x,
                       units="g * cm**2 / s",
                       validators=[ValidateParameter('center')])

    def _angular_momentum_y(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_y"]
    registry.add_field((ftype, "angular_momentum_y"),
                       function=_angular_momentum_y,
                       units="g * cm**2 / s",
                       validators=[ValidateParameter('center')])

    def _angular_momentum_z(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_z"]
    registry.add_field((ftype, "angular_momentum_z"),
                       function=_angular_momentum_z,
                       units="g * cm**2 / s",
                       validators=[ValidateParameter('center')])

