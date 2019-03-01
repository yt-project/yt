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

import numpy as np

from .derived_field import \
    ValidateParameter

from .field_plugin_registry import \
    register_field_plugin

from .vector_operations import \
    create_magnitude_field

from yt.utilities.lib.misc_utilities import \
    obtain_position_vector, \
    obtain_relative_velocity_vector

@register_field_plugin
def setup_angular_momentum(registry, ftype = "gas", slice_info = None):
    unit_system = registry.ds.unit_system
    def _specific_angular_momentum_x(field, data):
        xv, yv, zv = obtain_relative_velocity_vector(data)
        rv = obtain_position_vector(data)
        units = rv.units
        rv = np.rollaxis(rv, 0, len(rv.shape))
        rv = data.ds.arr(rv, input_units=units)
        return yv * rv[...,2] - zv * rv[...,1]

    def _specific_angular_momentum_y(field, data):
        xv, yv, zv = obtain_relative_velocity_vector(data)
        rv = obtain_position_vector(data)
        units = rv.units
        rv = np.rollaxis(rv, 0, len(rv.shape))
        rv = data.ds.arr(rv, input_units=units)
        return - (xv * rv[...,2] - zv * rv[...,0])

    def _specific_angular_momentum_z(field, data):
        xv, yv, zv = obtain_relative_velocity_vector(data)
        rv = obtain_position_vector(data)
        units = rv.units
        rv = np.rollaxis(rv, 0, len(rv.shape))
        rv = data.ds.arr(rv, input_units=units)
        return xv * rv[...,1] - yv * rv[...,0]

    registry.add_field((ftype, "specific_angular_momentum_x"),
                       sampling_type="local",
                       function=_specific_angular_momentum_x,
                       units=unit_system["specific_angular_momentum"],
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("bulk_velocity")])

    registry.add_field((ftype, "specific_angular_momentum_y"),
                       sampling_type="local",
                       function=_specific_angular_momentum_y,
                       units=unit_system["specific_angular_momentum"],
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("bulk_velocity")])

    registry.add_field((ftype, "specific_angular_momentum_z"),
                       sampling_type="local",
                       function=_specific_angular_momentum_z,
                       units=unit_system["specific_angular_momentum"],
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("bulk_velocity")])

    create_magnitude_field(registry, "specific_angular_momentum",
                           unit_system["specific_angular_momentum"], ftype=ftype)

    def _angular_momentum_x(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_x"]

    registry.add_field((ftype, "angular_momentum_x"),
                       sampling_type="cell",
                       function=_angular_momentum_x,
                       units=unit_system["angular_momentum"],
                       validators=[ValidateParameter('center'),
                                   ValidateParameter("bulk_velocity")])

    def _angular_momentum_y(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_y"]

    registry.add_field((ftype, "angular_momentum_y"),
                       sampling_type="local",
                       function=_angular_momentum_y,
                       units=unit_system["angular_momentum"],
                       validators=[ValidateParameter('center'),
                                   ValidateParameter("bulk_velocity")])

    def _angular_momentum_z(field, data):
        return data[ftype, "cell_mass"] \
             * data[ftype, "specific_angular_momentum_z"]

    registry.add_field((ftype, "angular_momentum_z"),
                       sampling_type="local",
                       function=_angular_momentum_z,
                       units=unit_system["angular_momentum"],
                       validators=[ValidateParameter('center'),
                                   ValidateParameter("bulk_velocity")])

    create_magnitude_field(registry, "angular_momentum",
                           unit_system["angular_momentum"], ftype=ftype)
