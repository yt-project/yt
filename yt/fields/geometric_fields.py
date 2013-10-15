"""
Geometric fields live here.  Not coordinate ones, but ones that perform
transformations.



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

from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component, \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, get_sph_r, \
    get_sph_theta, get_sph_phi, \
    periodic_dist, euclidean_dist

@register_field_plugin
def setup_geometric_fields(registry, ftype = "gas", slice_info = None):
    ### spherical coordinates: r (radius)
    def _spherical_r(field, data):
        center = data.get_field_parameter("center")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_sph_r(coords)

    registry.add_field("spherical_r",
             function=_spherical_r,
             validators=[ValidateParameter("center")],
             units="cm")

    ### spherical coordinates: theta (angle with respect to normal)
    def _spherical_theta(field, data):
        center = data.get_field_parameter("center")
        normal = data.get_field_parameter("normal")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_sph_theta(coords, normal)

    registry.add_field("spherical_theta",
             function=_spherical_theta,
             validators=[ValidateParameter("center"),
             ValidateParameter("normal")])

    ### spherical coordinates: phi (angle in the plane perpendicular to the normal)
    def _spherical_phi(field, data):
        center = data.get_field_parameter("center")
        normal = data.get_field_parameter("normal")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_sph_phi(coords, normal)

    registry.add_field("spherical_phi",
             function=_spherical_phi,
             validators=[ValidateParameter("center"),
             ValidateParameter("normal")])

    ### cylindrical coordinates: R (radius in the cylinder's plane)
    def _cylindrical_r(field, data):
        center = data.get_field_parameter("center")
        normal = data.get_field_parameter("normal")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_cyl_r(coords, normal)

    registry.add_field("cylindrical_r",
             function=_cylindrical_r,
             validators=[ValidateParameter("center"),
                        ValidateParameter("normal")],
             units="cm")

    ### cylindrical coordinates: z (height above the cylinder's plane)
    def _cylindrical_z(field, data):
        center = data.get_field_parameter("center")
        normal = data.get_field_parameter("normal")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_cyl_z(coords, normal)

    registry.add_field("cylindrical_z",
             function=_cylindrical_z,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units="cm")

    ### cylindrical coordinates: theta (angle in the cylinder's plane)
    def _cylindrical_theta(field, data):
        center = data.get_field_parameter("center")
        normal = data.get_field_parameter("normal")
        coords = obtain_rvec(data)
        coords[0,...] -= center[0]
        coords[1,...] -= center[1]
        coords[2,...] -= center[2]
        return get_cyl_theta(coords, normal)

    registry.add_field("cylindrical_theta",
             function=_cylindrical_theta,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units = "")

    ### The old field DiskAngle is the same as the spherical coordinates'
    ### 'theta' angle. I'm keeping DiskAngle for backwards compatibility.
    # @todo: remove in 3.0?
    def _disk_angle(field, data):
        return data["spherical_theta"]

    registry.add_field("disk_angle",
              function=_disk_angle, take_log=False,
              validators=[ValidateParameter("center"),
                          ValidateParameter("normal")],
              display_field=False,
              units = "")

    ### The old field Height is the same as the cylindrical coordinates' z
    ### field. I'm keeping Height for backwards compatibility.
    # @todo: remove in 3.0?
    def _height(field, data):
        return data["cylindrical_z"]

    registry.add_field("height", function=_height, 
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units="cm",
             display_field=False)

