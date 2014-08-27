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
    ValidateParameter, \
    ValidateGridType, \
    ValidateSpatial

from .field_functions import \
     get_periodic_rvec, \
     get_radius

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

    def _radius(field, data):
        return get_radius(data, "")

    registry.add_field(("index", "radius"), function=_radius,
              validators=[ValidateParameter("center")],
              units="cm")

    def _grid_level(field, data):
        return np.ones(data.ActiveDimensions)*(data.Level)
    registry.add_field(("index", "grid_level"),
              function=_grid_level, units = "",
              validators=[ValidateGridType(),
                          ValidateSpatial(0)])

    def _grid_indices(field, data):
        return np.ones(data["index", "ones"].shape)*(data.id-data._id_offset)
    registry.add_field(("index", "grid_indices"),
              function=_grid_indices, units = "",
              validators=[ValidateGridType(),
                          ValidateSpatial(0)], take_log=False)
    def _ones_over_dx(field, data):
        return np.ones(data["index", "ones"].shape,
                       dtype="float64")/data["index", "dx"]
    registry.add_field(("index", "ones_over_dx"), function=_ones_over_dx,
              units = "1 / cm",
              display_field=False)

    def _zeros(field, data):
        arr = np.zeros(data["index", "ones"].shape, dtype='float64')
        return data.apply_units(arr, field.units)

    registry.add_field(("index", "zeros"), function=_zeros,
              units = "",
              display_field=False)

    def _ones(field, data):
        arr = np.ones(data.ires.shape, dtype="float64")
        if data._spatial:
            return data._reshape_vals(arr)
        return data.apply_units(arr, field.units)

    registry.add_field(("index", "ones"), function=_ones,
              units = "",
              display_field=False)

    ### spherical coordinates: r (radius)
    def _spherical_r(field, data):
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_sph_r(coords), "code_length").in_cgs()

    registry.add_field(("index", "spherical_r"),
             function=_spherical_r,
             validators=[ValidateParameter("center")],
             units="cm")

    ### spherical coordinates: theta (angle with respect to normal)
    def _spherical_theta(field, data):
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_theta(coords, normal)

    registry.add_field(("index", "spherical_theta"),
             function=_spherical_theta,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")])

    ### spherical coordinates: phi (angle in the plane perpendicular to the normal)
    def _spherical_phi(field, data):
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_phi(coords, normal)

    registry.add_field(("index", "spherical_phi"),
             function=_spherical_phi,
             validators=[ValidateParameter("center"),
             ValidateParameter("normal")])

    ### cylindrical coordinates: R (radius in the cylinder's plane)
    def _cylindrical_r(field, data):
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_r(coords, normal), "code_length").in_cgs()

    registry.add_field(("index", "cylindrical_r"),
             function=_cylindrical_r,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units="cm")

    ### cylindrical coordinates: z (height above the cylinder's plane)
    def _cylindrical_z(field, data):
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_z(coords, normal), "code_length").in_cgs()

    registry.add_field(("index", "cylindrical_z"),
             function=_cylindrical_z,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units="cm")

    ### cylindrical coordinates: theta (angle in the cylinder's plane)
    def _cylindrical_theta(field, data):
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_cyl_theta(coords, normal)

    registry.add_field(("index", "cylindrical_theta"),
             function=_cylindrical_theta,
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units = "")

    ### The old field DiskAngle is the same as the spherical coordinates'
    ### 'theta' angle. I'm keeping DiskAngle for backwards compatibility.
    # @todo: remove in 3.0?
    def _disk_angle(field, data):
        return data["index", "spherical_theta"]

    registry.add_field(("index", "disk_angle"),
              function=_disk_angle, take_log=False,
              validators=[ValidateParameter("center"),
                          ValidateParameter("normal")],
              display_field=False,
              units = "")

    ### The old field Height is the same as the cylindrical coordinates' z
    ### field. I'm keeping Height for backwards compatibility.
    # @todo: remove in 3.0?
    def _height(field, data):
        return data["index", "cylindrical_z"]

    registry.add_field(("index", "height"), function=_height, 
             validators=[ValidateParameter("center"),
                         ValidateParameter("normal")],
             units="cm",
             display_field=False)

