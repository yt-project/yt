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
    ValidateSpatial

from .field_functions import \
    get_periodic_rvec, \
    get_radius

from .field_plugin_registry import \
    register_field_plugin

from yt.utilities.math_utils import \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, get_sph_r, \
    get_sph_theta, get_sph_phi

from yt.utilities.lib.geometry_utils import \
    compute_morton

@register_field_plugin
def setup_geometric_fields(registry, ftype="gas", slice_info=None):
    unit_system = registry.ds.unit_system
    def _radius(field, data):
        """The spherical radius component of the mesh cells.

        Relative to the coordinate system defined by the *center* field
        parameter.
        """
        return get_radius(data, "")

    registry.add_field(("index", "radius"),
                       function=_radius,
                       validators=[ValidateParameter("center")],
                       units=unit_system["length"])

    def _grid_level(field, data):
        """The AMR refinement level"""
        arr = np.ones(data.ires.shape, dtype="float64") 
        arr *= data.ires
        if data._spatial:
            return data._reshape_vals(arr)
        return arr

    registry.add_field(("index", "grid_level"),
                       function=_grid_level,
                       units="",
                       validators=[ValidateSpatial(0)])

    def _grid_indices(field, data):
        """The index of the leaf grid the mesh cells exist on"""
        if hasattr(data, "domain_id"):
            this_id = data.domain_id
        else:
            this_id = data.id - data._id_offset
        arr = np.ones(data["index", "ones"].shape)
        arr *= this_id
        if data._spatial:
            return data._reshape_vals(arr)
        return arr

    registry.add_field(("index", "grid_indices"),
                       function=_grid_indices,
                       units="",
                       validators=[ValidateSpatial(0)],
                       take_log=False)

    def _ones_over_dx(field, data):
        """The inverse of the local cell spacing"""
        return np.ones(data["index", "ones"].shape,
                       dtype="float64")/data["index", "dx"]

    registry.add_field(("index", "ones_over_dx"),
                       function=_ones_over_dx,
                       units=unit_system["length"]**-1,
                       display_field=False)

    def _zeros(field, data):
        """Returns zero for all cells"""
        arr = np.zeros(data["index", "ones"].shape, dtype='float64')
        return data.apply_units(arr, field.units)

    registry.add_field(("index", "zeros"),
                       function=_zeros,
                       units="",
                       display_field=False)

    def _ones(field, data):
        """Returns one for all cells"""
        arr = np.ones(data.ires.shape, dtype="float64")
        if data._spatial:
            return data._reshape_vals(arr)
        return data.apply_units(arr, field.units)

    registry.add_field(("index", "ones"),
                       function=_ones,
                       units="",
                       display_field=False)

    def _morton_index(field, data):
        """This is the morton index, which is properly a uint64 field.  Because
        we make some assumptions that the fields returned by derived fields are
        float64, this returns a "view" on the data that is float64.  To get
        back the original uint64, you need to call .view("uint64") on it;
        however, it should be true that if you sort the uint64, you will get
        the same order as if you sort the float64 view.
        """
        eps = np.finfo("f8").eps
        uq = data.ds.domain_left_edge.uq
        LE = data.ds.domain_left_edge - eps * uq
        RE = data.ds.domain_right_edge + eps * uq
        # .ravel() only copies if it needs to
        morton = compute_morton(data["index", "x"].ravel(),
                                data["index", "y"].ravel(),
                                data["index", "z"].ravel(), LE, RE)
        morton.shape = data["index", "x"].shape
        return morton.view("f8")
    registry.add_field(("index", "morton_index"), function=_morton_index,
                       units = "")
        
    def _spherical_radius(field, data):
        """The spherical radius component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *center* field
        parameter.
        """
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_sph_r(coords), "code_length").in_base(unit_system.name)

    registry.add_field(("index", "spherical_radius"),
                       function=_spherical_radius,
                       validators=[ValidateParameter("center")],
                       units=unit_system["length"])

    def _spherical_r(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data['index', 'spherical_radius']

    registry.add_field(("index", "spherical_r"),
                       function=_spherical_r,
                       validators=[ValidateParameter("center")],
                       units=unit_system["length"])

    def _spherical_theta(field, data):
        """The spherical theta component of the positions of the mesh cells.

        theta is the poloidal position angle in the plane parallel to the
        *normal* vector

        Relative to the coordinate system defined by the *center* and *normal*
        field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_theta(coords, normal)

    registry.add_field(("index", "spherical_theta"),
                       function=_spherical_theta,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units='')

    def _spherical_phi(field, data):
        """The spherical phi component of the positions of the mesh cells.

        phi is the azimuthal position angle in the plane perpendicular to the
        *normal* vector

        Relative to the coordinate system defined by the *center* and *normal*
        field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_phi(coords, normal)

    registry.add_field(("index", "spherical_phi"),
                       function=_spherical_phi,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units='')

    def _cylindrical_radius(field, data):
        """The cylindrical radius component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_r(coords, normal), "code_length").in_base(unit_system.name)

    registry.add_field(("index", "cylindrical_radius"),
                       function=_cylindrical_radius,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units=unit_system["length"])

    def _cylindrical_r(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data['index', 'cylindrical_radius']

    registry.add_field(("index", "cylindrical_r"),
                       function=_cylindrical_r,
                       validators=[ValidateParameter("center")],
                       units=unit_system["length"])

    def _cylindrical_z(field, data):
        """The cylindrical z component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_z(coords, normal), "code_length").in_base(unit_system.name)

    registry.add_field(("index", "cylindrical_z"),
                       function=_cylindrical_z,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units=unit_system["length"])

    def _cylindrical_theta(field, data):
        """The cylindrical z component of the positions of the mesh cells.

        theta is the azimuthal position angle in the plane perpendicular to the
        *normal* vector.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_cyl_theta(coords, normal)

    registry.add_field(("index", "cylindrical_theta"),
                       function=_cylindrical_theta,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units="")

    def _disk_angle(field, data):
        """This field is dprecated and will be removed in a future release"""
        return data["index", "spherical_theta"]

    registry.add_field(("index", "disk_angle"),
                       function=_disk_angle,
                       take_log=False,
                       display_field=False,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units="")

    def _height(field, data):
        """This field is deprecated and will be removed in a future release"""
        return data["index", "cylindrical_z"]

    registry.add_field(("index", "height"),
                       function=_height,
                       validators=[ValidateParameter("center"),
                                   ValidateParameter("normal")],
                       units=unit_system["length"],
                       display_field=False)
