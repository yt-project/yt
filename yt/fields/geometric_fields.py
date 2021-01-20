from typing import Any, Optional, Union

import numpy as np
from numpy import ndarray
from unyt.array import unyt_array

from yt.fields.derived_field import DerivedField
from yt.fields.field_detector import FieldDetector
from yt.frontends.ramses.data_structures import RAMSESDomainSubset
from yt.frontends.ramses.fields import RAMSESFieldInfo
from yt.utilities.lib.geometry_utils import compute_morton
from yt.utilities.math_utils import (
    get_cyl_r,
    get_cyl_theta,
    get_cyl_z,
    get_sph_phi,
    get_sph_r,
    get_sph_theta,
)

from .derived_field import ValidateParameter, ValidateSpatial
from .field_functions import get_periodic_rvec, get_radius
from .field_plugin_registry import register_field_plugin


@register_field_plugin
def setup_geometric_fields(
    registry: RAMSESFieldInfo, ftype: str = "gas", slice_info: Optional[Any] = None
) -> None:
    unit_system = registry.ds.unit_system

    def _radius(field: DerivedField, data: FieldDetector) -> unyt_array:
        """The spherical radius component of the mesh cells.

        Relative to the coordinate system defined by the *center* field
        parameter.
        """
        return get_radius(data, "", field.name[0])

    registry.add_field(
        ("index", "radius"),
        sampling_type="cell",
        function=_radius,
        validators=[ValidateParameter("center")],
        units=unit_system["length"],
    )

    def _grid_level(
        field: DerivedField, data: Union[FieldDetector, RAMSESDomainSubset]
    ) -> ndarray:
        """The AMR refinement level"""
        arr = np.ones(data.ires.shape, dtype="float64")
        arr *= data.ires
        if data._spatial:
            return data._reshape_vals(arr)
        return arr

    registry.add_field(
        ("index", "grid_level"),
        sampling_type="cell",
        function=_grid_level,
        units="",
        validators=[ValidateSpatial(0)],
    )

    def _grid_indices(field: DerivedField, data: FieldDetector) -> ndarray:
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

    registry.add_field(
        ("index", "grid_indices"),
        sampling_type="cell",
        function=_grid_indices,
        units="",
        validators=[ValidateSpatial(0)],
        take_log=False,
    )

    def _ones_over_dx(field: DerivedField, data: FieldDetector) -> unyt_array:
        """The inverse of the local cell spacing"""
        return (
            np.ones(data["index", "ones"].shape, dtype="float64") / data["index", "dx"]
        )

    registry.add_field(
        ("index", "ones_over_dx"),
        sampling_type="cell",
        function=_ones_over_dx,
        units=unit_system["length"] ** -1,
        display_field=False,
    )

    def _zeros(field: DerivedField, data: FieldDetector) -> unyt_array:
        """Returns zero for all cells"""
        arr = np.zeros(data["index", "ones"].shape, dtype="float64")
        return data.apply_units(arr, field.units)

    registry.add_field(
        ("index", "zeros"),
        sampling_type="cell",
        function=_zeros,
        units="",
        display_field=False,
    )

    def _ones(field: DerivedField, data: FieldDetector) -> unyt_array:
        """Returns one for all cells"""
        tmp = np.ones(data.ires.shape, dtype="float64")
        arr = data.apply_units(tmp, field.units)
        if data._spatial:
            return data._reshape_vals(arr)
        return arr

    registry.add_field(
        ("index", "ones"),
        sampling_type="cell",
        function=_ones,
        units="",
        display_field=False,
    )

    def _morton_index(field: DerivedField, data: FieldDetector) -> ndarray:
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
        morton = compute_morton(
            data["index", "x"].ravel(),
            data["index", "y"].ravel(),
            data["index", "z"].ravel(),
            LE,
            RE,
        )
        morton.shape = data["index", "x"].shape
        return morton.view("f8")

    registry.add_field(
        ("index", "morton_index"),
        sampling_type="cell",
        function=_morton_index,
        units="",
    )

    def _spherical_radius(field: DerivedField, data: FieldDetector) -> unyt_array:
        """The spherical radius component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *center* field
        parameter.
        """
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_sph_r(coords), "code_length").in_base(unit_system.name)

    registry.add_field(
        ("index", "spherical_radius"),
        sampling_type="cell",
        function=_spherical_radius,
        validators=[ValidateParameter("center")],
        units=unit_system["length"],
    )

    def _spherical_r(field: DerivedField, data: FieldDetector) -> unyt_array:
        """This field is deprecated and will be removed in a future release"""
        return data["index", "spherical_radius"]

    registry.add_field(
        ("index", "spherical_r"),
        sampling_type="cell",
        function=_spherical_r,
        validators=[ValidateParameter("center")],
        units=unit_system["length"],
    )

    def _spherical_theta(field: DerivedField, data: FieldDetector) -> ndarray:
        """The spherical theta component of the positions of the mesh cells.

        theta is the poloidal position angle in the plane parallel to the
        *normal* vector

        Relative to the coordinate system defined by the *center* and *normal*
        field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_theta(coords, normal)

    registry.add_field(
        ("index", "spherical_theta"),
        sampling_type="cell",
        function=_spherical_theta,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units="",
    )

    def _spherical_phi(field: DerivedField, data: FieldDetector) -> ndarray:
        """The spherical phi component of the positions of the mesh cells.

        phi is the azimuthal position angle in the plane perpendicular to the
        *normal* vector

        Relative to the coordinate system defined by the *center* and *normal*
        field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_sph_phi(coords, normal)

    registry.add_field(
        ("index", "spherical_phi"),
        sampling_type="cell",
        function=_spherical_phi,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units="",
    )

    def _cylindrical_radius(field: DerivedField, data: FieldDetector) -> unyt_array:
        """The cylindrical radius component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_r(coords, normal), "code_length").in_base(
            unit_system.name
        )

    registry.add_field(
        ("index", "cylindrical_radius"),
        sampling_type="cell",
        function=_cylindrical_radius,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units=unit_system["length"],
    )

    def _cylindrical_r(field: DerivedField, data: FieldDetector) -> unyt_array:
        """This field is deprecated and will be removed in a future release"""
        return data["index", "cylindrical_radius"]

    registry.add_field(
        ("index", "cylindrical_r"),
        sampling_type="cell",
        function=_cylindrical_r,
        validators=[ValidateParameter("center")],
        units=unit_system["length"],
    )

    def _cylindrical_z(field: DerivedField, data: FieldDetector) -> unyt_array:
        """The cylindrical z component of the positions of the mesh cells.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return data.ds.arr(get_cyl_z(coords, normal), "code_length").in_base(
            unit_system.name
        )

    registry.add_field(
        ("index", "cylindrical_z"),
        sampling_type="cell",
        function=_cylindrical_z,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units=unit_system["length"],
    )

    def _cylindrical_theta(field: DerivedField, data: FieldDetector) -> ndarray:
        """The cylindrical z component of the positions of the mesh cells.

        theta is the azimuthal position angle in the plane perpendicular to the
        *normal* vector.

        Relative to the coordinate system defined by the *normal* vector and
        *center* field parameters.
        """
        normal = data.get_field_parameter("normal")
        coords = get_periodic_rvec(data)
        return get_cyl_theta(coords, normal)

    registry.add_field(
        ("index", "cylindrical_theta"),
        sampling_type="cell",
        function=_cylindrical_theta,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units="",
    )

    def _disk_angle(field: DerivedField, data: FieldDetector) -> ndarray:
        """This field is dprecated and will be removed in a future release"""
        return data["index", "spherical_theta"]

    registry.add_field(
        ("index", "disk_angle"),
        sampling_type="cell",
        function=_disk_angle,
        take_log=False,
        display_field=False,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units="",
    )

    def _height(field: DerivedField, data: FieldDetector) -> unyt_array:
        """This field is deprecated and will be removed in a future release"""
        return data["index", "cylindrical_z"]

    registry.add_field(
        ("index", "height"),
        sampling_type="cell",
        function=_height,
        validators=[ValidateParameter("center"), ValidateParameter("normal")],
        units=unit_system["length"],
        display_field=False,
    )
