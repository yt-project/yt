from functools import cached_property
from typing import Dict

import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.utilities.lib.pixelization_routines import pixelize_cylinder

from ._axes_transforms import AxesTransform
from .coordinate_handler import (
    CoordinateHandler,
    DefaultProperties,
    _get_coord_fields,
    _get_polar_bounds,
    _setup_dummy_cartesian_coords_and_widths,
    _setup_polar_coordinates,
    cartesian_to_cylindrical,
    cylindrical_to_cartesian,
)

#
# Cylindrical fields
#


class CylindricalCoordinateHandler(CoordinateHandler):
    name = "cylindrical"
    _default_axis_order = ("r", "z", "theta")
    _default_axes_transforms: Dict[str, AxesTransform] = {
        "r": AxesTransform.GEOMETRY_NATIVE,
        "theta": AxesTransform.GEOMETRY_NATIVE,
        "z": AxesTransform.POLAR,
    }

    def __init__(self, ds, ordering=None):
        super().__init__(ds, ordering)

    def setup_fields(self, registry):
        # Missing implementation for x and y coordinates.
        _setup_dummy_cartesian_coords_and_widths(registry, axes=("x", "y"))
        _setup_polar_coordinates(registry, self.axis_id)

        f1, f2 = _get_coord_fields(self.axis_id["z"])
        registry.add_field(
            ("index", "dz"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="code_length",
        )

        registry.add_field(
            ("index", "z"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="code_length",
        )

        def _CylindricalVolume(field, data):
            r = data["index", "r"]
            dr = data["index", "dr"]
            vol = 0.5 * ((r + 0.5 * dr) ** 2 - (r - 0.5 * dr) ** 2)
            vol *= data["index", "dtheta"]
            vol *= data["index", "dz"]
            return vol

        registry.add_field(
            ("index", "cell_volume"),
            sampling_type="cell",
            function=_CylindricalVolume,
            units="code_length**3",
        )
        registry.alias(("index", "volume"), ("index", "cell_volume"))

        def _path_z(field, data):
            return data["index", "dz"]

        registry.add_field(
            ("index", "path_element_z"),
            sampling_type="cell",
            function=_path_z,
            units="code_length",
        )

    def pixelize(
        self,
        dimension,
        data_source,
        field,
        bounds,
        size,
        antialias=True,
        periodic=False,
        *,
        axes_transform=AxesTransform.DEFAULT,
    ):
        # Note that above, we set periodic by default to be *false*.  This is
        # because our pixelizers, at present, do not handle periodicity
        # correctly, and if you change the "width" of a cylindrical plot, it
        # double-counts in the edge buffers.  See, for instance, issue 1669.
        if axes_transform is AxesTransform.GEOMETRY_NATIVE:
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )

        name = self.axis_name[dimension]
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[name]

        if name in ("r", "theta"):
            if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        elif name == "z":
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            # This is admittedly a very hacky way to resolve a bug
            # it's very likely that the *right* fix would have to
            # be applied upstream of this function, *but* this case
            # has never worked properly so maybe it's still preferable to
            # not having a solution ?
            # see https://github.com/yt-project/yt/pull/3533
            bounds = (*bounds[2:4], *bounds[:2])
            return self._cyl_pixelize(data_source, field, bounds, size, antialias)
        else:
            # Pixelizing along a cylindrical surface is a bit tricky
            raise NotImplementedError

    def pixelize_line(self, field, start_point, end_point, npoints):
        raise NotImplementedError

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias):
        buff = np.full((size[1], size[0]), np.nan, dtype="f8")
        pixelize_cylinder(
            buff,
            data_source["px"],
            data_source["pdx"],
            data_source["py"],
            data_source["pdy"],
            data_source[field],
            bounds,
        )
        return buff

    _x_pairs = (("r", "theta"), ("z", "r"), ("theta", "r"))  # deprecated
    _y_pairs = (("r", "z"), ("z", "theta"), ("theta", "z"))  # deprecated

    _image_axis_name = None  # deprecated

    @property
    def image_axis_name(self):
        issue_deprecation_warning(
            "The image_axis_name property isn't used "
            "internally in yt anymore and is deprecated",
            since="4.2.0",
        )
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {
            self.axis_id["r"]: ("\\theta", "z"),
            self.axis_id["z"]: ("x", "y"),
            self.axis_id["theta"]: ("r", "z"),
        }
        for i in list(rv.keys()):
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].upper()] = rv[i]
        self._image_axis_name = rv
        return rv

    def convert_from_cartesian(self, coord):
        return cartesian_to_cylindrical(coord)

    def convert_to_cartesian(self, coord):
        return cylindrical_to_cartesian(coord)

    def convert_to_cylindrical(self, coord):
        return coord

    def convert_from_cylindrical(self, coord):
        return coord

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    @property
    def period(self):
        return np.array([0.0, 0.0, 2.0 * np.pi])

    @cached_property
    def _polar_bounds(self):
        return _get_polar_bounds(self, axes=("r", "theta"))

    def _get_display_center(self, center, axes_transform):
        if axes_transform is AxesTransform.GEOMETRY_NATIVE:
            return center
        else:
            return super()._get_display_center(center, axes_transform)

    def sanitize_center(self, center, axis, *, axes_transform=AxesTransform.DEFAULT):
        if axes_transform is AxesTransform.GEOMETRY_NATIVE:
            return super().sanitize_center(center, axis, axes_transform=axes_transform)

        name = self.axis_name[axis]
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[name]
        center, display_center = super().sanitize_center(
            center, axis, axes_transform=AxesTransform.GEOMETRY_NATIVE
        )
        display_center = [
            0.0 * display_center[0],
            0.0 * display_center[1],
            0.0 * display_center[2],
        ]
        r_ax = self.axis_id["r"]
        theta_ax = self.axis_id["theta"]
        z_ax = self.axis_id["z"]
        if name == "r":
            if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            display_center[theta_ax] = self.ds.domain_center[theta_ax]
            display_center[z_ax] = self.ds.domain_center[z_ax]
        elif name == "theta":
            if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            # use existing center value
            for idx in (r_ax, z_ax):
                display_center[idx] = center[idx]
        elif name == "z":
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xxmin, xxmax, yymin, yymax = self._polar_bounds
            xc = (xxmin + xxmax) / 2
            yc = (yymin + yymax) / 2
            display_center = (xc, yc, 0 * xc)
        else:
            RuntimeError(f"Unknown axis name {name!r} for cylindrical coordinates")
        return center, display_center

    def sanitize_width(
        self, axis, width, depth, *, axes_transform=AxesTransform.DEFAULT
    ):
        if axes_transform is AxesTransform.GEOMETRY_NATIVE:
            return super().sanitize_width(
                axis, width, depth, axes_transform=axes_transform
            )

        name = self.axis_name[axis]
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[name]

        r_ax, theta_ax, z_ax = (
            self.ds.coordinates.axis_id[ax] for ax in ("r", "theta", "z")
        )
        if width is not None:
            width = super().sanitize_width(
                axis, width, depth, axes_transform=axes_transform
            )
        # Note: regardless of axes, these are set up to give consistent plots
        # when plotted, which is not strictly a "right hand rule" for axes.
        elif name == "r":  # soup can label
            if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            width = [self.ds.domain_width[theta_ax], self.ds.domain_width[z_ax]]
        elif name == "theta":
            if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            width = [self.ds.domain_width[r_ax], self.ds.domain_width[z_ax]]
        elif name == "z":
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xxmin, xxmax, yymin, yymax = self._polar_bounds
            xw = xxmax - xxmin
            yw = yymax - yymin
            width = [xw, yw]
        return width

    def _get_plot_axes_default_properties(
        self, normal_axis_name: str, axes_transform: AxesTransform
    ) -> DefaultProperties:
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[normal_axis_name]

        if normal_axis_name == "r":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return dict(
                    x_axis_label=r"\theta",
                    y_axis_label="z",
                    x_axis_units="rad",
                    y_axis_units=None,
                )
            else:
                raise NotImplementedError
        elif normal_axis_name == "theta":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return dict(
                    x_axis_label="R",
                    y_axis_label="z",
                    x_axis_units=None,
                    y_axis_units=None,
                )
            else:
                raise NotImplementedError
        elif normal_axis_name == "z":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return dict(
                    x_axis_label="R",
                    y_axis_label=r"\theta",
                    x_axis_units=None,
                    y_axis_units="rad",
                )
            elif axes_transform is AxesTransform.POLAR:
                return dict(
                    x_axis_label="x",
                    y_axis_label="y",
                    x_axis_units=None,
                    y_axis_units=None,
                )
            else:
                raise NotImplementedError
        else:
            raise ValueError(f"Unknown axis name {normal_axis_name!r}")
