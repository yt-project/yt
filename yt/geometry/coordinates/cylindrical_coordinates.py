import numpy as np

from yt.utilities.lib.pixelization_routines import pixelize_cartesian, pixelize_cylinder

from .coordinate_handler import (
    CoordinateHandler,
    _get_coord_fields,
    _unknown_coord,
    cartesian_to_cylindrical,
    cylindrical_to_cartesian,
)

#
# Cylindrical fields
#


class CylindricalCoordinateHandler(CoordinateHandler):
    name = "cylindrical"

    def __init__(self, ds, ordering=("r", "z", "theta")):
        super(CylindricalCoordinateHandler, self).__init__(ds, ordering)
        self.image_units = {}
        self.image_units[self.axis_id["r"]] = ("rad", None)
        self.image_units[self.axis_id["theta"]] = (None, None)
        self.image_units[self.axis_id["z"]] = (None, None)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(
            ("index", "dx"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "dy"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "x"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "y"), sampling_type="cell", function=_unknown_coord
        )
        f1, f2 = _get_coord_fields(self.axis_id["r"])
        registry.add_field(
            ("index", "dr"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="code_length",
        )

        registry.add_field(
            ("index", "r"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="code_length",
        )

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

        f1, f2 = _get_coord_fields(self.axis_id["theta"], "")
        registry.add_field(
            ("index", "dtheta"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="",
        )

        registry.add_field(
            ("index", "theta"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="",
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

        def _path_r(field, data):
            return data["index", "dr"]

        registry.add_field(
            ("index", "path_element_r"),
            sampling_type="cell",
            function=_path_r,
            units="code_length",
        )

        def _path_theta(field, data):
            # Note: this already assumes cell-centered
            return data["index", "r"] * data["index", "dtheta"]

        registry.add_field(
            ("index", "path_element_theta"),
            sampling_type="cell",
            function=_path_theta,
            units="code_length",
        )

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
    ):
        # Note that above, we set periodic by default to be *false*.  This is
        # because our pixelizers, at present, do not handle periodicity
        # correctly, and if you change the "width" of a cylindrical plot, it
        # double-counts in the edge buffers.  See, for instance, issue 1669.
        ax_name = self.axis_name[dimension]
        if ax_name in ("r", "theta"):
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        elif ax_name == "z":
            return self._cyl_pixelize(data_source, field, bounds, size, antialias)
        else:
            # Pixelizing along a cylindrical surface is a bit tricky
            raise NotImplementedError

    def _ortho_pixelize(
        self, data_source, field, bounds, size, antialias, dim, periodic
    ):
        period = self.period[:2].copy()  # dummy here
        period[0] = self.period[self.x_axis[dim]]
        period[1] = self.period[self.y_axis[dim]]
        if hasattr(period, "in_units"):
            period = period.in_units("code_length").d
        buff = np.zeros(size, dtype="f8")
        pixelize_cartesian(
            buff,
            data_source["px"],
            data_source["py"],
            data_source["pdx"],
            data_source["pdy"],
            data_source[field],
            bounds,
            int(antialias),
            period,
            int(periodic),
        )
        return buff

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
        self.sanitize_buffer_fill_values(buff)
        return buff

    _x_pairs = (("r", "theta"), ("z", "r"), ("theta", "r"))
    _y_pairs = (("r", "z"), ("z", "theta"), ("theta", "z"))

    _image_axis_name = None

    @property
    def image_axis_name(self):
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {
            self.axis_id["r"]: ("theta", "z"),
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

    def sanitize_center(self, center, axis):
        center, display_center = super(
            CylindricalCoordinateHandler, self
        ).sanitize_center(center, axis)
        display_center = [
            0.0 * display_center[0],
            0.0 * display_center[1],
            0.0 * display_center[2],
        ]
        ax_name = self.axis_name[axis]
        r_ax = self.axis_id["r"]
        theta_ax = self.axis_id["theta"]
        z_ax = self.axis_id["z"]
        if ax_name == "r":
            display_center[theta_ax] = self.ds.domain_center[theta_ax]
            display_center[z_ax] = self.ds.domain_center[z_ax]
        elif ax_name == "theta":
            # use existing center value
            for idx in (r_ax, z_ax):
                display_center[idx] = center[idx]
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        r_ax, theta_ax, z_ax = (
            self.ds.coordinates.axis_id[ax] for ax in ("r", "theta", "z")
        )
        if width is not None:
            width = super(CylindricalCoordinateHandler, self).sanitize_width(
                axis, width, depth
            )
        # Note: regardless of axes, these are set up to give consistent plots
        # when plotted, which is not strictly a "right hand rule" for axes.
        elif name == "r":  # soup can label
            width = [2.0 * np.pi * self.ds.domain_width.uq, self.ds.domain_width[z_ax]]
        elif name == "theta":
            width = [self.ds.domain_width[r_ax], self.ds.domain_width[z_ax]]
        elif name == "z":
            width = [
                2.0 * self.ds.domain_right_edge[r_ax],
                2.0 * self.ds.domain_right_edge[r_ax],
            ]
        return width
