import numpy as np

from yt.utilities.lib.pixelization_routines import pixelize_aitoff, pixelize_cylinder

from .coordinate_handler import CoordinateHandler, _get_coord_fields, _unknown_coord


class SphericalCoordinateHandler(CoordinateHandler):
    name = "spherical"

    def __init__(self, ds, ordering=("r", "theta", "phi")):
        super().__init__(ds, ordering)
        # Generate
        self.image_units = {}
        self.image_units[self.axis_id["r"]] = ("rad", "rad")
        self.image_units[self.axis_id["theta"]] = (None, None)
        self.image_units[self.axis_id["phi"]] = (None, None)

    def setup_fields(self, registry):
        # return the fields for r, z, theta
        registry.add_field(
            ("index", "dx"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "dy"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "dz"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "x"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "y"), sampling_type="cell", function=_unknown_coord
        )

        registry.add_field(
            ("index", "z"), sampling_type="cell", function=_unknown_coord
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

        f1, f2 = _get_coord_fields(self.axis_id["phi"], "")
        registry.add_field(
            ("index", "dphi"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="",
        )

        registry.add_field(
            ("index", "phi"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="",
        )

        def _SphericalVolume(field, data):
            # Here we compute the spherical volume element exactly
            r = data["index", "r"]
            dr = data["index", "dr"]
            theta = data["index", "theta"]
            dtheta = data["index", "dtheta"]
            vol = ((r + 0.5 * dr) ** 3 - (r - 0.5 * dr) ** 3) / 3.0
            vol *= np.cos(theta - 0.5 * dtheta) - np.cos(theta + 0.5 * dtheta)
            vol *= data["index", "dphi"]
            return vol

        registry.add_field(
            ("index", "cell_volume"),
            sampling_type="cell",
            function=_SphericalVolume,
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

        def _path_phi(field, data):
            # Note: this already assumes cell-centered
            return (
                data["index", "r"]
                * data["index", "dphi"]
                * np.sin(data["index", "theta"])
            )

        registry.add_field(
            ("index", "path_element_phi"),
            sampling_type="cell",
            function=_path_phi,
            units="code_length",
        )

    def pixelize(
        self, dimension, data_source, field, bounds, size, antialias=True, periodic=True
    ):
        self.period
        name = self.axis_name[dimension]
        if name == "r":
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        elif name in ("theta", "phi"):
            return self._cyl_pixelize(
                data_source, field, bounds, size, antialias, dimension
            )
        else:
            raise NotImplementedError

    def _ortho_pixelize(
        self, data_source, field, bounds, size, antialias, dim, periodic
    ):
        buff = pixelize_aitoff(
            data_source["py"],
            data_source["pdy"],
            data_source["px"],
            data_source["pdx"],
            size,
            data_source[field],
            None,
            None,
            theta_offset=0,
            phi_offset=0,
        ).transpose()
        return buff

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias, dimension):
        name = self.axis_name[dimension]
        buff = np.full((size[1], size[0]), np.nan, dtype="f8")
        if name == "theta":
            pixelize_cylinder(
                buff,
                data_source["px"],
                data_source["pdx"],
                data_source["py"],
                data_source["pdy"],
                data_source[field],
                bounds,
            )
        elif name == "phi":
            # Note that we feed in buff.T here
            pixelize_cylinder(
                buff.T,
                data_source["px"],
                data_source["pdx"],
                data_source["py"],
                data_source["pdy"],
                data_source[field],
                bounds,
            )
        else:
            raise RuntimeError
        self.sanitize_buffer_fill_values(buff)
        return buff

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            ri = self.axis_id["r"]
            thetai = self.axis_id["theta"]
            phii = self.axis_id["phi"]
            r = coord[:, ri]
            theta = coord[:, thetai]
            phi = coord[:, phii]
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:, ri] = np.cos(phi) * np.sin(theta) * r
            nc[:, thetai] = np.sin(phi) * np.sin(theta) * r
            nc[:, phii] = np.cos(theta) * r
        else:
            r, theta, phi = coord
            nc = (
                np.cos(phi) * np.sin(theta) * r,
                np.sin(phi) * np.sin(theta) * r,
                np.cos(theta) * r,
            )
        return nc

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    def convert_to_spherical(self, coord):
        raise NotImplementedError

    def convert_from_spherical(self, coord):
        raise NotImplementedError

    _image_axis_name = None

    @property
    def image_axis_name(self):
        if self._image_axis_name is not None:
            return self._image_axis_name
        # This is the x and y axes labels that get displayed.  For
        # non-Cartesian coordinates, we usually want to override these for
        # Cartesian coordinates, since we transform them.
        rv = {
            self.axis_id["r"]: ("theta", "phi"),
            self.axis_id["theta"]: ("x / \\sin(\\theta)", "y / \\sin(\\theta)"),
            self.axis_id["phi"]: ("R", "z"),
        }
        for i in list(rv.keys()):
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].capitalize()] = rv[i]
        self._image_axis_name = rv
        return rv

    _x_pairs = (("r", "theta"), ("theta", "r"), ("phi", "r"))
    _y_pairs = (("r", "phi"), ("theta", "phi"), ("phi", "theta"))

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super().sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == "r":
            display_center = center
        elif name == "theta":
            display_center = (
                0.0 * display_center[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            )
        elif name == "phi":
            display_center = [
                self.ds.domain_width[0] / 2.0 + self.ds.domain_left_edge[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            ]
            ri = self.axis_id["r"]
            c = self.ds.domain_width[ri] / 2.0 + self.ds.domain_left_edge[ri]
            display_center[ri] = c
            display_center = tuple(display_center)
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super().sanitize_width(axis, width, depth)
        elif name == "r":
            width = [
                self.ds.domain_width[self.x_axis["r"]],
                self.ds.domain_width[self.y_axis["r"]],
            ]
        elif name == "theta":
            ri = self.axis_id["r"]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            width = [2.0 * self.ds.domain_width[ri], 2.0 * self.ds.domain_width[ri]]
        elif name == "phi":
            ri = self.axis_id["r"]
            width = [self.ds.domain_right_edge[ri], 2.0 * self.ds.domain_width[ri]]
        return width

    def _sanity_check(self):
        """This prints out a handful of diagnostics that help verify the
        dataset is well formed."""
        # We just check a few things here.
        dd = self.ds.all_data()
        r0 = self.ds.domain_left_edge[self.axis_id["r"]]
        r1 = self.ds.domain_right_edge[self.axis_id["r"]]
        v1 = 4.0 * np.pi / 3.0 * (r1 ** 3 - r0 ** 3)
        print(f"Total volume should be 4*pi*r**3 = {v1:0.16e}")
        v2 = dd.quantities.total_quantity("cell_volume")
        print(f"Actual volume is                   {v2:0.16e}")
        print(f"Relative difference: {np.abs(v2 - v1) / (v2 + v1):0.16e}")
