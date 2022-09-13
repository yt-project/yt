import numpy as np

from yt.utilities.lib.pixelization_routines import pixelize_cartesian, pixelize_cylinder

from .coordinate_handler import (
    CoordinateHandler,
    _get_coord_fields,
    _setup_dummy_cartesian_coords_and_widths,
)


class GeographicCoordinateHandler(CoordinateHandler):
    radial_axis = "altitude"
    name = "geographic"

    def __init__(self, ds, ordering=None):
        if not ordering:
            ordering = ("latitude", "longitude", self.radial_axis)
        super().__init__(ds, ordering)
        self.image_units = {}
        self.image_units[self.axis_id["latitude"]] = (None, None)
        self.image_units[self.axis_id["longitude"]] = (None, None)
        self.image_units[self.axis_id[self.radial_axis]] = ("deg", "deg")

    def setup_fields(self, registry):
        # Missing implementation for x, y and z coordinates.
        _setup_dummy_cartesian_coords_and_widths(registry, axes=("x", "y", "z"))

        f1, f2 = _get_coord_fields(self.axis_id["latitude"], "")
        registry.add_field(
            ("index", "dlatitude"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="",
        )

        registry.add_field(
            ("index", "latitude"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="",
        )

        f1, f2 = _get_coord_fields(self.axis_id["longitude"], "")
        registry.add_field(
            ("index", "dlongitude"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="",
        )

        registry.add_field(
            ("index", "longitude"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="",
        )

        f1, f2 = _get_coord_fields(self.axis_id[self.radial_axis])
        registry.add_field(
            ("index", f"d{self.radial_axis}"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="code_length",
        )

        registry.add_field(
            ("index", self.radial_axis),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="code_length",
        )

        def _SphericalVolume(field, data):
            # We can use the transformed coordinates here.
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

        def _path_radial_axis(field, data):
            return data["index", f"d{self.radial_axis}"]

        registry.add_field(
            ("index", f"path_element_{self.radial_axis}"),
            sampling_type="cell",
            function=_path_radial_axis,
            units="code_length",
        )

        def _path_latitude(field, data):
            # We use r here explicitly
            return data["index", "r"] * data["index", "dlatitude"] * np.pi / 180.0

        registry.add_field(
            ("index", "path_element_latitude"),
            sampling_type="cell",
            function=_path_latitude,
            units="code_length",
        )

        def _path_longitude(field, data):
            # We use r here explicitly
            return (
                data["index", "r"]
                * data["index", "dlongitude"]
                * np.pi
                / 180.0
                * np.sin((data["index", "latitude"] + 90.0) * np.pi / 180.0)
            )

        registry.add_field(
            ("index", "path_element_longitude"),
            sampling_type="cell",
            function=_path_longitude,
            units="code_length",
        )

        def _latitude_to_theta(field, data):
            # latitude runs from -90 to 90
            return (data[("index", "latitude")] + 90) * np.pi / 180.0

        registry.add_field(
            ("index", "theta"),
            sampling_type="cell",
            function=_latitude_to_theta,
            units="",
        )

        def _dlatitude_to_dtheta(field, data):
            return data[("index", "dlatitude")] * np.pi / 180.0

        registry.add_field(
            ("index", "dtheta"),
            sampling_type="cell",
            function=_dlatitude_to_dtheta,
            units="",
        )

        def _longitude_to_phi(field, data):
            # longitude runs from -180 to 180
            return (data[("index", "longitude")] + 180) * np.pi / 180.0

        registry.add_field(
            ("index", "phi"), sampling_type="cell", function=_longitude_to_phi, units=""
        )

        def _dlongitude_to_dphi(field, data):
            return data[("index", "dlongitude")] * np.pi / 180.0

        registry.add_field(
            ("index", "dphi"),
            sampling_type="cell",
            function=_dlongitude_to_dphi,
            units="",
        )

        self._setup_radial_fields(registry)

    def _setup_radial_fields(self, registry):
        # This stays here because we don't want to risk the field detector not
        # properly getting the data_source, etc.
        def _altitude_to_radius(field, data):
            surface_height = data.get_field_parameter("surface_height")
            if surface_height is None:
                if hasattr(data.ds, "surface_height"):
                    surface_height = data.ds.surface_height
                else:
                    surface_height = data.ds.quan(0.0, "code_length")
            return data[("index", "altitude")] + surface_height

        registry.add_field(
            ("index", "r"),
            sampling_type="cell",
            function=_altitude_to_radius,
            units="code_length",
        )
        registry.alias(("index", "dr"), ("index", "daltitude"))

    def _retrieve_radial_offset(self, data_source=None):
        # This returns the factor by which the radial field is multiplied and
        # the scalar its offset by.  Typically the "factor" will *only* be
        # either 1.0 or -1.0.  The order will be factor * r + offset.
        # Altitude is the radius from the central zone minus the radius of the
        # surface.  Depth to radius is negative value of depth plus the
        # outermost radius.
        surface_height = None
        if data_source is not None:
            surface_height = data_source.get_field_parameter("surface_height")
        if surface_height is None:
            if hasattr(self.ds, "surface_height"):
                surface_height = self.ds.surface_height
            else:
                surface_height = self.ds.quan(0.0, "code_length")
        return surface_height, 1.0

    def pixelize(
        self, dimension, data_source, field, bounds, size, antialias=True, periodic=True
    ):
        if self.axis_name[dimension] in ("latitude", "longitude"):
            return self._cyl_pixelize(
                data_source, field, bounds, size, antialias, dimension
            )
        elif self.axis_name[dimension] == self.radial_axis:
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        else:
            raise NotImplementedError

    def pixelize_line(self, field, start_point, end_point, npoints):
        raise NotImplementedError

    def _ortho_pixelize(
        self, data_source, field, bounds, size, antialias, dimension, periodic
    ):

        period = self.period[:2].copy()
        period[0] = self.period[self.x_axis[dimension]]
        period[1] = self.period[self.y_axis[dimension]]
        if hasattr(period, "in_units"):
            period = period.in_units("code_length").d

        # For a radial axis, px will correspond to longitude and py will
        # correspond to latitude.
        px = data_source["px"]
        pdx = data_source["pdx"]
        py = data_source["py"]
        pdy = data_source["pdy"]
        buff = np.full((size[1], size[0]), np.nan, dtype="float64")
        pixelize_cartesian(
            buff,
            px,
            py,
            pdx,
            pdy,
            data_source[field],
            bounds,
            int(antialias),
            period,
            int(periodic),
        )
        return buff

    def _cyl_pixelize(self, data_source, field, bounds, size, antialias, dimension):
        offset, factor = self._retrieve_radial_offset(data_source)
        r = factor * data_source["py"] + offset
        # Because of the axis-ordering, dimensions 0 and 1 both have r as py
        # and the angular coordinate as px.  But we need to figure out how to
        # convert our coordinate back to an actual angle, based on which
        # dimension we're in.
        pdx = data_source["pdx"].d * np.pi / 180
        if self.axis_name[self.x_axis[dimension]] == "latitude":
            px = (data_source["px"].d + 90) * np.pi / 180
            do_transpose = True
        elif self.axis_name[self.x_axis[dimension]] == "longitude":
            px = (data_source["px"].d + 180) * np.pi / 180
            do_transpose = False
        else:
            # We should never get here!
            raise NotImplementedError
        buff = np.full((size[1], size[0]), np.nan, dtype="f8")
        pixelize_cylinder(
            buff, r, data_source["pdy"], px, pdx, data_source[field], bounds
        )
        if do_transpose:
            buff = buff.transpose()
        self.sanitize_buffer_fill_values(buff)
        return buff

    def convert_from_cartesian(self, coord):
        raise NotImplementedError

    def convert_to_cartesian(self, coord):
        offset, factor = self._retrieve_radial_offset()
        if isinstance(coord, np.ndarray) and len(coord.shape) > 1:
            rad = self.axis_id[self.radial_axis]
            lon = self.axis_id["longitude"]
            lat = self.axis_id["latitude"]
            r = factor * coord[:, rad] + offset
            theta = coord[:, lon] * np.pi / 180
            phi = coord[:, lat] * np.pi / 180
            nc = np.zeros_like(coord)
            # r, theta, phi
            nc[:, lat] = np.cos(phi) * np.sin(theta) * r
            nc[:, lon] = np.sin(phi) * np.sin(theta) * r
            nc[:, rad] = np.cos(theta) * r
        else:
            a, b, c = coord
            theta = b * np.pi / 180
            phi = a * np.pi / 180
            r = factor * c + offset
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
            self.axis_id["latitude"]: (
                "x / \\sin(\\mathrm{latitude})",
                "y / \\sin(\\mathrm{latitude})",
            ),
            self.axis_id["longitude"]: ("R", "z"),
            self.axis_id[self.radial_axis]: ("longitude", "latitude"),
        }
        for i in list(rv.keys()):
            rv[self.axis_name[i]] = rv[i]
            rv[self.axis_name[i].capitalize()] = rv[i]
        self._image_axis_name = rv
        return rv

    _x_pairs = (
        ("latitude", "longitude"),
        ("longitude", "latitude"),
        ("altitude", "longitude"),
    )

    _y_pairs = (
        ("latitude", "altitude"),
        ("longitude", "altitude"),
        ("altitude", "latitude"),
    )

    _data_projection = None

    @property
    def data_projection(self):
        # this will control the default projection to use when displaying data
        if self._data_projection is not None:
            return self._data_projection
        dpj = {}
        for ax in self.axis_order:
            if ax == self.radial_axis:
                dpj[ax] = "Mollweide"
            else:
                dpj[ax] = None
        self._data_projection = dpj
        return dpj

    _data_transform = None

    @property
    def data_transform(self):
        # this is the coordinate system on which the data is defined (the crs).
        if self._data_transform is not None:
            return self._data_transform
        dtx = {}
        for ax in self.axis_order:
            if ax == self.radial_axis:
                dtx[ax] = "PlateCarree"
            else:
                dtx[ax] = None
        self._data_transform = dtx
        return dtx

    @property
    def period(self):
        return self.ds.domain_width

    def sanitize_center(self, center, axis):
        center, display_center = super().sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == self.radial_axis:
            display_center = center
        elif name == "latitude":
            display_center = (
                0.0 * display_center[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            )
        elif name == "longitude":
            ri = self.axis_id[self.radial_axis]
            c = (self.ds.domain_right_edge[ri] + self.ds.domain_left_edge[ri]) / 2.0
            display_center = [
                0.0 * display_center[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            ]
            display_center[self.axis_id["latitude"]] = c
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super().sanitize_width(axis, width, depth)
        elif name == self.radial_axis:
            rax = self.radial_axis
            width = [
                self.ds.domain_width[self.x_axis[rax]],
                self.ds.domain_width[self.y_axis[rax]],
            ]
        elif name == "latitude":
            ri = self.axis_id[self.radial_axis]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            width = [2.0 * self.ds.domain_width[ri], 2.0 * self.ds.domain_width[ri]]
        elif name == "longitude":
            ri = self.axis_id[self.radial_axis]
            width = [self.ds.domain_width[ri], 2.0 * self.ds.domain_width[ri]]
        return width


class InternalGeographicCoordinateHandler(GeographicCoordinateHandler):
    radial_axis = "depth"
    name = "internal_geographic"

    def _setup_radial_fields(self, registry):
        # Altitude is the radius from the central zone minus the radius of the
        # surface.
        def _depth_to_radius(field, data):
            outer_radius = data.get_field_parameter("outer_radius")
            if outer_radius is None:
                if hasattr(data.ds, "outer_radius"):
                    outer_radius = data.ds.outer_radius
                else:
                    # Otherwise, we assume that the depth goes to full depth,
                    # so we can look at the domain right edge in depth.
                    rax = self.axis_id[self.radial_axis]
                    outer_radius = data.ds.domain_right_edge[rax]
            return -1.0 * data[("index", "depth")] + outer_radius

        registry.add_field(
            ("index", "r"),
            sampling_type="cell",
            function=_depth_to_radius,
            units="code_length",
        )
        registry.alias(("index", "dr"), ("index", "ddepth"))

    def _retrieve_radial_offset(self, data_source=None):
        # Depth means switching sign and adding to full radius
        outer_radius = None
        if data_source is not None:
            outer_radius = data_source.get_field_parameter("outer_radius")
        if outer_radius is None:
            if hasattr(self.ds, "outer_radius"):
                outer_radius = self.ds.outer_radius
            else:
                # Otherwise, we assume that the depth goes to full depth,
                # so we can look at the domain right edge in depth.
                rax = self.axis_id[self.radial_axis]
                outer_radius = self.ds.domain_right_edge[rax]
        return outer_radius, -1.0

    _x_pairs = (
        ("latitude", "longitude"),
        ("longitude", "latitude"),
        ("depth", "longitude"),
    )

    _y_pairs = (("latitude", "depth"), ("longitude", "depth"), ("depth", "latitude"))

    def sanitize_center(self, center, axis):
        center, display_center = super(
            GeographicCoordinateHandler, self
        ).sanitize_center(center, axis)
        name = self.axis_name[axis]
        if name == self.radial_axis:
            display_center = center
        elif name == "latitude":
            display_center = (
                0.0 * display_center[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            )
        elif name == "longitude":
            ri = self.axis_id[self.radial_axis]
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            display_center = [
                0.0 * display_center[0],
                0.0 * display_center[1],
                0.0 * display_center[2],
            ]
            display_center[self.axis_id["latitude"]] = outermost / 2.0
        return center, display_center

    def sanitize_width(self, axis, width, depth):
        name = self.axis_name[axis]
        if width is not None:
            width = super(GeographicCoordinateHandler, self).sanitize_width(
                axis, width, depth
            )
        elif name == self.radial_axis:
            rax = self.radial_axis
            width = [
                self.ds.domain_width[self.x_axis[rax]],
                self.ds.domain_width[self.y_axis[rax]],
            ]
        elif name == "latitude":
            ri = self.axis_id[self.radial_axis]
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            width = [2.0 * outermost, 2.0 * outermost]
        elif name == "longitude":
            ri = self.axis_id[self.radial_axis]
            offset, factor = self._retrieve_radial_offset()
            outermost = factor * self.ds.domain_left_edge[ri] + offset
            width = [outermost, 2.0 * outermost]
        return width
