from functools import cached_property
from typing import Dict

import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.utilities.lib.pixelization_routines import pixelize_aitoff, pixelize_cylinder

from ._axes_transforms import AxesTransform
from .coordinate_handler import (
    CoordinateHandler,
    DefaultProperties,
    _get_coord_fields,
    _get_polar_bounds,
    _setup_dummy_cartesian_coords_and_widths,
    _setup_polar_coordinates,
)


class SphericalCoordinateHandler(CoordinateHandler):
    name = "spherical"
    _default_axis_order = ("r", "theta", "phi")
    _default_axes_transforms: Dict[str, AxesTransform] = {
        "r": AxesTransform.AITOFF_HAMMER,
        "theta": AxesTransform.POLAR,
        "phi": AxesTransform.POLAR,
    }

    def __init__(self, ds, ordering=None):
        super().__init__(ds, ordering)

    def setup_fields(self, registry):
        # Missing implementation for x, y and z coordinates.
        _setup_dummy_cartesian_coords_and_widths(registry, axes=("x", "y", "z"))
        _setup_polar_coordinates(registry, self.axis_id)

        f1, f2 = _get_coord_fields(self.axis_id["phi"], "dimensionless")
        registry.add_field(
            ("index", "dphi"),
            sampling_type="cell",
            function=f1,
            display_field=False,
            units="dimensionless",
        )

        registry.add_field(
            ("index", "phi"),
            sampling_type="cell",
            function=f2,
            display_field=False,
            units="dimensionless",
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
        self,
        dimension,
        data_source,
        field,
        bounds,
        size,
        antialias=True,
        periodic=True,
        *,
        axes_transform=AxesTransform.DEFAULT,
    ):
        if axes_transform is AxesTransform.GEOMETRY_NATIVE:
            return self._ortho_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )

        self.period
        name = self.axis_name[dimension]
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[name]

        if name == "r":
            if axes_transform is not AxesTransform.AITOFF_HAMMER:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            return self._aitoff_hammer_pixelize(
                data_source, field, bounds, size, antialias, dimension, periodic
            )
        elif name in ("theta", "phi"):
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            if name == "theta":
                # This is admittedly a very hacky way to resolve a bug
                # it's very likely that the *right* fix would have to
                # be applied upstream of this function, *but* this case
                # has never worked properly so maybe it's still preferable to
                # not having a solution ?
                # see https://github.com/yt-project/yt/pull/3533
                bounds = (*bounds[2:4], *bounds[:2])
            return self._cyl_pixelize(
                data_source, field, bounds, size, antialias, dimension
            )
        else:
            raise NotImplementedError

    def pixelize_line(self, field, start_point, end_point, npoints):
        raise NotImplementedError

    def _aitoff_hammer_pixelize(
        self, data_source, field, bounds, size, antialias, dim, periodic
    ):
        # use Aitoff projection
        # http://paulbourke.net/geometry/transformationprojection/
        bounds = tuple(_.ndview for _ in self._aitoff_bounds)
        buff = pixelize_aitoff(
            azimuth=data_source["py"],
            dazimuth=data_source["pdy"],
            colatitude=data_source["px"],
            dcolatitude=data_source["pdx"],
            buff_size=size,
            field=data_source[field],
            bounds=bounds,
            input_img=None,
            azimuth_offset=0,
            colatitude_offset=0,
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
            self.axis_id["r"]: (
                # these are the Hammer-Aitoff normalized coordinates
                # conventions:
                #   - theta is the colatitude, from 0 to PI
                #   - bartheta is the latitude, from -PI/2 to +PI/2 (bartheta = PI/2 - theta)
                #   - phi is the azimuth, from 0 to 2PI
                #   - lambda is the longitude, from -PI to PI (lambda = phi - PI)
                r"\frac{2\cos(\mathrm{\bar{\theta}})\sin(\lambda/2)}{\sqrt{1 + \cos(\bar{\theta}) \cos(\lambda/2)}}",
                r"\frac{sin(\bar{\theta})}{\sqrt{1 + \cos(\bar{\theta}) \cos(\lambda/2)}}",
                "\\theta",
            ),
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

    def _get_plot_axes_default_properties(
        self, normal_axis_name: str, axes_transform: AxesTransform
    ) -> DefaultProperties:
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[normal_axis_name]

        if normal_axis_name == "r":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return {
                    "x_axis_label": r"\theta",
                    "y_axis_label": r"\phi",
                    "x_axis_units": "rad",
                    "y_axis_units": "rad",
                }
            elif axes_transform is AxesTransform.AITOFF_HAMMER:
                return {
                    "x_axis_label": r"\frac{2\cos(\mathrm{\bar{\theta}})\sin(\lambda/2)}{\sqrt{1 + \cos(\bar{\theta}) \cos(\lambda/2)}}",
                    "y_axis_label": r"\frac{sin(\bar{\theta})}{\sqrt{1 + \cos(\bar{\theta}) \cos(\lambda/2)}}",
                    "x_axis_units": "dimensionless",
                    "y_axis_units": "dimensionless",
                }
            else:
                raise NotImplementedError
        elif normal_axis_name == "theta":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return {
                    "x_axis_label": "r",
                    "y_axis_label": r"\phi",
                    "x_axis_units": None,
                    "y_axis_units": "rad",
                }
            elif axes_transform is AxesTransform.POLAR:
                return {
                    "x_axis_label": r"x / \sin(\theta)",
                    "y_axis_label": r"y / \sin(\theta)",
                    "x_axis_units": None,
                    "y_axis_units": None,
                }
            else:
                raise NotImplementedError
        elif normal_axis_name == "phi":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                return {
                    "x_axis_label": "r",
                    "y_axis_label": r"\theta",
                    "x_axis_units": None,
                    "y_axis_units": "rad",
                }
            elif axes_transform is AxesTransform.POLAR:
                return {
                    "x_axis_label": "R",
                    "y_axis_label": "z",
                    "x_axis_units": None,
                    "y_axis_units": None,
                }
            else:
                raise NotImplementedError
        else:
            raise ValueError(f"Unknown axis name {normal_axis_name!r}")

    @property
    def period(self):
        return self.ds.domain_width

    @cached_property
    def _poloidal_bounds(self):
        rmin, rmax, thetamin, thetamax, _, _ = self._r_theta_phi_bounds
        corners = [
            (rmin, thetamin),
            (rmin, thetamax),
            (rmax, thetamin),
            (rmax, thetamax),
        ]

        def to_poloidal_plane(r, theta):
            # take a r, theta position and return
            # cylindrical R, z coordinates
            R = r * np.sin(theta)
            z = r * np.cos(theta)
            return R, z

        cylindrical_corner_coords = [to_poloidal_plane(*corner) for corner in corners]

        thetamin = thetamin.d
        thetamax = thetamax.d

        Rmin = min(R for R, z in cylindrical_corner_coords)

        if thetamin <= np.pi / 2 <= thetamax:
            Rmax = rmax
        else:
            Rmax = max(R for R, z in cylindrical_corner_coords)

        zmin = min(z for R, z in cylindrical_corner_coords)
        zmax = max(z for R, z in cylindrical_corner_coords)

        return Rmin, Rmax, zmin, zmax

    @cached_property
    def _conic_bounds(self):
        return _get_polar_bounds(self, axes=("r", "phi"))

    @cached_property
    def _r_theta_phi_bounds(self):
        # radius
        ri = self.axis_id["r"]
        rmin = self.ds.domain_left_edge[ri]
        rmax = self.ds.domain_right_edge[ri]

        # colatitude
        ti = self.axis_id["theta"]
        thetamin = self.ds.domain_left_edge[ti]
        thetamax = self.ds.domain_right_edge[ti]

        # azimuth
        pi = self.axis_id["phi"]
        phimin = self.ds.domain_left_edge[pi]
        phimax = self.ds.domain_right_edge[pi]

        return rmin, rmax, thetamin, thetamax, phimin, phimax

    @cached_property
    def _aitoff_bounds(self):
        # at the time of writing this function, yt's support for curvilinear
        # coordinates is a bit hacky, as many components of the system still
        # expect to receive coordinates with a length dimension. Ultimately
        # this is not needed but calls for a large refactor.
        ONE = self.ds.quan(1, "code_length")

        _, _, thetamin, thetamax, phimin, phimax = self._r_theta_phi_bounds

        # latitude
        latmin = ONE * np.pi / 2 - thetamax
        latmax = ONE * np.pi / 2 - thetamin

        # longitude
        lonmin = phimin - ONE * np.pi
        lonmax = phimax - ONE * np.pi

        corners = [
            (latmin, lonmin),
            (latmin, lonmax),
            (latmax, lonmin),
            (latmax, lonmax),
        ]

        def aitoff_z(latitude, longitude):
            return np.sqrt(1 + np.cos(latitude) * np.cos(longitude / 2))

        def aitoff_x(latitude, longitude):
            return (
                2
                * np.cos(latitude)
                * np.sin(longitude / 2)
                / aitoff_z(latitude, longitude)
            )

        def aitoff_y(latitude, longitude):
            return np.sin(latitude) / aitoff_z(latitude, longitude)

        def to_aitoff_plane(latitude, longitude):
            return aitoff_x(latitude, longitude), aitoff_y(latitude, longitude)

        aitoff_corner_coords = [to_aitoff_plane(*corner) for corner in corners]

        xmin = ONE * min(x for x, y in aitoff_corner_coords)
        xmax = ONE * max(x for x, y in aitoff_corner_coords)

        # theta is the colatitude
        # What this branch is meant to do is check whether the equator (latitude = 0)
        # is included in the domain.

        if latmin < 0 < latmax:
            xmin = min(xmin, ONE * aitoff_x(0, lonmin))
            xmax = max(xmax, ONE * aitoff_x(0, lonmax))

        # the y direction is more straightforward because aitoff-projected parallels (y)
        # draw a convex shape, while aitoff-projected meridians (x) draw a concave shape
        ymin = ONE * min(y for x, y in aitoff_corner_coords)
        ymax = ONE * max(y for x, y in aitoff_corner_coords)

        return xmin, xmax, ymin, ymax

    def sanitize_center(self, center, axis, *, axes_transform=AxesTransform.DEFAULT):
        center, display_center = super().sanitize_center(center, axis)
        name = self.axis_name[axis]
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = self._default_axes_transforms[name]

        if name == "r":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                _, _, xxmin, xxmax, yymin, yymax = self._r_theta_phi_bounds
            elif axes_transform is AxesTransform.AITOFF_HAMMER:
                xxmin, xxmax, yymin, yymax = self._aitoff_bounds
            else:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xc = (xxmin + xxmax) / 2
            yc = (yymin + yymax) / 2
            display_center = (0 * xc, xc, yc)

        elif name == "theta":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                xxmin, xxmax, _, _, yymin, yymax = self._r_theta_phi_bounds
            elif axes_transform is AxesTransform.POLAR:
                xxmin, xxmax, yymin, yymax = self._conic_bounds
            else:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xc = (xxmin + xxmax) / 2
            yc = (yymin + yymax) / 2
            display_center = (xc, 0 * xc, yc)
        elif name == "phi":
            if axes_transform is AxesTransform.GEOMETRY_NATIVE:
                xxmin, xxmax, yymin, yymax, _, _ = self._r_theta_phi_bounds
            elif axes_transform is AxesTransform.POLAR:
                xxmin, xxmax, yymin, yymax = self._poloidal_bounds
            else:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xc = (xxmin + xxmax) / 2
            yc = (yymin + yymax) / 2
            display_center = (xc, yc)
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

        if width is not None:
            width = super().sanitize_width(axis, width, depth)
        elif name == "r":
            if axes_transform is not AxesTransform.AITOFF_HAMMER:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            xxmin, xxmax, yymin, yymax = self._aitoff_bounds
            xw = xxmax - xxmin
            yw = yymax - yymin
            width = [xw, yw]
        elif name == "theta":
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            # Remember, in spherical coordinates when we cut in theta,
            # we create a conic section
            xxmin, xxmax, yymin, yymax = self._conic_bounds
            xw = xxmax - xxmin
            yw = yymax - yymin
            width = [xw, yw]
        elif name == "phi":
            if axes_transform is not AxesTransform.POLAR:
                raise NotImplementedError(
                    f"{axes_transform} is not implemented for normal axis {name!r}"
                )
            Rmin, Rmax, zmin, zmax = self._poloidal_bounds
            xw = Rmax - Rmin
            yw = zmax - zmin
            width = [xw, yw]
        return width
