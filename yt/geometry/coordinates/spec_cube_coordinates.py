from typing import Dict

from yt._maintenance.deprecation import issue_deprecation_warning

from ._axes_transforms import AxesTransform
from .cartesian_coordinates import (
    CartesianCoordinateHandler,
    DefaultProperties,
)
from .coordinate_handler import _get_coord_fields


class SpectralCubeCoordinateHandler(CartesianCoordinateHandler):
    name = "spectral_cube"
    _default_axes_transforms: Dict[str, AxesTransform] = {
        "x": AxesTransform.GEOMETRY_NATIVE,
        "y": AxesTransform.GEOMETRY_NATIVE,
        "z": AxesTransform.GEOMETRY_NATIVE,
    }

    def __init__(self, ds, ordering=None):
        if ordering is None:
            ordering = tuple(
                "xyz"[axis] for axis in (ds.lon_axis, ds.lat_axis, ds.spec_axis)
            )
        super().__init__(ds, ordering)

        def _spec_axis(ax, x, y):
            p = (x, y)[ax]
            return [self.ds.pixel2spec(pp).v for pp in p]

        self.axis_field = {}
        self.axis_field[self.ds.spec_axis] = _spec_axis

    def setup_fields(self, registry):
        if not self.ds.no_cgs_equiv_length:
            return super().setup_fields(registry)
        for axi, ax in enumerate("xyz"):
            f1, f2 = _get_coord_fields(axi)

            def _get_length_func():
                def _length_func(field, data):
                    # Just use axis 0
                    rv = data.ds.arr(data.fcoords[..., 0].copy(), field.units)
                    rv[:] = 1.0
                    return rv

                return _length_func

            registry.add_field(
                ("index", f"d{ax}"),
                sampling_type="cell",
                function=f1,
                display_field=False,
                units="code_length",
            )

            registry.add_field(
                ("index", f"path_element_{ax}"),
                sampling_type="cell",
                function=_get_length_func(),
                display_field=False,
                units="",
            )

            registry.add_field(
                ("index", f"{ax}"),
                sampling_type="cell",
                function=f2,
                display_field=False,
                units="code_length",
            )

        self._register_volume(registry)
        self._check_fields(registry)

    _x_pairs = (("x", "y"), ("y", "x"), ("z", "x"))
    _y_pairs = (("x", "z"), ("y", "z"), ("z", "y"))

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    @property
    def image_axis_name(self):
        issue_deprecation_warning(
            "The image_axis_name property isn't used "
            "internally in yt anymore and is deprecated",
            since="4.2.0",
        )
        return self._image_axis_name

    def _get_plot_axes_default_properties(
        self, normal_axis_name: str, axes_transform: AxesTransform
    ) -> DefaultProperties:
        if axes_transform is AxesTransform.DEFAULT:
            axes_transform = AxesTransform.GEOMETRY_NATIVE
        if axes_transform is not AxesTransform.GEOMETRY_NATIVE:
            raise NotImplementedError(
                f"spectral cube coordinates don't implement {axes_transform} yet"
            )

        if normal_axis_name == "x":
            if self.ds.lon_name == "X" and self.ds.lat_name == "Y":
                return {
                    "x_axis_label": "x",
                    "y_axis_label": self.ds.spec_name,
                    "x_axis_units": None,
                    "y_axis_units": None,
                }
            else:
                return {
                    "x_axis_label": r"Image\ x",
                    "y_axis_label": self.ds.spec_name,
                    "x_axis_units": "pixel",
                    "y_axis_units": None,
                }
        elif normal_axis_name == "y":
            {
                "x_axis_label": "x",
                "y_axis_label": self.ds.spec_name,
                "x_axis_units": None,
                "y_axis_units": None,
            }
            if self.ds.lon_name == "X" and self.ds.lat_name == "Y":
                return {
                    "x_axis_label": "x",
                    "y_axis_label": self.ds.spec_name,
                    "x_axis_units": None,
                    "y_axis_units": None,
                }
            else:
                return {
                    "x_axis_label": r"Image\ x",
                    "y_axis_label": self.ds.spec_name,
                    "x_axis_units": "pixel",
                    "y_axis_units": None,
                }
        elif normal_axis_name == "z":
            if self.ds.lon_name == "X" and self.ds.lat_name == "Y":
                return {
                    "x_axis_label": "x",
                    "y_axis_label": "y",
                    "x_axis_units": None,
                    "y_axis_units": None,
                }
            else:
                return {
                    "x_axis_label": r"Image\ x",
                    "y_axis_label": r"Image\ y",
                    "x_axis_units": "pixel",
                    "y_axis_units": "pixel",
                }
        else:
            raise ValueError(f"Unknown axis {normal_axis_name!r}")
