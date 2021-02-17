from .cartesian_coordinates import CartesianCoordinateHandler
from .coordinate_handler import _get_coord_fields


class SpectralCubeCoordinateHandler(CartesianCoordinateHandler):
    name = "spectral_cube"

    def __init__(self, ds, ordering=("x", "y", "z")):
        ordering = tuple(
            "xyz"[axis] for axis in (ds.lon_axis, ds.lat_axis, ds.spec_axis)
        )
        super().__init__(ds, ordering)

        self.default_unit_label = {}
        names = {}
        if ds.lon_name != "X" or ds.lat_name != "Y":
            names["x"] = r"Image\ x"
            names["y"] = r"Image\ y"
            # We can just use ds.lon_axis here
            self.default_unit_label[ds.lon_axis] = "pixel"
            self.default_unit_label[ds.lat_axis] = "pixel"
        names["z"] = ds.spec_name
        # Again, can use spec_axis here
        self.default_unit_label[ds.spec_axis] = ds.spec_unit

        self._image_axis_name = ian = {}
        for ax in "xyz":
            axi = self.axis_id[ax]
            xax = self.axis_name[self.x_axis[ax]]
            yax = self.axis_name[self.y_axis[ax]]
            ian[axi] = ian[ax] = ian[ax.upper()] = (
                names.get(xax, xax),
                names.get(yax, yax),
            )

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

        def _cell_volume(field, data):
            rv = data["index", "dx"].copy(order="K")
            rv *= data["index", "dy"]
            rv *= data["index", "dz"]
            return rv

        registry.add_field(
            ("index", "cell_volume"),
            sampling_type="cell",
            function=_cell_volume,
            display_field=False,
            units="code_length**3",
        )
        registry.alias(("index", "volume"), ("index", "cell_volume"))

        registry.check_derived_fields(
            [
                ("index", "dx"),
                ("index", "dy"),
                ("index", "dz"),
                ("index", "x"),
                ("index", "y"),
                ("index", "z"),
                ("index", "cell_volume"),
            ]
        )

    _x_pairs = (("x", "y"), ("y", "x"), ("z", "x"))
    _y_pairs = (("x", "z"), ("y", "z"), ("z", "y"))

    def convert_to_cylindrical(self, coord):
        raise NotImplementedError

    def convert_from_cylindrical(self, coord):
        raise NotImplementedError

    @property
    def image_axis_name(self):
        return self._image_axis_name
