from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer


class FITSFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    def __init__(self, ds, field_list, slice_info=None):
        super().__init__(ds, field_list, slice_info=slice_info)
        for field in ds.field_list:
            if field[0] == "fits":
                self[field].take_log = False


class YTFITSFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("density", ("code_mass/code_length**3", ["density"], None)),
        (
            "dark_matter_density",
            ("code_mass/code_length**3", ["dark_matter_density"], None),
        ),
        ("number_density", ("1/code_length**3", ["number_density"], None)),
        ("pressure", ("dyne/code_length**2", ["pressure"], None)),
        ("thermal_energy", ("erg / g", ["specific_thermal_energy"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("velocity_x", ("code_length/code_time", ["velocity_x"], None)),
        ("velocity_y", ("code_length/code_time", ["velocity_y"], None)),
        ("velocity_z", ("code_length/code_time", ["velocity_z"], None)),
        ("magnetic_field_x", ("gauss", [], None)),
        ("magnetic_field_y", ("gauss", [], None)),
        ("magnetic_field_z", ("gauss", [], None)),
        ("metallicity", ("Zsun", ["metallicity"], None)),
        # We need to have a bunch of species fields here, too
        ("metal_density", ("code_mass/code_length**3", ["metal_density"], None)),
        ("hi_density", ("code_mass/code_length**3", ["hi_density"], None)),
        ("hii_density", ("code_mass/code_length**3", ["hii_density"], None)),
        ("h2i_density", ("code_mass/code_length**3", ["h2i_density"], None)),
        ("h2ii_density", ("code_mass/code_length**3", ["h2ii_density"], None)),
        ("h2m_density", ("code_mass/code_length**3", ["h2m_density"], None)),
        ("hei_density", ("code_mass/code_length**3", ["hei_density"], None)),
        ("heii_density", ("code_mass/code_length**3", ["heii_density"], None)),
        ("heiii_density", ("code_mass/code_length**3", ["heiii_density"], None)),
        ("hdi_density", ("code_mass/code_length**3", ["hdi_density"], None)),
        ("di_density", ("code_mass/code_length**3", ["di_density"], None)),
        ("dii_density", ("code_mass/code_length**3", ["dii_density"], None)),
    )

    def __init__(self, ds, field_list, slice_info=None):
        super().__init__(ds, field_list, slice_info=slice_info)


class WCSFITSFieldInfo(FITSFieldInfo):
    def setup_fluid_fields(self):
        wcs_2d = getattr(self.ds, "wcs_2d", self.ds.wcs)

        def _pixel(field, data):
            return data.ds.arr(data[("index", "ones")], "pixel")

        self.add_field(
            ("fits", "pixel"), sampling_type="cell", function=_pixel, units="pixel"
        )

        def _get_2d_wcs(data, axis):
            w_coords = wcs_2d.wcs_pix2world(
                data[("index", "x")], data[("index", "y")], 1
            )
            return w_coords[axis]

        def world_f(axis, unit):
            def _world_f(field, data):
                return data.ds.arr(_get_2d_wcs(data, axis), unit)

            return _world_f

        for (i, axis), name in zip(
            enumerate([self.ds.lon_axis, self.ds.lat_axis]),
            [self.ds.lon_name, self.ds.lat_name],
        ):
            unit = str(wcs_2d.wcs.cunit[i])
            if unit.lower() == "deg":
                unit = "degree"
            if unit.lower() == "rad":
                unit = "radian"
            self.add_field(
                ("fits", name),
                sampling_type="cell",
                function=world_f(axis, unit),
                units=unit,
            )

        if self.ds.dimensionality == 3:

            def _spec(field, data):
                axis = "xyz"[data.ds.spec_axis]
                sp = (
                    data[("fits", axis)].ndarray_view() - self.ds._p0
                ) * self.ds._dz + self.ds._z0
                return data.ds.arr(sp, data.ds.spec_unit)

            self.add_field(
                ("fits", "spectral"),
                sampling_type="cell",
                function=_spec,
                units=self.ds.spec_unit,
                display_name=self.ds.spec_name,
            )
