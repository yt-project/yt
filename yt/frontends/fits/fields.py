"""
FITS-specific fields
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

class FITSFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    def __init__(self, ds, field_list, slice_info=None):
        super(FITSFieldInfo, self).__init__(ds, field_list, slice_info=slice_info)
        for field in ds.field_list:
            if field[0] == "fits": self[field].take_log = False

    def _setup_spec_cube_fields(self):

        def _get_2d_wcs(data, axis):
            w_coords = data.ds.wcs_2d.wcs_pix2world(data["x"], data["y"], 1)
            return w_coords[axis]

        def world_f(axis, unit):
            def _world_f(field, data):
                return data.ds.arr(_get_2d_wcs(data, axis), unit)
            return _world_f

        for (i, axis), name in zip(enumerate([self.ds.lon_axis, self.ds.lat_axis]),
                                             [self.ds.lon_name, self.ds.lat_name]):
            unit = str(self.ds.wcs_2d.wcs.cunit[i])
            if unit.lower() == "deg": unit = "degree"
            if unit.lower() == "rad": unit = "radian"
            self.add_field(("fits",name), function=world_f(axis, unit), units=unit)

        if self.ds.dimensionality == 3:
            def _spec(field, data):
                axis = "xyz"[data.ds.spec_axis]
                sp = (data[axis].ndarray_view()-self.ds._p0)*self.ds._dz + self.ds._z0
                return data.ds.arr(sp, data.ds.spec_unit)
            self.add_field(("fits","spectral"), function=_spec,
                           units=self.ds.spec_unit, display_name=self.ds.spec_name)

    def setup_fluid_fields(self):

        if self.ds.spec_cube:
            def _pixel(field, data):
                return data.ds.arr(data["ones"], "pixel")
            self.add_field(("fits","pixel"), function=_pixel, units="pixel")
            self._setup_spec_cube_fields()
            return
