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

import numpy as np
from yt.utilities.exceptions import *
from yt.fields.field_info_container import \
    FieldInfoContainer

class FITSFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    def __init__(self, pf, field_list, slice_info=None):
        super(FITSFieldInfo, self).__init__(pf, field_list, slice_info=slice_info)
        for field in pf.field_list:
            if field[0] == "fits": self[field].take_log = False

    def _setup_ppv_fields(self):

        def _get_2d_wcs(data, axis):
            w_coords = data.pf.wcs_2d.wcs_pix2world(data["x"], data["y"], 1)
            return w_coords[axis]

        def world_f(axis, unit):
            def _world_f(field, data):
                return data.pf.arr(_get_2d_wcs(data, axis), unit)
            return _world_f

        for (i, axis), name in zip(enumerate([self.pf.lon_axis, self.pf.lat_axis]),
                                             [self.pf.lon_name, self.pf.lat_name]):
            unit = str(self.pf.wcs_2d.wcs.cunit[i])
            if unit.lower() == "deg": unit = "degree"
            if unit.lower() == "rad": unit = "radian"
            self.add_field(("fits",name), function=world_f(axis, unit), units=unit)

        if self.pf.dimensionality == 3:
            def _vel_los(field, data):
                axis = "xyz"[data.pf.vel_axis]
                return data.pf.arr(data[axis].ndarray_view(),data.pf.vel_unit)
            self.add_field(("fits",self.pf.vel_name), function=_vel_los,
                           units=self.pf.vel_unit)

    def setup_fluid_fields(self):

        if self.pf.ppv_data:
            def _pixel(field, data):
                return data.pf.arr(data["ones"], "pixel")
            self.add_field(("fits","pixel"), function=_pixel, units="pixel")
            self._setup_ppv_fields()
            return
