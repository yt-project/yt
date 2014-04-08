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
    def _get_wcs(self, data, axis):
        if data.pf.dimensionality == 2:
            w_coords = data.pf.wcs.wcs_pix2world(data["x"], data["y"], 1)
        else:
            w_coords = data.pf.wcs.wcs_pix2world(data["x"], data["y"],
                                                 data["z"], 1)
        return w_coords[axis]
    def setup_fluid_fields(self):
        def world_f(axis):
            def _world_f(field, data):
                return self._get_wcs(data, axis)
            return _world_f
        for i in range(self.pf.dimensionality):
            if self.pf.wcs.wcs.cname[i] == '':
                name = str(self.pf.wcs.wcs.ctype[i])
            else:
                name = str(self.pf.wcs.wcs.cname[i])
            unit = str(self.pf.wcs.wcs.cunit[i])
            if name != '' and unit != '':
                if unit.lower() == "deg": unit = "degree"
                if unit.lower() == "rad": unit = "radian"
                self.add_field(("fits",name), function=world_f(i), units=unit)



