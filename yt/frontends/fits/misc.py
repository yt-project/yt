"""
Miscellaneous FITS routines
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.fields.api import add_field
from yt.fields.derived_field import ValidateSpatial
from yt.funcs import mylog
import os

def _make_counts(emin, emax):
    def _counts(field, data):
        e = data["event_energy"].in_units("keV")
        mask = np.logical_and(e >= emin, e < emax)
        x = data["event_x"][mask]
        y = data["event_y"][mask]
        z = np.ones(x.shape)
        pos = np.array([x,y,z]).transpose()
        img = data.deposit(pos, method="count")
        return data.pf.arr(img, "counts/pixel")
    return _counts

def setup_counts_fields(ebounds):
    for (emin, emax) in ebounds:
        cfunc = _make_counts(emin, emax)
        fname = "counts_%s-%s" % (emin, emax)
        mylog.info("Creating counts field %s." % fname)
        add_field(("gas",fname), function=cfunc,
                  units="counts/pixel",
                  validators = [ValidateSpatial()],
                  display_name="Counts (%s-%s keV)" % (emin, emax))

def ds9_region(ds, reg, obj=None):
    import pyregion
    r = pyregion.open(reg)
    reg_name = reg.split(".")[0]
    nx = ds.domain_dimensions[ds.lon_axis]
    ny = ds.domain_dimensions[ds.lat_axis]
    mask = r.get_mask(header=ds.wcs_2d.to_header(),
                      shape=[nx,ny])
    def _reg_field(field, data):
        i = data["x"].ndarray_view().astype("int")-1
        j = data["y"].ndarray_view().astype("int")-1
        new_mask = mask[i,j]
        ret = data["zeros"].copy()
        ret[new_mask] = 1.
        return ret
    ds.index
    ds.field_info.add_field(("gas",reg_name), function=_reg_field)
    if obj is None:
        obj = ds.all_data()
    return obj.cut_region(["obj['%s'] > 0" % (reg_name)])