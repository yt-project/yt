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
from yt.utilities.on_demand_imports import ap

def _make_counts(emin, emax):
    def _counts(field, data):
        e = data["event_energy"].in_units("keV")
        mask = np.logical_and(e >= emin, e < emax)
        x = data["event_x"][mask]
        y = data["event_y"][mask]
        z = np.ones(x.shape)
        pos = np.array([x,y,z]).transpose()
        img = data.deposit(pos, method="count")
        if data.has_field_parameter("sigma"):
            sigma = data.get_field_parameter("sigma")
        else:
            sigma = None
        if sigma is not None and sigma > 0.0:
            kern = ap.conv.Gaussian2DKernel(stddev=sigma)
            img[:,:,0] = ap.conv.convolve(img[:,:,0], kern)
        return data.pf.arr(img, "counts/pixel")
    return _counts

def setup_counts_fields(ds, ebounds):
    r"""
    Create deposited image fields from X-ray count data in energy bands.

    Parameters
    ----------
    ds : Dataset
        The FITS events file dataset to add the counts fields to.
    ebounds : list of tuples
        A list of tuples, one for each field, with (emin, emax) as the
        energy bounds for the image.

    Examples
    --------
    >>> ds = yt.load("evt.fits")
    >>> ebounds = [(0.1,2.0),(2.0,3.0)]
    >>> setup_counts_fields(ds, ebounds)
    """
    for (emin, emax) in ebounds:
        cfunc = _make_counts(emin, emax)
        fname = "counts_%s-%s" % (emin, emax)
        mylog.info("Creating counts field %s." % fname)
        ds.add_field(("gas",fname), function=cfunc,
                     units="counts/pixel",
                     validators = [ValidateSpatial()],
                     display_name="Counts (%s-%s keV)" % (emin, emax))