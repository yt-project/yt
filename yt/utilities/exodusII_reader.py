import string
from itertools import takewhile
from netCDF4 import Dataset
import numpy as np
from yt.config import ytcfg
import os


def sanitize_string(s):
    s = "".join(_ for _ in takewhile(lambda a: a in string.printable, s))
    return s


def get_data(fn):
    try:
        f = Dataset(fn)
    except RuntimeError:
        f = Dataset(os.path.join(ytcfg.get("yt", "test_data_dir"), fn))
    fvars = f.variables
    # Is this correct?
    etypes = fvars["eb_status"][:]
    nelem = etypes.shape[0]
    varnames = [sanitize_string(v.tostring()) for v in
                fvars["name_elem_var"][:]]
    nodnames = [sanitize_string(v.tostring()) for v in
                fvars["name_nod_var"][:]]
    coord = np.array([fvars["coord%s" % ax][:]
                     for ax in 'xyz']).transpose().copy()
    coords = []
    connects = []
    data = []
    for i in range(nelem):
        connects.append(fvars["connect%s" % (i+1)][:].astype("i8"))
        ci = connects[-1]
        coords.append(coord)  # Same for all
        vals = {}
        for j, v in enumerate(varnames):
            values = fvars["vals_elem_var%seb%s" % (j+1, i+1)][:]
            vals['gas', v] = values.astype("f8")[-1, :]
        for j, v in enumerate(nodnames):
            # We want just for this set of nodes all the node variables
            # Use (ci - 1) to get these values
            values = fvars["vals_nod_var%s" % (j+1)][:]
            vals['gas', v] = values.astype("f8")[-1, ci - 1, ...]
        data.append(vals)
    return coords, connects, data
