"""
Some convenience functions, objects, and iterators

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import glob
import numpy as na
import os, os.path, inspect, types
from functools import wraps

# Named imports
from yt.funcs import *
from yt.config import ytcfg
from yt.utilities.parameter_file_storage import \
    output_type_registry, \
    EnzoRunDatabase

def all_pfs(basedir='.', skip=None, max_depth=1, name_spec="*.hierarchy", **kwargs):
    """
    This function searchs a directory and its sub-directories, up to a
    depth of *max_depth*, for parameter files.  It looks for the
    *name_spec* and then instantiates an EnzoStaticOutput from
    each. You can skip every *skip* parameter files, if *skip* is not
    None; otherwise it will return all files.  All subsequent *kwargs*
    are passed on to the EnzoStaticOutput constructor.
    """
    list_of_names = []
    basedir = os.path.expanduser(basedir)
    for i in range(max_depth):
        bb = list('*' * i) + [name_spec]
        list_of_names += glob.glob(os.path.join(basedir,*bb))
    list_of_names.sort(key=lambda b: os.path.basename(b))
    for fn in list_of_names[::skip]:
        yield load(fn[:-10], **kwargs)

def max_spheres(width, unit, **kwargs):
    """
    This calls :func:`~yt.convenience.all_pfs` and then for each parameter file
    creates a :class:`~yt.data_objects.api.AMRSphereBase` for each one,
    centered on the point of highest density, with radius *width* in units of
    *unit*.
    """
    for pf in all_pfs(**kwargs):
        v, c = pf.h.find_max("Density")
        yield pf.h.sphere(c, width/pf[unit])

def load(*args ,**kwargs):
    """
    This function attempts to determine the base data type of a filename or
    other set of arguments by calling
    :meth:`yt.data_objects.api.StaticOutput._is_valid` until it finds a
    match, at which point it returns an instance of the appropriate
    :class:`yt.data_objects.api.StaticOutput` subclass.
    """
    if len(args) == 0:
        try:
            import Tkinter, tkFileDialog
        except ImportError:
            return None
        root = Tkinter.Tk()
        filename = tkFileDialog.askopenfilename(parent=root,title='Choose a file')
        if filename != None:
            return load(filename)
        else:
            return None
    candidates = []
    args = [os.path.expanduser(arg) if isinstance(arg, types.StringTypes)
            else arg for arg in args]
    valid_file = [os.path.isfile(arg) if isinstance(arg, types.StringTypes) 
            else False for arg in args]
    if not any(valid_file):
        mylog.error("None of the arguments provided to load() is a valid file")
        mylog.error("Please check that you have used a correct path")
        return None
    for n, c in output_type_registry.items():
        if n is None: continue
        if c._is_valid(*args, **kwargs): candidates.append(n)
    if len(candidates) == 1:
        return output_type_registry[candidates[0]](*args, **kwargs)
    if len(candidates) == 0:
        if ytcfg.get("yt", "enzo_db") != '' \
           and len(args) == 1 \
           and isinstance(args[0], types.StringTypes):
            erdb = EnzoRunDatabase()
            fn = erdb.find_uuid(args[0])
            n = "EnzoStaticOutput"
            if n in output_type_registry \
               and output_type_registry[n]._is_valid(fn):
                return output_type_registry[n](fn)
        mylog.error("Couldn't figure out output type for %s", args[0])
        return None
    mylog.error("Multiple output type candidates for %s:", args[0])
    for c in candidates:
        mylog.error("    Possible: %s", c)
    return None

def projload(pf, axis, weight_field = None):
    # This is something of a hack, so that we can just get back a projection
    # and not utilize any of the intermediate hierarchy objects.
    class ProjMock(dict):
        pass
    import h5py
    f = h5py.File(os.path.join(pf.fullpath, pf.parameter_filename + ".yt"))
    b = f["/Projections/%s/" % (axis)]
    wf = "weight_field_%s" % weight_field
    if wf not in b: raise KeyError(wf)
    fields = []
    for k in b:
        if k.startswith("weight_field"): continue
        if k.endswith("_%s" % weight_field):
            fields.append(k)
    proj = ProjMock()
    for f in ["px","py","pdx","pdy"]:
        proj[f] = b[f][:]
    for f in fields:
        new_name = f[:-(len(weight_field) + 1)]
        proj[new_name] = b[f][:]
    proj.axis = axis
    proj.pf = pf
    f.close()
    return proj

def _chunk(arrlike, chunksize = 800000):
    total_size = arrlike.shape[0]
    pbar = get_pbar("Transferring %s " % (arrlike.name), total_size)
    start = 0; end = 0
    bits = []
    while start < total_size:
        bits.append(arrlike[start:start+chunksize])
        pbar.update(start)
        start += chunksize
    pbar.finish()
    return na.concatenate(bits)

def dapload(p, axis, weight_field = None):
    r"""Load a projection dataset from a DAP server.

    If you have projections stored externally on a DAP server, this function
    can load them (transferring in chunks to avoid overloading) locally and
    display them.

    Parameters
    ----------
    p : string
        URL for the dataset on the DAP server
    axis : int
        The axis of projection to load (0, 1, 2)
    weight_field : string
        The weight_field used in the projection

    Returns
    -------
    projmock : ProjMock
        This is a mockup of a projection that mostly fills the API.  It can be
        used with `yt.visualization.image_panner.api.VariableMeshPanner`
        objects.

    See Also
    --------
    http://www.opendap.org/ and http://pydap.org/2.x/ . (Note that HDF5 is not
    supported on PyDAP 3.x servers.)

    Examples
    --------

    >>> p = "http://datasets-r-us.org/output_0013.h5"
    >>> proj = dapload(p, 0, "Density")
    >>> vmp = VariableMeshPanner(proj, (512, 512), "Density", ImageSaver(0))
    >>> vmp.zoom(1.0)
    """
    class PFMock(dict):
        domain_left_edge = na.zeros(3, dtype='float64')
        domain_right_edge = na.ones(3, dtype='float64')
    pf = PFMock()
    class ProjMock(dict):
        pass
    import dap.client
    f = dap.client.open(p)
    b = f["Projections"]["%s" % (axis)]
    wf = "weight_field_%s" % weight_field
    if wf not in b: raise KeyError(wf)
    fields = []
    for k in b:
        if k.name.startswith("weight_field"): continue
        if k.name.endswith("_%s" % weight_field):
            fields.append(k.name)
    proj = ProjMock()
    for f in ["px","py","pdx","pdy"]:
        proj[f] = _chunk(b[f])
    for f in fields:
        new_name = f[:-(len(str(weight_field)) + 1)]
        proj[new_name] = _chunk(b[f])
    proj.axis = axis
    proj.pf = pf
    return proj

