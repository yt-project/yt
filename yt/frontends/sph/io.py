"""
Gadget-specific data-file handling function

Author: Christopher E Moody <juxtaposicion@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Christopher E Moody, Matthew Turk.  All Rights Reserved.

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

import h5py
import numpy as np
from yt.funcs import *
from yt.utilities.exceptions import *

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.fortran_utils import read_record

_vector_fields = ("Coordinates", "Velocity")

class IOHandlerOWLS(BaseIOHandler):
    _data_style = "OWLS"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)
        psize = defaultdict(lambda: 0)
        chunks = list(chunks)
        for ftype, fname in fields:
            ptf[ftype].append(fname)
        for chunk in chunks: # Will be OWLS domains
            for subset in chunk.objs:
                f = h5py.File(subset.domain.domain_filename, "r")
                for ptype, field_list in sorted(ptf.items()):
                    coords = f["/%s/Coordinates" % ptype][:].astype("float64")
                    psize[ptype] += selector.count_points(
                        coords[:,0], coords[:,1], coords[:,2])
                    del coords
                f.close()
        # Now we have all the sizes, and we can allocate
        ind = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", psize[field[0]], field)
            if field[1] in _vector_fields:
                shape = (psize[field[0]], 3)
            else:
                shape = psize[field[0]]
            rv[field] = np.empty(shape, dtype="float64")
            ind[field] = 0
        for chunk in chunks: # Will be OWLS domains
            for subset in chunk.objs:
                f = h5py.File(subset.domain.domain_filename, "r")
                for ptype, field_list in sorted(ptf.items()):
                    g = f["/%s" % ptype]
                    coords = g["Coordinates"][:].astype("float64")
                    mask = selector.select_points(
                                coords[:,0], coords[:,1], coords[:,2])
                    del coords
                    if mask is None: continue
                    for field in field_list:
                        data = g[field][mask,...]
                        my_ind = ind[ptype, field]
                        mylog.debug("Filling from %s to %s with %s",
                            my_ind, my_ind+data.shape[0], field)
                        rv[ptype, field][my_ind:my_ind + data.shape[0],...] = data
                        ind[ptype, field] += data.shape[0]
                f.close()
        return rv

    def _initialize_octree(self, domain, octree):
        f = h5py.File(domain.domain_filename, "r")
        for key in f.keys():
            if not key.startswith("PartType"): continue
            pos = f[key]["Coordinates"][:].astype("float64")
            octree.add(pos, domain.domain_number)
        f.close()

    def _count_particles(self, domain):
        f = h5py.File(domain.domain_filename, "r")
        np = f["/Header"].attrs["NumPart_ThisFile"][:]
        f.close()
        npart = dict(("PartType%s" % (i), v) for i, v in enumerate(np)) 
        return npart

    def _identify_fields(self, domain):
        f = h5py.File(domain.domain_filename, "r")
        fields = []
        for key in f.keys():
            if not key.startswith("PartType"): continue
            g = f[key]
            #ptype = int(key[8:])
            ptype = str(key)
            for k in g.keys():
                if not hasattr(g[k], "shape"): continue
                # str => not unicode!
                fields.append((ptype, str(k)))
        f.close()
        return fields


ZeroMass = object()

class IOHandlerGadgetBinary(BaseIOHandler):
    _data_style = "gadget_binary"

    # Particle types (Table 3 in GADGET-2 user guide)
    _ptypes = ( "Gas",
                "Halo",
                "Disk",
                "Bulge",
                "Stars",
                "Bndry" )
    #
    # Blocks in the file:
    #   HEAD
    #   POS
    #   VEL
    #   ID
    #   MASS    (variable mass only)
    #   U       (gas only)
    #   RHO     (gas only)
    #   HSML    (gas only)
    #   POT     (only if enabled in makefile)
    #   ACCE    (only if enabled in makefile)
    #   ENDT    (only if enabled in makefile)
    #   TSTP    (only if enabled in makefile)

    _fields = ( "Coordinates",
                "Velocities",
                "ParticleIDs",
                ("Masses", ZeroMass),
                ("InternalEnergy", "Gas"),
                ("Density", "Gas"),
                ("SmoothingLength", "Gas"),
    )

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_selection(self, chunks, selector, fields):
        raise NotImplementedError

    def _initialize_octree(self, domain, octree):
        count = sum(domain.total_particles.values())
        dt = [("px", "float32"), ("py", "float32"), ("pz", "float32")]
        with open(domain.domain_filename, "rb") as f:
            # The first total_particles * 3 values are positions
            pp = np.fromfile(f, dtype = dt, count = count)
        pos = np.empty((count, 3), dtype="float64")
        pos[:,0] = pp['px']
        pos[:,1] = pp['py']
        pos[:,2] = pp['pz']
        del pp
        octree.add(pos, domain.domain_number)

    def _count_particles(self, domain):
        npart = dict((self._ptypes[i], v)
            for i, v in enumerate(domain.header["Npart"])) 
        return npart

    def _identify_fields(self, domain):
        # We can just look at the particle counts.
        field_list = []
        tp = domain.total_particles
        for i, ptype in enumerate(self._ptypes):
            count = tp[ptype]
            if count == 0: continue
            m = domain.header["Massarr"][i]
            for field in self._fields:
                if isinstance(field, types.TupleType):
                    field, req = field
                    if req is ZeroMass:
                        if m > 0.0 : continue
                    elif req != field:
                        continue
                field_list.append((ptype, field))
        return field_list
