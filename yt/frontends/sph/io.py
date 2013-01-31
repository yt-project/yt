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

_vector_fields = ("Coordinates", "Velocity", "Velocities")

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
                # This double-reads
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
            octree.add(pos, domain.domain_id)
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
        rv = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)
        psize = defaultdict(lambda: 0)
        chunks = list(chunks)
        ptypes = set()
        for ftype, fname in fields:
            ptf[ftype].append(fname)
            ptypes.add(ftype)
        ptypes = list(ptypes)
        ptypes.sort(key = lambda a: self._ptypes.index(a))
        for chunk in chunks:
            for subset in chunk.objs:
                poff = subset.domain.field_offsets
                tp = subset.domain.total_particles
                f = open(subset.domain.domain_filename, "rb")
                for ptype in ptypes:
                    f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                    pos = self._read_field_from_file(f,
                                tp[ptype], "Coordinates")
                    psize[ptype] += selector.count_points(
                        pos[:,0], pos[:,1], pos[:,2])
                    del pos
                f.close()
        ind = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", psize[field[0]], field)
            if field[1] in _vector_fields:
                shape = (psize[field[0]], 3)
            else:
                shape = psize[field[0]]
            rv[field] = np.empty(shape, dtype="float64")
            ind[field] = 0
        for chunk in chunks: 
            for subset in chunk.objs:
                poff = subset.domain.field_offsets
                tp = subset.domain.total_particles
                f = open(subset.domain.domain_filename, "rb")
                for ptype, field_list in sorted(ptf.items()):
                    f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                    pos = self._read_field_from_file(f,
                                tp[ptype], "Coordinates")
                    mask = selector.select_points(
                        pos[:,0], pos[:,1], pos[:,2])
                    del pos
                    if mask is None: continue
                    for field in field_list:
                        f.seek(poff[ptype, field], os.SEEK_SET)
                        data = self._read_field_from_file(f, tp[ptype], field)
                        data = data[mask]
                        my_ind = ind[ptype, field]
                        mylog.debug("Filling from %s to %s with %s",
                            my_ind, my_ind+data.shape[0], field)
                        rv[ptype, field][my_ind:my_ind + data.shape[0],...] = data
                        ind[ptype, field] += data.shape[0]
                f.close()
        return rv

    def _read_field_from_file(self, f, count, name):
        if count == 0: return
        if name == "ParticleIDs":
            dt = "int32"
        else:
            dt = "float32"
        if name in _vector_fields:
            count *= 3
        arr = np.fromfile(f, dtype=dt, count = count)
        if name in _vector_fields:
            arr = arr.reshape((count/3, 3), order="C")
        return arr.astype("float64")

    def _initialize_octree(self, domain, octree):
        count = sum(domain.total_particles.values())
        dt = [("px", "float32"), ("py", "float32"), ("pz", "float32")]
        with open(domain.domain_filename, "rb") as f:
            f.seek(self._header_offset)
            # The first total_particles * 3 values are positions
            pp = np.fromfile(f, dtype = dt, count = count)
        pos = np.empty((count, 3), dtype="float64")
        pos[:,0] = pp['px']
        pos[:,1] = pp['py']
        pos[:,2] = pp['pz']
        del pp
        octree.add(pos, domain.domain_id)

    def _count_particles(self, domain):
        npart = dict((self._ptypes[i], v)
            for i, v in enumerate(domain.header["Npart"])) 
        return npart

    _header_offset = 256
    _field_size = 4
    def _calculate_field_offsets(self, field_list, pcount):
        # field_list is (ftype, fname) but the blocks are ordered
        # (fname, ftype) in the file.
        pos = self._header_offset # 256 bytes for the header
        fs = self._field_size
        offsets = {}
        for field in self._fields:
            if not isinstance(field, types.StringTypes):
                field = field[0]
            for ptype in self._ptypes:
                if (ptype, field) not in field_list: continue
                offsets[(ptype, field)] = pos
                if field in _vector_fields:
                    pos += 3 * pcount[ptype] * fs
                else:
                    pos += pcount[ptype] * fs
        return offsets

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

class IOHandlerTipsyBinary(BaseIOHandler):
    _data_style = "tipsy"

    _pdtypes = None # dtypes, to be filled in later

    _ptypes = ( "Gas",
                "DarkMatter",
                "Stars" )

    _fields = ( ("Gas", "Mass"),
                ("Gas", "Coordinates"),
                ("Gas", "Velocities"),
                ("Gas", "Density"),
                ("Gas", "Temperature"),
                ("Gas", "Epsilon"),
                ("Gas", "Metals"),
                ("Gas", "Phi"),
                ("DarkMatter", "Mass"),
                ("DarkMatter", "Coordinates"),
                ("DarkMatter", "Velocities"),
                ("DarkMatter", "Epsilon"),
                ("DarkMatter", "Phi"),
                ("Stars", "Mass"),
                ("Stars", "Coordinates"),
                ("Stars", "Velocities"),
                ("Stars", "Metals"),
                ("Stars", "FormationTime"),
                ("Stars", "Epsilon"),
                ("Stars", "Phi")
              )

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _fill_fields(self, fields, vals, mask):
        if mask is None:
            size = 0
        else:
            size = mask.sum()
        rv = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", size, field)
            if field in _vector_fields:
                rv[field] = np.empty((size, 3), dtype="float64")
                if size == 0: continue
                rv[field][:,0] = vals[field]['x'][mask]
                rv[field][:,1] = vals[field]['y'][mask]
                rv[field][:,2] = vals[field]['z'][mask]
            else:
                rv[field] = np.empty(size, dtype="float64")
                if size == 0: continue
                rv[field][:] = vals[field]
        return rv

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)
        ptypes = set()
        for ftype, fname in fields:
            ptf[ftype].append(fname)
            ptypes.add(ftype)
        ptypes = list(ptypes)
        ptypes.sort(key = lambda a: self._ptypes.index(a))
        for chunk in chunks:
            for subset in chunk.objs:
                poff = subset.domain.field_offsets
                tp = subset.domain.total_particles
                f = open(subset.domain.domain_filename, "rb")
                for ptype in ptypes:
                    f.seek(poff[ptype], os.SEEK_SET)
                    p = np.fromfile(f, self._pdtypes[ptype], count=tp[ptype])
                    mask = selector.select_points(
                        p['Coordinates']['x'].astype("float64"),
                        p['Coordinates']['y'].astype("float64"),
                        p['Coordinates']['z'].astype("float64"))
                    tf = self._fill_fields(ptf[ptype], p, mask)
                    for field in tf:
                        rv[ptype, field] = tf[field]
                    del p, tf
                f.close()
        return rv

    def _initialize_octree(self, domain, octree):
        with open(domain.domain_filename, "rb") as f:
            f.seek(domain.pf._header_offset)
            for ptype in self._ptypes:
                # We'll just add the individual types separately
                count = domain.total_particles[ptype]
                if count == 0: continue
                pp = np.fromfile(f, dtype = self._pdtypes[ptype],
                                 count = count)
                pos = np.empty((count, 3), dtype="float64")
                mylog.info("Adding %0.3e %s particles", count, ptype)
                pos[:,0] = pp['Coordinates']['x']
                pos[:,1] = pp['Coordinates']['y']
                pos[:,2] = pp['Coordinates']['z']
                mylog.debug("Spanning: %0.3e .. %0.3e in x",
                            pos[:,0].min(), pos[:,0].max())
                mylog.debug("Spanning: %0.3e .. %0.3e in y",
                            pos[:,1].min(), pos[:,1].max())
                mylog.debug("Spanning: %0.3e .. %0.3e in z",
                            pos[:,2].min(), pos[:,2].max())
                del pp
                octree.add(pos, domain.domain_id)

    def _count_particles(self, domain):
        npart = {
            "Gas": domain.pf.parameters['nsph'],
            "Stars": domain.pf.parameters['nstar'],
            "DarkMatter": domain.pf.parameters['ndark']
        }
        return npart

    def _create_dtypes(self, domain):
        # We can just look at the particle counts.
        self._header_offset = domain.pf._header_offset
        self._pdtypes = {}
        pds = {}
        field_list = []
        tp = domain.total_particles
        for ptype, field in self._fields:
            pfields = []
            if tp[ptype] == 0: continue
            if field in _vector_fields:
                dt = (field, [('x', '>f'), ('y', '>f'), ('z', '>f')])
            else:
                dt = (field, '>f')
            pds.setdefault(ptype, []).append(dt)
            field_list.append((ptype, field))
        for ptype in pds:
            self._pdtypes[ptype] = np.dtype(pds[ptype])
        self._field_list = field_list
        return self._field_list

    def _identify_fields(self, domain):
        return self._field_list

    def _calculate_particle_offsets(self, domain):
        field_offsets = {}
        pos = domain.pf._header_offset
        for ptype in self._ptypes:
            field_offsets[ptype] = pos
            if domain.total_particles[ptype] == 0: continue
            size = self._pdtypes[ptype].itemsize
            pos += domain.total_particles[ptype] * size
        return field_offsets
