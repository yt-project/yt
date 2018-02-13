"""
Gadget data-file handling functions




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os

from yt.extern.six import string_types
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .data_structures import \
    _get_gadget_format

from .definitions import \
    gadget_hdf5_ptypes, \
    SNAP_FORMAT_2_OFFSET


class IOHandlerGadgetHDF5(BaseIOHandler):
    _dataset_type = "gadget_hdf5"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")
    _known_ptypes = gadget_hdf5_ptypes
    _var_mass = None
    _element_names = ('Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                      'Neon', 'Magnesium', 'Silicon', 'Iron')

    @property
    def var_mass(self):
        if self._var_mass is None:
            vm = []
            for i, v in enumerate(self.ds["Massarr"]):
                if v == 0:
                    vm.append(self._known_ptypes[i])
            self._var_mass = tuple(vm)
        return self._var_mass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: x.filename):
            f = h5py.File(data_file.filename, "r")
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                x = f["/%s/Coordinates" % ptype][:, 0].astype("float64")
                y = f["/%s/Coordinates" % ptype][:, 1].astype("float64")
                z = f["/%s/Coordinates" % ptype][:, 2].astype("float64")
                yield ptype, (x, y, z)
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: x.filename):
            f = h5py.File(data_file.filename, "r")
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                g = f["/%s" % ptype]
                coords = g["Coordinates"][:].astype("float64")
                mask = selector.select_points(
                    coords[:, 0], coords[:, 1], coords[:, 2], 0.0)
                del coords
                if mask is None:
                    continue
                for field in field_list:

                    if field in ("Mass", "Masses") and \
                            ptype not in self.var_mass:
                        data = np.empty(mask.sum(), dtype="float64")
                        ind = self._known_ptypes.index(ptype)
                        data[:] = self.ds["Massarr"][ind]

                    elif field in self._element_names:
                        rfield = 'ElementAbundance/' + field
                        data = g[rfield][:][mask, ...]
                    elif field.startswith("Metallicity_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["Metallicity"][:, col][mask]
                    elif field.startswith("Chemistry_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["ChemistryAbundances"][:, col][mask]
                    else:
                        data = g[field][:][mask, ...]

                    yield (ptype, field), data
            f.close()

    def _initialize_index(self, data_file, regions):
        index_ptype = self.index_ptype
        f = h5py.File(data_file.filename, "r")
        if index_ptype == "all":
            pcount = f["/Header"].attrs["NumPart_ThisFile"][:].sum()
            keys = f.keys()
        else:
            pt = int(index_ptype[-1])
            pcount = f["/Header"].attrs["NumPart_ThisFile"][pt]
            keys = [index_ptype]
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        for key in keys:
            if not key.startswith("PartType"):
                continue
            if "Coordinates" not in f[key]:
                continue
            ds = f[key]["Coordinates"]
            dt = ds.dtype.newbyteorder("N")  # Native
            pos = np.empty(ds.shape, dtype=dt)
            pos[:] = ds
            regions.add_data_file(pos, data_file.file_id,
                                  data_file.ds.filter_bbox)
            morton[ind:ind + pos.shape[0]] = compute_morton(
                pos[:, 0], pos[:, 1], pos[:, 2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge,
                data_file.ds.filter_bbox)
            ind += pos.shape[0]
        f.close()
        return morton

    def _count_particles(self, data_file):
        f = h5py.File(data_file.filename, "r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:]
        f.close()
        npart = dict(("PartType%s" % (i), v) for i, v in enumerate(pcount))
        return npart

    def _identify_fields(self, data_file):
        f = h5py.File(data_file.filename, "r")
        fields = []
        cname = self.ds._particle_coordinates_name  # Coordinates
        mname = self.ds._particle_mass_name  # Mass

        # loop over all keys in OWLS hdf5 file
        #--------------------------------------------------
        for key in f.keys():

            # only want particle data
            #--------------------------------------
            if not key.startswith("PartType"):
                continue

            # particle data group
            #--------------------------------------
            g = f[key]
            if cname not in g:
                continue

            # note str => not unicode!
            ptype = str(key)
            if ptype not in self.var_mass:
                fields.append((ptype, mname))

            # loop over all keys in PartTypeX group
            #----------------------------------------
            for k in g.keys():

                if k == 'ElementAbundance':
                    gp = g[k]
                    for j in gp.keys():
                        kk = j
                        fields.append((ptype, str(kk)))
                elif k == 'Metallicity' and len(g[k].shape) > 1:
                    # Vector of metallicity
                    for i in range(g[k].shape[1]):
                        fields.append((ptype, "Metallicity_%02i" % i))
                elif k == "ChemistryAbundances" and len(g[k].shape) > 1:
                    for i in range(g[k].shape[1]):
                        fields.append((ptype, "Chemistry_%03i" % i))
                else:
                    kk = k
                    if not hasattr(g[kk], "shape"):
                        continue
                    if len(g[kk].shape) > 1:
                        self._vector_fields[kk] = g[kk].shape[1]
                    fields.append((ptype, str(kk)))

        f.close()
        return fields, {}


ZeroMass = object()


class IOHandlerGadgetBinary(BaseIOHandler):
    _dataset_type = "gadget_binary"
    _vector_fields = (("Coordinates", 3),
                      ("Velocity", 3),
                      ("Velocities", 3),
                      ("FourMetalFractions", 4))

    # Particle types (Table 3 in GADGET-2 user guide)
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

    _var_mass = None
    _format = None

    def __init__(self, ds, *args, **kwargs):
        self._vector_fields = dict(self._vector_fields)
        self._fields = ds._field_spec
        self._ptypes = ds._ptype_spec
        self.data_files = set([])
        gformat = _get_gadget_format(ds.parameter_filename)
        # gadget format 1 original, 2 with block name
        self._format = gformat[0]
        self._endian = gformat[1]
        super(IOHandlerGadgetBinary, self).__init__(ds, *args, **kwargs)

    @property
    def var_mass(self):
        if self._var_mass is None:
            vm = []
            for i, v in enumerate(self.ds["Massarr"]):
                if v == 0:
                    vm.append(self._ptypes[i])
            self._var_mass = tuple(vm)
        return self._var_mass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype in ptf:
                # This is where we could implement sub-chunking
                f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                pos = self._read_field_from_file(f,
                                                 tp[ptype], "Coordinates")
                yield ptype, (pos[:, 0], pos[:, 1], pos[:, 2])
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items()):
                f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                pos = self._read_field_from_file(f,
                                                 tp[ptype], "Coordinates")
                mask = selector.select_points(
                    pos[:, 0], pos[:, 1], pos[:, 2], 0.0)
                del pos
                if mask is None:
                    continue
                for field in field_list:
                    if field == "Mass" and ptype not in self.var_mass:
                        data = np.empty(mask.sum(), dtype="float64")
                        m = self.ds.parameters["Massarr"][
                            self._ptypes.index(ptype)]
                        data[:] = m
                        yield (ptype, field), data
                        continue
                    f.seek(poff[ptype, field], os.SEEK_SET)
                    data = self._read_field_from_file(f, tp[ptype], field)
                    data = data[mask, ...]
                    yield (ptype, field), data
            f.close()

    def _read_field_from_file(self, f, count, name):
        if count == 0:
            return
        if name == "ParticleIDs":
            dt = self._endian + "u4"
        else:
            dt = self._endian + self._float_type
        if name in self._vector_fields:
            count *= self._vector_fields[name]
        arr = np.fromfile(f, dtype=dt, count=count)
        if name in self._vector_fields:
            factor = self._vector_fields[name]
            arr = arr.reshape((count // factor, factor), order="C")
        return arr

    def _get_morton_from_position(self, data_file, count, offset_count,
                                  regions, DLE, DRE):
        with open(data_file.filename, "rb") as f:
            # We add on an additionally 4 for the first record.
            f.seek(data_file._position_offset + 4 + offset_count * 12)
            # The first total_particles * 3 values are positions
            pp = np.fromfile(f, dtype=self._endian + self._float_type,
                             count=count * 3)
            pp.shape = (count, 3)
            pp = pp.astype(self._float_type)
        regions.add_data_file(pp, data_file.file_id,
                                  data_file.ds.filter_bbox)
        morton = compute_morton(pp[:, 0], pp[:, 1], pp[:, 2], DLE, DRE,
                                data_file.ds.filter_bbox)
        return morton

    def _initialize_index(self, data_file, regions):
        DLE = data_file.ds.domain_left_edge
        DRE = data_file.ds.domain_right_edge
        self._float_type = data_file.ds._validate_header(data_file.filename)[1]
        if self.index_ptype == "all":
            count = sum(data_file.total_particles.values())
            return self._get_morton_from_position(
                data_file, count, 0, regions, DLE, DRE)
        else:
            idpos = self._ptypes.index(self.index_ptype)
            count = data_file.total_particles.get(self.index_ptype)
            account = [0] + [data_file.total_particles.get(ptype)
                             for ptype in self._ptypes]
            account = np.cumsum(account)
            return self._get_morton_from_position(
                data_file, account, account[idpos], regions, DLE, DRE)

    def _count_particles(self, data_file):
        npart = dict((self._ptypes[i], v)
                     for i, v in enumerate(data_file.header["Npart"]))
        return npart

    # header is 256, but we have 4 at beginning and end for ints
    _field_size = 4
    def _calculate_field_offsets(self, field_list, pcount,
                                 offset, file_size=None):
        # field_list is (ftype, fname) but the blocks are ordered
        # (fname, ftype) in the file.
        if self._format == 2:
            # Need to subtract offset due to extra header block
            pos = offset - SNAP_FORMAT_2_OFFSET
        else:
            pos = offset
        fs = self._field_size
        offsets = {}

        for field in self._fields:
            if not isinstance(field, string_types):
                field = field[0]
            if not any((ptype, field) in field_list
                       for ptype in self._ptypes):
                continue
            if self._format == 2:
                pos += 20  # skip block header
            elif self._format == 1:
                pos += 4
            else:
                raise RuntimeError(
                    "incorrect Gadget format %s!" % str(self._format))
            any_ptypes = False
            for ptype in self._ptypes:
                if field == "Mass" and ptype not in self.var_mass:
                    continue
                if (ptype, field) not in field_list:
                    continue
                offsets[(ptype, field)] = pos
                any_ptypes = True
                if field in self._vector_fields:
                    pos += self._vector_fields[field] * pcount[ptype] * fs
                else:
                    pos += pcount[ptype] * fs
            pos += 4
            if not any_ptypes:
                pos -= 8
        if file_size is not None:
            if (file_size != pos) & (self._format == 1):  # ignore the rest of format 2
                mylog.warning("Your Gadget-2 file may have extra " +
                              "columns or different precision!" +
                              " (%s file vs %s computed)",
                              file_size, pos)
        return offsets

    def _identify_fields(self, domain):
        # We can just look at the particle counts.
        field_list = []
        tp = domain.total_particles
        for i, ptype in enumerate(self._ptypes):
            count = tp[ptype]
            if count == 0:
                continue
            m = domain.header["Massarr"][i]
            for field in self._fields:
                if isinstance(field, tuple):
                    field, req = field
                    if req is ZeroMass:
                        if m > 0.0:
                            continue
                    elif isinstance(req, tuple) and ptype in req:
                        pass
                    elif req != ptype:
                        continue
                field_list.append((ptype, field))
        return field_list, {}
