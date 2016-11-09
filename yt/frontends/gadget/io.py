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
from yt.frontends.owls.io import \
    IOHandlerOWLS
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton
from yt.utilities.logger import ytLogger as mylog

class IOHandlerGadgetHDF5(IOHandlerOWLS):
    _dataset_type = "gadget_hdf5"

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

    def __init__(self, ds, *args, **kwargs):
        self._vector_fields = dict(self._vector_fields)
        self._fields = ds._field_spec
        self._ptypes = ds._ptype_spec
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
                yield ptype, (pos[:,0], pos[:,1], pos[:,2])
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
                    pos[:,0], pos[:,1], pos[:,2], 0.0)
                del pos
                if mask is None: continue
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
                    data = data[mask,...]
                    yield (ptype, field), data
            f.close()

    def _read_field_from_file(self, f, count, name):
        if count == 0: return
        if name == "ParticleIDs":
            dt = "uint32"
        else:
            dt = "float32"
        if name in self._vector_fields:
            count *= self._vector_fields[name]
        arr = np.fromfile(f, dtype=dt, count = count)
        if name in self._vector_fields:
            factor = self._vector_fields[name]
            arr = arr.reshape((count//factor, factor), order="C")
        return arr.astype("float64")

    def _initialize_index(self, data_file, regions):
        count = sum(data_file.total_particles.values())
        DLE = data_file.ds.domain_left_edge
        DRE = data_file.ds.domain_right_edge
        with open(data_file.filename, "rb") as f:
            # We add on an additionally 4 for the first record.
            f.seek(data_file._position_offset + 4)
            # The first total_particles * 3 values are positions
            pp = np.fromfile(f, dtype = 'float32', count = count*3)
            pp.shape = (count, 3)
        regions.add_data_file(pp, data_file.file_id, data_file.ds.filter_bbox)
        morton = compute_morton(pp[:,0], pp[:,1], pp[:,2], DLE, DRE,
                                data_file.ds.filter_bbox)
        return morton

    def _count_particles(self, data_file):
        npart = dict((self._ptypes[i], v)
            for i, v in enumerate(data_file.header["Npart"]))
        return npart

    # header is 256, but we have 4 at beginning and end for ints
    _field_size = 4
    def _calculate_field_offsets(self, field_list, pcount,
                                 offset, file_size = None):
        # field_list is (ftype, fname) but the blocks are ordered
        # (fname, ftype) in the file.
        pos = offset
        fs = self._field_size
        offsets = {}
        for field in self._fields:
            if not isinstance(field, string_types):
                field = field[0]
            if not any( (ptype, field) in field_list
                        for ptype in self._ptypes):
                continue
            pos += 4
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
            if not any_ptypes: pos -= 8
        if file_size is not None:
            if file_size != pos:
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
            if count == 0: continue
            m = domain.header["Massarr"][i]
            for field in self._fields:
                if isinstance(field, tuple):
                    field, req = field
                    if req is ZeroMass:
                        if m > 0.0 : continue
                    elif isinstance(req, tuple) and ptype in req:
                        pass
                    elif req != ptype:
                        continue
                field_list.append((ptype, field))
        return field_list, {}
