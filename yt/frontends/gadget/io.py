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
from yt.frontends.sph.io import \
    IOHandlerSPH
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .definitions import \
    gadget_hdf5_ptypes, \
    SNAP_FORMAT_2_OFFSET


class IOHandlerGadgetHDF5(IOHandlerSPH):
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
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                x = f["/%s/Coordinates" % ptype][si:ei, 0].astype("float64")
                y = f["/%s/Coordinates" % ptype][si:ei, 1].astype("float64")
                z = f["/%s/Coordinates" % ptype][si:ei, 2].astype("float64")
                if ptype == self.ds._sph_ptype:
                    pdtype = f["/%s/Coordinates" % ptype].dtype
                    pshape = f["/%s/Coordinates" % ptype].shape
                    hsml = self._get_smoothing_length(data_file, pdtype, pshape)
                else:
                    hsml = 0.0
                yield ptype, (x, y, z), hsml
            f.close()

    def _yield_coordinates(self, data_file, needed_ptype=None):
        si, ei = data_file.start, data_file.end
        f = h5py.File(data_file.filename, "r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        np.clip(pcount - si, 0, ei - si, out=pcount)
        pcount = pcount.sum()
        for key in f.keys():
            if not key.startswith("PartType"):
                continue
            if "Coordinates" not in f[key]:
                continue
            if needed_ptype and key != needed_ptype:
                continue
            ds = f[key]["Coordinates"][si:ei,...]
            dt = ds.dtype.newbyteorder("N") # Native
            pos = np.empty(ds.shape, dtype=dt)
            pos[:] = ds
            yield key, pos
        f.close()

    def _get_smoothing_length(self, data_file, position_dtype, position_shape):
        ptype = self.ds._sph_ptype
        ind = int(ptype[-1])
        si, ei = data_file.start, data_file.end
        with h5py.File(data_file.filename, "r") as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][ind].astype("int")
            pcount = np.clip(pcount - si, 0, ei - si)
            if self._dataset_type == "arepo_hdf5":
                ds = f[ptype]["Masses"][si:ei,...]/f[ptype]["Density"][si:ei,...]
                ds *= 3.0/(4.0*np.pi)
                ds **= (1./3.)
                ds *= self.ds.smoothing_factor
            else:
                ds = f[ptype]["SmoothingLength"][si:ei,...]
            dt = ds.dtype.newbyteorder("N") # Native
            if position_dtype is not None and dt < position_dtype:
                # Sometimes positions are stored in double precision
                # but smoothing lengths are stored in single precision.
                # In these cases upcast smoothing length to double precision
                # to avoid ValueErrors when we pass these arrays to Cython.
                dt = position_dtype
            hsml = np.empty(ds.shape, dtype=dt)
            hsml[:] = ds
            return hsml

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                g = f["/%s" % ptype]
                coords = g["Coordinates"][si:ei].astype("float64")
                if ptype == 'PartType0':
                    hsmls = self._get_smoothing_length(data_file,
                                                       g["Coordinates"].dtype, 
                                                       g["Coordinates"].shape)
                else:
                    hsmls = 0.0
                if getattr(selector, 'is_all_data', False):
                    mask = slice(None, None, None)
                else:
                    mask = selector.select_points(
                        coords[:,0], coords[:,1], coords[:,2], hsmls)
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
                        data = g[rfield][si:ei][mask, ...]
                    elif field.startswith("Metallicity_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["Metallicity"][si:ei, col][mask]
                    elif field.startswith("Chemistry_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["ChemistryAbundances"][si:ei, col][mask]
                    elif field == "smoothing_length":
                        data = hsmls[mask].astype("float64")
                    else:
                        data = g[field][si:ei][mask, ...]

                    yield (ptype, field), data
            f.close()

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        f = h5py.File(data_file.filename, "r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        f.close()
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
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
        if self._dataset_type == "arepo_hdf5":
            fields.append(("PartType0", "smoothing_length"))
        f.close()
        return fields, {}

class IOHandlerArepoHDF5(IOHandlerGadgetHDF5):
    _dataset_type = "arepo_hdf5"


ZeroMass = object()


class IOHandlerGadgetBinary(IOHandlerSPH):
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
        gformat, endianswap = ds._header.gadget_format
        # gadget format 1 original, 2 with block name
        self._format = gformat
        self._endian = endianswap
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
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype in ptf:
                f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                pos = self._read_field_from_file(
                    f, tp[ptype], "Coordinates")
                if ptype == self.ds._sph_ptype:
                    f.seek(poff[ptype, "SmoothingLength"], os.SEEK_SET)
                    hsml = self._read_field_from_file(
                        f, tp[ptype], "SmoothingLength")
                else:
                    hsml = 0.0
                yield ptype, (pos[:, 0], pos[:, 1], pos[:, 2]), hsml
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items()):
                if getattr(selector, 'is_all_data', False):
                    mask = slice(None, None, None)
                else:
                    f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                    pos = self._read_field_from_file(
                        f, tp[ptype], "Coordinates")
                    if ptype == self.ds._sph_ptype:
                        f.seek(poff[ptype, "SmoothingLength"], os.SEEK_SET)
                        hsml = self._read_field_from_file(
                            f, tp[ptype], "SmoothingLength")
                    else:
                        hsml = 0.0
                    mask = selector.select_points(
                        pos[:, 0], pos[:, 1], pos[:, 2], hsml)
                    del pos
                    del hsml
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
            dt = self._endian + self.ds._id_dtype
        else:
            dt = self._endian + self._float_type
        dt = np.dtype(dt)
        if name in self._vector_fields:
            count *= self._vector_fields[name]
        arr = np.fromfile(f, dtype=dt, count=count)
        # ensure data are in native endianness to avoid errors
        # when field data are passed to cython
        dt = dt.newbyteorder('N')
        arr = arr.astype(dt)
        if name in self._vector_fields:
            factor = self._vector_fields[name]
            arr = arr.reshape((count // factor, factor), order="C")
        return arr

    def _yield_coordinates(self, data_file, needed_ptype=None):
        self._float_type = data_file.ds._header.float_type
        self._field_size = np.dtype(self._float_type).itemsize
        with open(data_file.filename, "rb") as f:
            # We add on an additionally 4 for the first record.
            f.seek(data_file._position_offset + 4)
            for ptype, count in data_file.total_particles.items():
                if count == 0:
                    continue
                if needed_ptype is not None and ptype != needed_ptype:
                    continue
                # The first total_particles * 3 values are positions
                pp = np.fromfile(f, dtype = self._float_type, count = count*3)
                pp.shape = (count, 3)
                yield ptype, pp

    def _get_smoothing_length(self, data_file, position_dtype, position_shape):
        ret = self._get_field(data_file, 'SmoothingLength', 'Gas')
        if position_dtype is not None and ret.dtype != position_dtype:
            # Sometimes positions are stored in double precision
            # but smoothing lengths are stored in single precision.
            # In these cases upcast smoothing length to double precision
            # to avoid ValueErrors when we pass these arrays to Cython.
            ret = ret.astype(position_dtype)
        return ret

    def _get_field(self, data_file, field, ptype):
        poff = data_file.field_offsets
        tp = data_file.total_particles
        with open(data_file.filename, "rb") as f:
            f.seek(poff[ptype, field], os.SEEK_SET)
            pp = self._read_field_from_file(
                f, tp[ptype], field)
        return pp

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        pcount = np.array(data_file.header["Npart"])
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
        npart = dict((self._ptypes[i], v) for i, v in enumerate(pcount))
        return npart

    # header is 256, but we have 4 at beginning and end for ints
    _field_size = 4
    def _calculate_field_offsets(self, field_list, pcount,
                                 offset, df_start, file_size=None):
        # field_list is (ftype, fname) but the blocks are ordered
        # (fname, ftype) in the file.
        if self._format == 2:
            # Need to subtract offset due to extra header block
            pos = offset - SNAP_FORMAT_2_OFFSET
        else:
            pos = offset
        fs = self._field_size
        offsets = {}
        pcount = dict(zip(self._ptypes, pcount))

        for field in self._fields:
            if field == "ParticleIDs" and self.ds.long_ids:
                fs = 8
            else:
                fs = 4
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
                start_offset = df_start * fs
                if field in self._vector_fields:
                    start_offset *= self._vector_fields[field]
                pos += start_offset
                offsets[(ptype, field)] = pos
                any_ptypes = True
                remain_offset = (pcount[ptype] - df_start) * fs
                if field in self._vector_fields:
                    remain_offset *= self._vector_fields[field]
                pos += remain_offset
            pos += 4
            if not any_ptypes:
                pos -= 8
        if file_size is not None:
            if (file_size != pos) & (self._format == 1): #ignore the rest of format 2 
                diff = file_size - pos
                possible = []
                for ptype, psize in sorted(pcount.items()):
                    if psize == 0: continue
                    if float(diff) / psize == int(float(diff)/psize):
                        possible.append(ptype)
                mylog.warning("Your Gadget-2 file may have extra " +
                              "columns or different precision! " +
                              "(%s diff => %s?)", diff, possible)
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
