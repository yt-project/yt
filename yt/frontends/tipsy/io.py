"""
Tipsy data-file handling function




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import numpy as np
from numpy.lib.recfunctions import append_fields
import os
import struct

from yt.geometry.particle_geometry_handler import \
    CHUNKSIZE
from yt.frontends.sph.io import \
    IOHandlerSPH
from yt.frontends.tipsy.definitions import \
    npart_mapping
from yt.utilities.lib.particle_kdtree_tools import \
    generate_smoothing_length
from yt.utilities.logger import ytLogger as \
    mylog


class IOHandlerTipsyBinary(IOHandlerSPH):
    _dataset_type = "tipsy"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")

    _pdtypes = None  # dtypes, to be filled in later
    _aux_pdtypes = None  # auxiliary files' dtypes

    _ptypes = ("Gas",
               "DarkMatter",
               "Stars")
    _chunksize = 64 * 64 * 64

    _aux_fields = None
    _fields = (("Gas", "Mass"),
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

    def __init__(self, *args, **kwargs):
        self._aux_fields = []
        super(IOHandlerTipsyBinary, self).__init__(*args, **kwargs)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _fill_fields(self, fields, vals, hsml, mask, data_file):
        if mask is None:
            size = 0
        elif isinstance(mask, slice):
            size = vals[fields[0]].size
        else:
            size = mask.sum()
        rv = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", size, field)
            if field in self._vector_fields:
                rv[field] = np.empty((size, 3), dtype="float64")
                if size == 0:
                    continue
                rv[field][:, 0] = vals[field]['x'][mask]
                rv[field][:, 1] = vals[field]['y'][mask]
                rv[field][:, 2] = vals[field]['z'][mask]
            elif field == 'smoothing_length':
                rv[field] = hsml[mask]
            else:
                rv[field] = np.empty(size, dtype="float64")
                if size == 0:
                    continue
                rv[field][:] = vals[field][mask]
            if field == "Coordinates":
                eps = np.finfo(rv[field].dtype).eps
                for i in range(3):
                    rv[field][:, i] = np.clip(
                        rv[field][:, i],
                        self.ds.domain_left_edge[i].v + eps,
                        self.ds.domain_right_edge[i].v - eps)
        return rv

    def _read_particle_coords(self, chunks, ptf):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items(),
                                            key=lambda a: poff[a[0]]):
                f.seek(poff[ptype])
                total = 0
                while total < tp[ptype]:
                    count = min(self._chunksize, tp[ptype] - total)
                    p = np.fromfile(f, self._pdtypes[ptype], count=count)
                    total += p.size
                    d = [p["Coordinates"][ax].astype("float64")
                         for ax in 'xyz']
                    del p
                    if ptype == self.ds._sph_ptype:
                        hsml = self._read_smoothing_length(data_file, count)
                    else:
                        hsml = 0.0
                    yield ptype, d, hsml

    @property
    def hsml_filename(self):
        return '%s-%s' % (self.ds.parameter_filename, 'hsml')

    def _generate_smoothing_length(self, data_files, kdtree):
        if os.path.exists(self.hsml_filename):
            with open(self.hsml_filename, 'rb') as f:
                file_hash = struct.unpack('q', f.read(struct.calcsize('q')))[0]
            if file_hash != self.ds._file_hash:
                os.remove(self.hsml_filename)
            else:
                return
        positions = []
        for data_file in data_files:
            for _, ppos in self._yield_coordinates(
                    data_file, needed_ptype=self.ds._sph_ptype):
                positions.append(ppos)
        if positions == []:
            return
        positions = np.concatenate(positions)[kdtree.idx]
        hsml = generate_smoothing_length(
            positions, kdtree, self.ds._num_neighbors)
        hsml = hsml[np.argsort(kdtree.idx)]
        dtype = self._pdtypes['Gas']['Coordinates'][0]
        with open(self.hsml_filename, 'wb') as f:
            f.write(struct.pack('q', self.ds._file_hash))
            f.write(hsml.astype(dtype).tostring())

    def _read_smoothing_length(self, data_file, count):
        dtype = self._pdtypes['Gas']['Coordinates'][0]
        with open(self.hsml_filename, 'rb') as f:
            f.seek(struct.calcsize('q') + data_file.start*dtype.itemsize)
            hsmls = np.fromfile(f, dtype, count=count)
        return hsmls.astype('float64')

    def _get_smoothing_length(self, data_file, dtype, shape):
        return self._read_smoothing_length(data_file, shape[0])

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            poff = data_file.field_offsets
            aux_fields_offsets = self._calculate_particle_offsets_aux(data_file)
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")

            # we need to open all aux files for chunking to work
            aux_fh = {}
            for afield in self._aux_fields:
                aux_fh[afield] = open(data_file.filename + '.' + afield, 'rb')

            for ptype, field_list in sorted(ptf.items(),
                                            key=lambda a: poff[a[0]]):
                f.seek(poff[ptype])
                afields = list(set(field_list).intersection(self._aux_fields))
                count = min(self._chunksize, tp[ptype])
                p = np.fromfile(f, self._pdtypes[ptype], count=count)
                auxdata = []
                if data_file.total_particles[ptype] == 0:
                    continue
                for afield in afields:
                    aux_fh[afield].seek(aux_fields_offsets[afield][ptype])
                    if isinstance(self._aux_pdtypes[afield], np.dtype):
                        auxdata.append(np.fromfile(
                            aux_fh[afield],
                            self._aux_pdtypes[afield],
                            count=count)
                        )
                    else:
                        par = self.ds.parameters
                        nlines = 1 + par['nsph'] + par['ndark'] + par['nstar']
                        aux_fh[afield].seek(0)
                        sh = aux_fields_offsets[afield][ptype]
                        sf = nlines - count - sh
                        if tp[ptype] > 0:
                            aux = np.genfromtxt(
                                aux_fh[afield], skip_header=sh,
                                skip_footer=sf
                            )
                            if aux.ndim < 1:
                                aux = np.array([aux])
                            auxdata.append(aux)
                if afields:
                    p = append_fields(p, afields, auxdata)
                if ptype == 'Gas':
                    hsml = self._read_smoothing_length(data_file, count)
                else:
                    hsml = 0.
                if getattr(selector, 'is_all_data', False):
                    mask = slice(None, None, None)
                else:
                    x = p["Coordinates"]['x'].astype("float64")
                    y = p["Coordinates"]['y'].astype("float64")
                    z = p["Coordinates"]['z'].astype("float64")
                    mask = selector.select_points(x, y, z, hsml)
                    del x, y, z
                if mask is None:
                    continue
                tf = self._fill_fields(field_list, p, hsml, mask, data_file)
                for field in field_list:
                    yield (ptype, field), tf.pop(field)

            # close all file handles
            f.close()
            for fh in list(aux_fh.values()):
                fh.close()

    def _update_domain(self, data_file):
        '''
        This method is used to determine the size needed for a box that will
        bound the particles.  It simply finds the largest position of the
        whole set of particles, and sets the domain to +/- that value.
        '''
        ds = data_file.ds
        ind = 0
        # Check to make sure that the domain hasn't already been set
        # by the parameter file
        if np.all(np.isfinite(ds.domain_left_edge)) and \
                np.all(np.isfinite(ds.domain_right_edge)):
            return
        with open(data_file.filename, "rb") as f:
            ds.domain_left_edge = 0
            ds.domain_right_edge = 0
            f.seek(ds._header_offset)
            mi = np.array([1e30, 1e30, 1e30], dtype="float64")
            ma = -np.array([1e30, 1e30, 1e30], dtype="float64")
            for iptype, ptype in enumerate(self._ptypes):
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0:
                    continue
                stop = ind + count
                while ind < stop:
                    c = min(CHUNKSIZE, stop - ind)
                    pp = np.fromfile(f, dtype=self._pdtypes[ptype],
                                     count=c)
                    np.minimum(mi, [pp["Coordinates"]["x"].min(),
                                    pp["Coordinates"]["y"].min(),
                                    pp["Coordinates"]["z"].min()], mi)
                    np.maximum(ma, [pp["Coordinates"]["x"].max(),
                                    pp["Coordinates"]["y"].max(),
                                    pp["Coordinates"]["z"].max()], ma)
                    ind += c
        # We extend by 1%.
        DW = ma - mi
        mi -= 0.01 * DW
        ma += 0.01 * DW
        ds.domain_left_edge = ds.arr(mi, 'code_length')
        ds.domain_right_edge = ds.arr(ma, 'code_length')
        ds.domain_width = DW = ds.domain_right_edge - ds.domain_left_edge
        ds.unit_registry.add("unitary", float(DW.max() * DW.units.base_value),
                             DW.units.dimensions)

    def _yield_coordinates(self, data_file, needed_ptype=None):
        with open(data_file.filename, "rb") as f:
            poff = data_file.field_offsets
            for iptype, ptype in enumerate(self._ptypes):
                if ptype not in poff:
                    continue
                f.seek(poff[ptype])
                if needed_ptype is not None and ptype != needed_ptype:
                    continue
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0:
                    continue
                pp = np.fromfile(f, dtype=self._pdtypes[ptype],
                                 count=count)
                mis = np.empty(3, dtype="float64")
                mas = np.empty(3, dtype="float64")
                for axi, ax in enumerate('xyz'):
                    mi = pp["Coordinates"][ax].min()
                    ma = pp["Coordinates"][ax].max()
                    mylog.debug(
                        "Spanning: %0.3e .. %0.3e in %s", mi, ma, ax)
                    mis[axi] = mi
                    mas[axi] = ma
                pos = np.empty((pp.size, 3), dtype="float64")
                for i, ax in enumerate("xyz"):
                    pos[:,i] = pp["Coordinates"][ax]
                yield ptype, pos

    def _count_particles(self, data_file):
        pcount = np.array([data_file.ds.parameters['nsph'],
                           data_file.ds.parameters['nstar'],
                           data_file.ds.parameters['ndark']])
        si, ei = data_file.start, data_file.end
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
        ptypes = ['Gas', 'Stars', 'DarkMatter']
        npart = dict((ptype, v) for ptype, v in zip(ptypes, pcount))
        return npart

    @classmethod
    def _compute_dtypes(cls, field_dtypes, endian="<"):
        pds = {}
        for ptype, field in cls._fields:
            dtbase = field_dtypes.get(field, 'f')
            ff = "%s%s" % (endian, dtbase)
            if field in cls._vector_fields:
                dt = (field, [('x', ff), ('y', ff), ('z', ff)])
            else:
                dt = (field, ff)
            pds.setdefault(ptype, []).append(dt)
        pdtypes = {}
        for ptype in pds:
            pdtypes[ptype] = np.dtype(pds[ptype])
        return pdtypes

    def _create_dtypes(self, data_file):
        # We can just look at the particle counts.
        self._header_offset = data_file.ds._header_offset
        self._pdtypes = self._compute_dtypes(data_file.ds._field_dtypes,
                                             data_file.ds.endian)
        self._field_list = []
        for ptype, field in self._fields:
            if data_file.total_particles[ptype] == 0:
                # We do not want out _pdtypes to have empty particles.
                self._pdtypes.pop(ptype, None)
                continue
            self._field_list.append((ptype, field))

        if 'Gas' in self._pdtypes.keys():
            self._field_list.append(('Gas', 'smoothing_length'))

        # Find out which auxiliaries we have and what is their format
        tot_parts = np.sum(
            [data_file.ds.parameters['nsph'],
             data_file.ds.parameters['nstar'],
             data_file.ds.parameters['ndark']])
        endian = data_file.ds.endian
        self._aux_pdtypes = {}
        self._aux_fields = []
        for f in glob.glob(data_file.filename + '.*'):
            afield = f.rsplit('.')[-1]
            filename = data_file.filename + '.' + afield
            if not os.path.exists(filename):
                continue
            if afield in ['log', 'parameter', 'kdtree']:
                # Amiga halo finder makes files like this we need to ignore
                continue
            self._aux_fields.append(afield)
        skip_afields = []
        for afield in self._aux_fields:
            filename = data_file.filename + '.' + afield
            # We need to do some fairly ugly detection to see what format the
            # auxiliary files are in.  They can be either ascii or binary, and
            # the binary files can be either floats, ints, or doubles.  We're
            # going to use a try-catch cascade to determine the format.
            filesize = os.stat(filename).st_size
            dtype = np.dtype(endian + 'i4')
            tot_parts_from_file = np.fromfile(filename, dtype, count=1)
            if tot_parts_from_file != tot_parts:
                with open(filename, 'rb') as f:
                    header_nparts = f.readline()
                    try:
                        header_nparts = int(header_nparts)
                    except ValueError:
                        skip_afields.append(afield)
                        continue
                    if int(header_nparts) != tot_parts:
                        raise RuntimeError
                self._aux_pdtypes[afield] = "ascii"
            elif (filesize - 4) / 8 == tot_parts:
                self._aux_pdtypes[afield] = np.dtype([('aux', endian + 'd')])
            elif (filesize - 4) / 4 == tot_parts:
                if afield.startswith("i"):
                    self._aux_pdtypes[afield] = np.dtype([('aux', endian + 'i')])
                else:
                    self._aux_pdtypes[afield] = np.dtype([('aux', endian + 'f')])
            else:
                skip_afields.append(afield)
        for afield in skip_afields:
            self._aux_fields.remove(afield)
        # Add the auxiliary fields to each ptype we have
        for ptype in self._ptypes:
            if any([ptype == field[0] for field in self._field_list]):
                self._field_list += \
                    [(ptype, afield) for afield in self._aux_fields]
        return self._field_list

    def _identify_fields(self, data_file):
        return self._field_list, {}

    def _calculate_particle_offsets(self, data_file, pcounts):
        field_offsets = {}
        pos = data_file.ds._header_offset
        for i, ptype in enumerate(self._ptypes):
            if ptype not in self._pdtypes:
                continue
            size = self._pdtypes[ptype].itemsize
            field_offsets[ptype] = pos
            params = self.ds.parameters
            if params[npart_mapping[ptype]] > CHUNKSIZE:
                 for j in range(i):
                     npart = params[npart_mapping[self._ptypes[j]]]
                     spart = self._pdtypes[self._ptypes[j]].itemsize
                     if npart > CHUNKSIZE:
                         field_offsets[ptype] += npart*spart
                 field_offsets[ptype] += data_file.start*size
            pos += data_file.total_particles[ptype] * size
        return field_offsets

    def _calculate_particle_offsets_aux(self, data_file):
        aux_fields_offsets = {}
        params = self.ds.parameters
        for afield in self._aux_fields:
            aux_fields_offsets[afield] = {}
            if isinstance(self._aux_pdtypes[afield], np.dtype):
                pos = 4  # i4
                size = np.dtype(self._aux_pdtypes[afield]).itemsize
            else:
                pos = 1
                size = 1
            for i, ptype in enumerate(self._ptypes):
                if data_file.total_particles[ptype] == 0:
                    continue
                elif params[npart_mapping[ptype]] > CHUNKSIZE:
                    for j in range(i):
                        npart = params[npart_mapping[self._ptypes[j]]]
                        if npart > CHUNKSIZE:
                            pos += npart*size
                    pos += data_file.start*size
                aux_fields_offsets[afield][ptype] = pos
                pos += data_file.total_particles[ptype] * size
        return aux_fields_offsets
