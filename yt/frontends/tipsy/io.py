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

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton
from yt.utilities.logger import ytLogger as \
    mylog

CHUNKSIZE = 10000000


class IOHandlerTipsyBinary(BaseIOHandler):
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

    def _fill_fields(self, fields, vals, mask, data_file):
        if mask is None:
            size = 0
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
            else:
                rv[field] = np.empty(size, dtype="float64")
                if size == 0:
                    continue
                rv[field][:] = vals[field][mask]
            if field == "Coordinates":
                eps = np.finfo(rv[field].dtype).eps
                for i in range(3):
                    rv[field][:, i] = np.clip(rv[field][:, i],
                                              self.domain_left_edge[i] + eps,
                                              self.domain_right_edge[i] - eps)
        return rv

    def _read_particle_coords(self, chunks, ptf):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items(),
                                            key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                total = 0
                while total < tp[ptype]:
                    count = min(self._chunksize, tp[ptype] - total)
                    p = np.fromfile(f, self._pdtypes[ptype], count=count)
                    total += p.size
                    d = [p["Coordinates"][ax].astype("float64")
                         for ax in 'xyz']
                    del p
                    yield ptype, d

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            aux_fields_offsets = \
                self._calculate_particle_offsets_aux(data_file)
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")

            # we need to open all aux files for chunking to work
            aux_fh = {}
            for afield in self._aux_fields:
                aux_fh[afield] = open(data_file.filename + '.' + afield, 'rb')

            for ptype, field_list in sorted(ptf.items(),
                                            key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                afields = list(set(field_list).intersection(self._aux_fields))
                for afield in afields:
                    aux_fh[afield].seek(
                        aux_fields_offsets[afield][ptype][0], os.SEEK_SET)

                total = 0
                while total < tp[ptype]:
                    count = min(self._chunksize, tp[ptype] - total)
                    p = np.fromfile(f, self._pdtypes[ptype], count=count)

                    auxdata = []
                    for afield in afields:
                        if isinstance(self._aux_pdtypes[afield], np.dtype):
                            auxdata.append(
                                np.fromfile(aux_fh[afield],
                                            self._aux_pdtypes[afield],
                                            count=count)
                            )
                        else:
                            aux_fh[afield].seek(0, os.SEEK_SET)
                            sh = aux_fields_offsets[afield][ptype][0] + total
                            sf = aux_fields_offsets[afield][ptype][1] + \
                                tp[ptype] - count - total
                            if tp[ptype] > 0:
                                aux = np.genfromtxt(
                                    aux_fh[afield], skip_header=sh,
                                    skip_footer=sf
                                )
                                if aux.ndim < 1:
                                    aux = np.array([aux])
                                auxdata.append(aux)

                    total += p.size
                    if afields:
                        p = append_fields(p, afields, auxdata)
                    mask = selector.select_points(
                        p["Coordinates"]['x'].astype("float64"),
                        p["Coordinates"]['y'].astype("float64"),
                        p["Coordinates"]['z'].astype("float64"), 0.0)
                    if mask is None:
                        continue
                    tf = self._fill_fields(field_list, p, mask, data_file)
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

    def _initialize_index(self, data_file, regions):
        ds = data_file.ds
        morton = np.empty(sum(list(data_file.total_particles.values())),
                          dtype="uint64")
        ind = 0
        DLE, DRE = ds.domain_left_edge, ds.domain_right_edge
        self.domain_left_edge = DLE.in_units("code_length").ndarray_view()
        self.domain_right_edge = DRE.in_units("code_length").ndarray_view()
        with open(data_file.filename, "rb") as f:
            f.seek(ds._header_offset)
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
                        pos[:, i] = pp["Coordinates"][ax]
                    regions.add_data_file(pos, data_file.file_id,
                                          data_file.ds.filter_bbox)
                    morton[ind:ind + c] = compute_morton(
                        pos[:, 0], pos[:, 1], pos[:, 2],
                        DLE, DRE, data_file.ds.filter_bbox)
                    ind += c
        mylog.info("Adding %0.3e particles", morton.size)
        return morton

    def _count_particles(self, data_file):
        npart = {
            "Gas": data_file.ds.parameters['nsph'],
            "Stars": data_file.ds.parameters['nstar'],
            "DarkMatter": data_file.ds.parameters['ndark']
        }
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

        # Find out which auxiliaries we have and what is their format
        tot_parts = np.sum(list(data_file.total_particles.values()))
        endian = data_file.ds.endian
        self._aux_pdtypes = {}
        self._aux_fields = []
        for f in glob.glob(data_file.filename + '.*'):
            afield = f.rsplit('.')[-1]
            filename = data_file.filename + '.' + afield
            if not os.path.exists(filename): continue
            self._aux_fields.append(afield)
        for afield in self._aux_fields:
            filename = data_file.filename + '.' + afield
            # We need to do some fairly ugly detection to see what format the
            # auxiliary files are in.  They can be either ascii or binary, and
            # the binary files can be either floats, ints, or doubles.  We're
            # going to use a try-catch cascade to determine the format.
            filesize = os.stat(filename).st_size
            if np.fromfile(filename, np.dtype(endian + 'i4'),
                           count=1) != tot_parts:
                with open(filename, 'rb') as f:
                    header_nparts = f.readline()
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
                raise RuntimeError

        # Add the auxiliary fields to each ptype we have
        for ptype in self._ptypes:
            if any([ptype == field[0] for field in self._field_list]):
                self._field_list += \
                    [(ptype, afield) for afield in self._aux_fields]
        return self._field_list

    def _identify_fields(self, data_file):
        return self._field_list, {}

    def _calculate_particle_offsets(self, data_file):
        field_offsets = {}
        pos = data_file.ds._header_offset
        for ptype in self._ptypes:
            field_offsets[ptype] = pos
            if data_file.total_particles[ptype] == 0:
                continue
            size = self._pdtypes[ptype].itemsize
            pos += data_file.total_particles[ptype] * size
        return field_offsets

    def _calculate_particle_offsets_aux(self, data_file):
        aux_fields_offsets = {}
        tp = data_file.total_particles
        for afield in self._aux_fields:
            aux_fields_offsets[afield] = {}
            if isinstance(self._aux_pdtypes[afield], np.dtype):
                pos = 4  # i4
                for ptype in self._ptypes:
                    aux_fields_offsets[afield][ptype] = (pos, 0)
                    if data_file.total_particles[ptype] == 0:
                        continue
                    size = np.dtype(self._aux_pdtypes[afield]).itemsize
                    pos += data_file.total_particles[ptype] * size
            else:
                aux_fields_offsets[afield].update(
                    {'Gas': (1, tp["DarkMatter"] + tp["Stars"]),
                     'DarkMatter': (1 + tp["Gas"], tp["Stars"]),
                     'Stars': (1 + tp["DarkMatter"] + tp["Gas"], 0)}
                )
        return aux_fields_offsets
