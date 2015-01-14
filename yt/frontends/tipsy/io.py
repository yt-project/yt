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
import os

from yt.geometry.oct_container import \
    _ORDER_MAX
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

    _pdtypes = None # dtypes, to be filled in later

    _ptypes = ( "Gas",
                "DarkMatter",
                "Stars" )
    _chunksize = 64*64*64

    _aux_fields = None
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

    def __init__(self, *args, **kwargs):
        self._aux_fields = []
        super(IOHandlerTipsyBinary, self).__init__(*args, **kwargs)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_aux_fields(self, field, mask, data_file):
        """
        Read in auxiliary files from gasoline/pkdgrav.
        This method will automatically detect the format of the file.
        """
        filename = data_file.filename+'.'+field
        dtype = None
        # We need to do some fairly ugly detection to see what format the auxiliary
        # files are in.  They can be either ascii or binary, and the binary files can be
        # either floats, ints, or doubles.  We're going to use a try-catch cascade to
        # determine the format.
        try:#ASCII
            auxdata = np.genfromtxt(filename, skip_header=1)
            if auxdata.size != np.sum(data_file.total_particles.values()):
                print("Error reading auxiliary tipsy file")
                raise RuntimeError 
        except ValueError:#binary/xdr
            f = open(filename, 'rb')
            l = struct.unpack(data_file.ds.endian+"i", f.read(4))[0]
            if l != np.sum(data_file.total_particles.values()):
                print("Error reading auxiliary tipsy file")
                raise RuntimeError
            dtype = 'd'
            if field in ('iord', 'igasorder', 'grp'):#These fields are integers
                dtype = 'i'
            try:# If we try loading doubles by default, we can catch an exception and try floats next
                auxdata = np.array(struct.unpack(data_file.ds.endian+(l*dtype), f.read()))
            except struct.error:
                f.seek(4)
                dtype = 'f'
                try:
                    auxdata = np.array(struct.unpack(data_file.ds.endian+(l*dtype), f.read()))
                except struct.error: # None of the binary attempts to read succeeded
                    print("Error reading auxiliary tipsy file")
                    raise RuntimeError

        # Use the mask to slice out the appropriate particle type data
        if mask.size == data_file.total_particles['Gas']:
            return auxdata[:data_file.total_particles['Gas']]
        elif mask.size == data_file.total_particles['DarkMatter']:
            return auxdata[data_file.total_particles['Gas']:-data_file.total_particles['DarkMatter']]
        else:
            return auxdata[-data_file.total_particles['Stars']:]

    def _fill_fields(self, fields, vals, mask, data_file):
        if mask is None:
            size = 0
        else:
            size = mask.sum()
        rv = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", size, field)
            if field in self._aux_fields: #Read each of the auxiliary fields
                rv[field] = self._read_aux_fields(field, mask, data_file)
            elif field in self._vector_fields:
                rv[field] = np.empty((size, 3), dtype="float64")
                if size == 0: continue
                rv[field][:,0] = vals[field]['x'][mask]
                rv[field][:,1] = vals[field]['y'][mask]
                rv[field][:,2] = vals[field]['z'][mask]
            else:
                rv[field] = np.empty(size, dtype="float64")
                if size == 0: continue
                rv[field][:] = vals[field][mask]
            if field == "Coordinates":
                eps = np.finfo(rv[field].dtype).eps
                for i in range(3):
                  rv[field][:,i] = np.clip(rv[field][:,i],
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
            for ptype, field_list in sorted(ptf.items(), key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                total = 0
                while total < tp[ptype]:
                    p = np.fromfile(f, self._pdtypes[ptype],
                            count=min(self._chunksize, tp[ptype] - total))
                    total += p.size
                    d = [p["Coordinates"][ax].astype("float64") for ax in 'xyz']
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
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items(), key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                total = 0
                while total < tp[ptype]:
                    p = np.fromfile(f, self._pdtypes[ptype],
                        count=min(self._chunksize, tp[ptype] - total))
                    total += p.size
                    mask = selector.select_points(
                        p["Coordinates"]['x'].astype("float64"),
                        p["Coordinates"]['y'].astype("float64"),
                        p["Coordinates"]['z'].astype("float64"), 0.0)
                    if mask is None: continue
                    tf = self._fill_fields(field_list, p, mask, data_file)
                    for field in field_list:
                        yield (ptype, field), tf.pop(field)
            f.close()

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
        if np.all(np.isfinite(ds.domain_left_edge)) and np.all(np.isfinite(ds.domain_right_edge)):
            return
        with open(data_file.filename, "rb") as f:
            ds.domain_left_edge = 0
            ds.domain_right_edge = 0
            f.seek(ds._header_offset)
            mi =   np.array([1e30, 1e30, 1e30], dtype="float64")
            ma =  -np.array([1e30, 1e30, 1e30], dtype="float64")
            for iptype, ptype in enumerate(self._ptypes):
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0: continue
                start, stop = ind, ind + count
                while ind < stop:
                    c = min(CHUNKSIZE, stop - ind)
                    pp = np.fromfile(f, dtype = self._pdtypes[ptype],
                                     count = c)
                    eps = np.finfo(pp["Coordinates"]["x"].dtype).eps
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
        ds.unit_registry.add("unitary", float(DW.max() * DW.units.cgs_value),
                                 DW.units.dimensions)

    def _initialize_index(self, data_file, regions):
        ds = data_file.ds
        morton = np.empty(sum(data_file.total_particles.values()),
                          dtype="uint64")
        ind = 0
        DLE, DRE = ds.domain_left_edge, ds.domain_right_edge
        dx = (DRE - DLE) / (2**_ORDER_MAX)
        self.domain_left_edge = DLE.in_units("code_length").ndarray_view()
        self.domain_right_edge = DRE.in_units("code_length").ndarray_view()
        with open(data_file.filename, "rb") as f:
            f.seek(ds._header_offset)
            for iptype, ptype in enumerate(self._ptypes):
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0: continue
                start, stop = ind, ind + count
                while ind < stop:
                    c = min(CHUNKSIZE, stop - ind)
                    pp = np.fromfile(f, dtype = self._pdtypes[ptype],
                                     count = c)
                    mis = np.empty(3, dtype="float64")
                    mas = np.empty(3, dtype="float64")
                    for axi, ax in enumerate('xyz'):
                        mi = pp["Coordinates"][ax].min()
                        ma = pp["Coordinates"][ax].max()
                        mylog.debug("Spanning: %0.3e .. %0.3e in %s", mi, ma, ax)
                        mis[axi] = mi
                        mas[axi] = ma
                    pos = np.empty((pp.size, 3), dtype="float64")
                    for i, ax in enumerate("xyz"):
                        eps = np.finfo(pp["Coordinates"][ax].dtype).eps
                        pos[:,i] = pp["Coordinates"][ax]
                    regions.add_data_file(pos, data_file.file_id,
                                          data_file.ds.filter_bbox)
                    morton[ind:ind+c] = compute_morton(
                        pos[:,0], pos[:,1], pos[:,2],
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
    def _compute_dtypes(cls, field_dtypes, endian = "<"):
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
        self._pdtypes = {}
        pds = {}
        field_list = []
        tp = data_file.total_particles
        aux_filenames = glob.glob(data_file.filename+'.*') # Find out which auxiliaries we have
        self._aux_fields = [f[1+len(data_file.filename):] for f in aux_filenames]
        self._pdtypes = self._compute_dtypes(data_file.ds._field_dtypes,
                                             data_file.ds.endian)
        for ptype, field in self._fields:
            if tp[ptype] == 0:
                # We do not want out _pdtypes to have empty particles.
                self._pdtypes.pop(ptype, None)
                continue
            field_list.append((ptype, field))
        if any(["Gas"==f[0] for f in field_list]): #Add the auxiliary fields to each ptype we have
            field_list += [("Gas",a) for a in self._aux_fields]
        if any(["DarkMatter"==f[0] for f in field_list]):
            field_list += [("DarkMatter",a) for a in self._aux_fields]
        if any(["Stars"==f[0] for f in field_list]):
            field_list += [("Stars",a) for a in self._aux_fields]
        self._field_list = field_list
        return self._field_list

    def _identify_fields(self, data_file):
        return self._field_list, {}

    def _calculate_particle_offsets(self, data_file):
        field_offsets = {}
        pos = data_file.ds._header_offset
        for ptype in self._ptypes:
            field_offsets[ptype] = pos
            if data_file.total_particles[ptype] == 0: continue
            size = self._pdtypes[ptype].itemsize
            pos += data_file.total_particles[ptype] * size
        return field_offsets
