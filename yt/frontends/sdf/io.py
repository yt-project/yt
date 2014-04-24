"""
SDF data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import h5py
import numpy as np
from yt.funcs import *
from yt.utilities.exceptions import *

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.fortran_utils import read_record
from yt.utilities.lib.geometry_utils import compute_morton

from yt.geometry.oct_container import _ORDER_MAX
CHUNKSIZE = 32**3

class IOHandlerSDF(BaseIOHandler):
    _dataset_type = "sdf_particles"

    @property
    def _handle(self):
        return self.pf.sdf_container

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in data_files:
            pcount = self._handle['x'].size
            yield "dark_matter", (
                self._handle['x'].astype("float64"),
                self._handle['y'].astype("float64"),
                self._handle['z'].astype("float64"))

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in data_files:
            pcount = self._handle['x'].size
            for ptype, field_list in sorted(ptf.items()):
                x = self._handle['x'].astype("float64")
                y = self._handle['y'].astype("float64")
                z = self._handle['z'].astype("float64")
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None: continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.pf.parameters["particle_mass"]
                    else:
                        data = self._handle[field][mask].astype("float64")
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        x, y, z = (self._handle[ax] for ax in 'xyz')
        pcount = x.size
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        pos = np.empty((CHUNKSIZE, 3), dtype=x.dtype)
        while ind < pcount:
            npart = min(CHUNKSIZE, pcount - ind)
            pos[:npart,0] = x[ind:ind+npart]
            pos[:npart,1] = y[ind:ind+npart]
            pos[:npart,2] = z[ind:ind+npart]
            if np.any(pos.min(axis=0) < self.pf.domain_left_edge) or \
               np.any(pos.max(axis=0) > self.pf.domain_right_edge):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       self.pf.domain_left_edge,
                                       self.pf.domain_right_edge)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+npart] = compute_morton(
                pos[:npart,0], pos[:npart,1], pos[:npart,2],
                data_file.pf.domain_left_edge,
                data_file.pf.domain_right_edge)
            ind += CHUNKSIZE
        return morton

    def _count_particles(self, data_file):
        return {'dark_matter': self._handle['x'].size}

    def _identify_fields(self, data_file):
        fields = [("dark_matter", v) for v in self._handle.keys()]
        fields.append(("dark_matter", "mass"))
        return fields, {}

import re
import os

_types = {
    'int': 'int32',
    'int64_t': 'int64',
    'float': 'float32',
    'double': 'float64',
    'unsigned int': 'I',
    'unsigned char': 'B',
}

def get_type(vtype, len=None):
    try:
        t = _types[vtype]
        if len is not None:
            t = np.dtype((t, len))
        else:
            t = np.dtype(t)
    except KeyError:
        t = eval("np."+vtype)
    return t

def lstrip(text_list):
    return [t.strip() for t in text_list]

def get_struct_vars(line):
    spl = lstrip(line.split(";"))
    multiv = lstrip(spl[0].split(","))
    ret = lstrip(multiv[0].split())
    ctype = ret[0]
    vnames = [ret[-1]] + multiv[1:]
    vnames = [v.strip() for v in vnames]
    for vtype in ret[1:-1]:
        ctype += ' ' + vtype
    num = None
    if len(vnames) == 1:
        if '[' in vnames[0]:
            num = int(vnames[0].split('[')[-1].strip(']'))
            #num = int(re.sub("\D", "", vnames[0]))
    ctype = get_type(ctype, len=num)
    return ctype, vnames

class DataStruct(object):
    """docstring for DataStruct"""

    _offset = 0

    def __init__(self, dtypes, num, filename):
        self.filename = filename
        self.dtype = np.dtype(dtypes)
        self.size = num
        self.itemsize = self.dtype.itemsize
        self.data = {}
        self.handle = None

    def set_offset(self, offset):
        self._offset = offset
        if self.size == -1:
            file_size = os.path.getsize(self.filename) 
            file_size -= offset
            self.size = float(file_size) / self.itemsize
            assert(int(self.size) == self.size)

    def build_memmap(self):
        assert(self.size != -1)
        self.handle = np.memmap(self.filename, dtype=self.dtype,
                        mode='r', shape=self.size, offset=self._offset)
        for k in self.dtype.names:
            self.data[k] = self.handle[k]

class SDFRead(dict):

    """docstring for SDFRead"""

    _eof = 'SDF-EOH'

    def __init__(self, filename):
        self.filename = filename
        self.parameters = {}
        self.structs = []
        self.comments = []
        self.parse_header()
        self.set_offsets()
        self.load_memmaps()

    def parse_header(self):
        """docstring for parse_header"""
        # Pre-process
        ascfile = open(self.filename, 'r')
        while True:
            l = ascfile.readline()
            if self._eof in l: break

            self.parse_line(l, ascfile)

        hoff = ascfile.tell()
        ascfile.close()
        self.parameters['header_offset'] = hoff

    def parse_line(self, line, ascfile):
        """Parse a line of sdf"""


        if 'struct' in line:
            self.parse_struct(line, ascfile)
            return

        if "#" in line:
            self.comments.append(line)
            return

        spl = lstrip(line.split("="))
        vtype, vname = lstrip(spl[0].split())
        vname = vname.strip("[]")
        vval = spl[-1].strip(";")
        if vtype == 'parameter':
            self.parameters[vname] = vval
            return
        elif vtype == "char":
            vtype = "str"

        try:
            vval = eval("np."+vtype+"(%s)" % vval)
        except AttributeError:
            vval = eval("np."+_types[vtype]+"(%s)" % vval)

        self.parameters[vname] = vval

    def parse_struct(self, line, ascfile):
        assert 'struct' in line

        str_types = []
        comments = []
        str_lines = []
        l = ascfile.readline()
        while "}" not in l:
            vtype, vnames = get_struct_vars(l)
            for v in vnames:
                str_types.append((v, vtype))
            l = ascfile.readline()
        num = l.strip("}[]")
        num = num.strip("\;\\\n]")
        if len(num) == 0:
            # We need to compute the number of records.  The DataStruct will
            # handle this.
            num = '-1'
        num = int(num)
        struct = DataStruct(str_types, num, self.filename)
        self.structs.append(struct)
        return

    def set_offsets(self):
        running_off = self.parameters['header_offset']
        for struct in self.structs:
            struct.set_offset(running_off)
            running_off += struct.size * struct.itemsize
        return

    def load_memmaps(self):
        for struct in self.structs:
            struct.build_memmap()
            self.update(struct.data)



