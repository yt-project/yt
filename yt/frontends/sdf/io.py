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
        return self.ds.sdf_container

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
                self._handle['x'], self._handle['y'], self._handle['z'])

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
                x = self._handle['x']
                y = self._handle['y']
                z = self._handle['z']
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None: continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.ds.parameters["particle_mass"]
                    else:
                        data = self._handle[field][mask]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        x, y, z = (self._handle[ax] for ax in 'xyz')
        pcount = x.size
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        while ind < pcount:
            npart = min(CHUNKSIZE, pcount - ind)
            pos = np.empty((npart, 3), dtype=x.dtype)
            pos[:,0] = x[ind:ind+npart]
            pos[:,1] = y[ind:ind+npart]
            pos[:,2] = z[ind:ind+npart]
            if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
               np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       self.ds.domain_left_edge,
                                       self.ds.domain_right_edge)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+npart] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge)
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

    def __init__(self, filename, header=None):
        self.filename = filename
        if header is None:
            header = filename
        self.header = header
        self.parameters = {}
        self.structs = []
        self.comments = []
        self.parse_header()
        self.set_offsets()
        self.load_memmaps()

    def parse_header(self):
        """docstring for parse_header"""
        # Pre-process
        ascfile = open(self.header, 'r')
        while True:
            l = ascfile.readline()
            if self._eof in l: break

            self.parse_line(l, ascfile)

        hoff = ascfile.tell()
        ascfile.close()
        if self.header != self.filename:
            hoff = 0
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


class SDFIndex(object):

    """docstring for SDFIndex

    This provides an index mechanism into the full SDF Dataset.

    Most useful class methods:
        get_cell_data(level, cell_iarr, fields)
        iter_bbox_data(left, right, fields)
        iter_bbox_data(left, right, fields)

    """
    def __init__(self, sdfdata, indexdata, level=9):
        super(SDFIndex, self).__init__()
        self.sdfdata = sdfdata
        self.indexdata = indexdata
        self.level = level
        self.rmin = None
        self.rmax = None
        self.domain_width = None
        self.domain_buffer = 0
        self.domain_dims = 0
        self.domain_active_dims = 0
        self.masks = {
            "p" : int("011"*level, 2),
            "t" : int("101"*level, 2),
            "r" : int("110"*level, 2),
            "z" : int("011"*level, 2),
            "y" : int("101"*level, 2),
            "x" : int("110"*level, 2),
            2 : int("011"*level, 2),
            1 : int("101"*level, 2),
            0 : int("110"*level, 2),
        }
        self.dim_slices = {
            "p" : slice(0, None, 3),
            "t" : slice(1, None, 3),
            "r" : slice(2, None, 3),
            "z" : slice(0, None, 3),
            "y" : slice(1, None, 3),
            "x" : slice(2, None, 3),
            2 : slice(0, None, 3),
            1 : slice(1, None, 3),
            0 : slice(2, None, 3),
        }
        self.set_bounds()

    def set_bounds(self):
        r_0 = self.sdfdata.parameters['R0']
        DW = 2.0 * r_0

        self.rmin = np.zeros(3)
        self.rmax = np.zeros(3)
        sorted_rtp = self.sdfdata.parameters.get("sorted_rtp", False)
        if sorted_rtp:
            self.rmin[:] = [0.0, 0.0, -np.pi]
            self.rmax[:] = [r_0*1.01, 2*np.pi, np.pi]
        else:
            self.rmin[0] -= self.sdfdata.parameters.get('Rx', 0.0)
            self.rmin[1] -= self.sdfdata.parameters.get('Ry', 0.0)
            self.rmin[2] -= self.sdfdata.parameters.get('Rz', 0.0)
            self.rmax[0] += self.sdfdata.parameters.get('Rx', r_0)
            self.rmax[1] += self.sdfdata.parameters.get('Ry', r_0)
            self.rmax[2] += self.sdfdata.parameters.get('Rz', r_0)

        #/* expand root for non-power-of-two */
        expand_root = 0.0
        ic_Nmesh = self.sdfdata.parameters.get('ic_Nmesh',0)
        if ic_Nmesh != 0:
            f2 = 1<<int(np.log2(ic_Nmesh-1)+1)
            if (f2 != ic_Nmesh):
                expand_root = 1.0*f2/ic_Nmesh - 1.0;
            print 'Expanding: ', f2, ic_Nmesh, expand_root
        self.rmin *= 1.0 + expand_root
        self.rmax *= 1.0 + expand_root
        self.domain_width = self.rmax - self.rmin
        self.domain_dims = 1 << self.level
        self.domain_buffer = (self.domain_dims - int(self.domain_dims/(1.0 + expand_root)))/2
        self.domain_active_dims = self.domain_dims - 2*self.domain_buffer
        print 'Domain stuff:', self.domain_width, self.domain_dims, self.domain_active_dims

    def get_key(self, iarr, level=None):
        if level is None:
            level = self.level
        i1, i2, i3 = iarr
        rep1 = np.binary_repr(i1, width=self.level)
        rep2 = np.binary_repr(i2, width=self.level)
        rep3 = np.binary_repr(i3, width=self.level)
        inter = np.zeros(self.level*3, dtype='c')
        inter[self.dim_slices[0]] = rep1
        inter[self.dim_slices[1]] = rep2
        inter[self.dim_slices[2]] = rep3
        return int(inter.tostring(), 2)

    def get_key_ijk(self, i1, i2, i3, level=None):
        return self.get_key(np.array([i1, i2, i3]), level=level)

    def get_slice_key(self, ind, dim='r'):
        slb = np.binary_repr(ind, width=self.level)
        expanded = np.array([0]*self.level*3, dtype='c')
        expanded[self.dim_slices[dim]] = slb
        return int(expanded.tostring(), 2)

    def get_slice_chunks(self, slice_dim, slice_index):
        sl_key = self.get_slice_key(slice_index, dim=slice_dim)
        mask = (self.indexdata['index'] & ~self.masks[slice_dim]) == sl_key
        offsets = self.indexdata['base'][mask]
        lengths = self.indexdata['len'][mask]
        return mask, offsets, lengths

    def get_ibbox_slow(self, ileft, iright):
        """
        Given left and right indicies, return a mask and
        set of offsets+lengths into the sdf data.
        """
        mask = np.zeros(self.indexdata['index'].shape, dtype='bool')
        ileft = np.array(ileft)
        iright = np.array(iright)
        for i in range(3):
            left_key = self.get_slice_key(ileft[i], dim=i)
            right_key= self.get_slice_key(iright[i], dim=i)
            dim_inds = (self.indexdata['index'] & ~self.masks[i])
            mask *= (dim_inds >= left_key) * (dim_inds <= right_key)
            del dim_inds

        offsets = self.indexdata['base'][mask]
        lengths = self.indexdata['len'][mask]
        return mask, offsets, lengths

    def get_ibbox(self, ileft, iright):
        """
        Given left and right indicies, return a mask and
        set of offsets+lengths into the sdf data.
        """
        mask = np.zeros(self.indexdata['index'].shape, dtype='bool')

        print 'Getting data from ileft to iright:',  ileft, iright

        X, Y, Z = np.mgrid[ileft[0]:iright[0]+1,
                           ileft[1]:iright[1]+1,
                           ileft[2]:iright[2]+1]

        X = X.ravel()
        Y = Y.ravel()
        Z = Z.ravel()
        # Correct For periodicity
        X[X < self.domain_buffer] += self.domain_active_dims
        X[X >= self.domain_dims -  self.domain_buffer] -= self.domain_active_dims
        Y[Y < self.domain_buffer] += self.domain_active_dims
        Y[Y >= self.domain_dims -  self.domain_buffer] -= self.domain_active_dims
        Z[Z < self.domain_buffer] += self.domain_active_dims
        Z[Z >= self.domain_dims -  self.domain_buffer] -= self.domain_active_dims

        print 'periodic:',  X.min(), X.max(), Y.min(), Y.max(), Z.min(), Z.max()

        indices = np.array([self.get_key_ijk(x, y, z) for x, y, z in zip(X, Y, Z)])
        indices = indices[indices < self.indexdata['index'].shape[0]]
        return indices

    def get_bbox(self, left, right):
        """
        Given left and right indicies, return a mask and
        set of offsets+lengths into the sdf data.
        """
        ileft = np.floor((left - self.rmin) / self.domain_width *  self.domain_dims)
        iright = np.floor((right - self.rmin) / self.domain_width * self.domain_dims)

        return self.get_ibbox(ileft, iright)

    def get_data(self, chunk, fields):
        data = {}
        for field in fields:
            data[field] = self.sdfdata[field][chunk]
        return data

    def iter_data(self, inds, fields):
        num_inds = len(inds)
        num_reads = 0
        print 'Reading %i chunks' % num_inds
        i = 0
        while (i < num_inds):
            ind = inds[i]
            base = self.indexdata['base'][ind]
            length = self.indexdata['len'][ind]
            # Concatenate aligned reads
            nexti = i+1
            combined = 0
            while nexti < len(inds):
                nextind = inds[nexti]
                #        print 'b: %i l: %i end: %i  next: %i' % ( base, length, base + length, self.indexdata['base'][nextind] )
                if base + length == self.indexdata['base'][nextind]:
                    length += self.indexdata['len'][nextind]
                    i += 1
                    nexti += 1
                    combined += 1
                else:
                    break

            chunk = slice(base, base+length)
            print 'Reading chunk %i of length %i after catting %i' % (i, length, combined)
            num_reads += 1
            data = self.get_data(chunk, fields)
            yield data
            del data
            i += 1
        print 'Read %i chunks, batched into %i reads' % (num_inds, num_reads)

    def iter_bbox_data(self, left, right, fields):
        print 'Loading region from ', left, 'to', right
        inds = self.get_bbox(left, right)
        return self.iter_data(inds, fields)

    def iter_ibbox_data(self, left, right, fields):
        print 'Loading region from ', left, 'to', right
        inds = self.get_ibbox(left, right)
        return self.iter_data(inds, fields)

    def get_contiguous_chunk(self, left_key, right_key, fields):
        max_key = self.indexdata['index'][-1]
        if left_key > max_key:
            raise RuntimeError("Left key is too large. Key: %i Max Key: %i" % (left_key, max_key))
        base = self.indexdata['base'][left_key]
        right_key = min(right_key, self.indexdata['index'][-1])
        length = self.indexdata['base'][right_key] + \
            self.indexdata['len'][right_key] - base
        print 'Getting contiguous chunk of size %i starting at %i' % (length, base)
        return self.get_data(slice(base, base + length), fields)

    def iter_slice_data(self, slice_dim, slice_index, fields):
        mask, offsets, lengths = self.get_slice_chunks(slice_dim, slice_index)
        for off, l in zip(offsets, lengths):
            data = {}
            chunk = slice(off, off+l)
            for field in fields:
                data[field] = self.sdfdata[field][chunk]
            yield data
            del data

    def get_key_bounds(self, level, cell_iarr):
        """
        Get index keys for index file supplied.

        level: int
            Requested level
        cell_iarr: array-like, length 3
            Requested cell from given level.

        Returns:
            lmax_lk, lmax_rk
        """
        shift = self.level-level
        level_buff = 0
        level_lk = self.get_key(cell_iarr + level_buff)
        level_rk = self.get_key(cell_iarr + level_buff) + 1
        lmax_lk = (level_lk << shift*3)
        lmax_rk = (((level_rk) << shift*3) -1)
        #print "Level ", level, np.binary_repr(level_lk, width=self.level*3), np.binary_repr(level_rk, width=self.level*3)
        #print "Level ", self.level, np.binary_repr(lmax_lk, width=self.level*3), np.binary_repr(lmax_rk, width=self.level*3)
        return lmax_lk, lmax_rk

    def get_cell_data(self, level, cell_iarr, fields):
        """
        Get data from requested cell

        This uses the raw cell index, and doesn't account for periodicity or
        an expanded domain (non-power of 2).

        level: int
            Requested level
        cell_iarr: array-like, length 3
            Requested cell from given level.         fields: list
            Requested fields

        Returns:
            cell_data: dict
                Dictionary of field_name, field_data
        """
        cell_iarr = np.array(cell_iarr)
        lk, rk =self.get_key_bounds(level, cell_iarr)
        return self.get_contiguous_chunk(lk, rk, fields)
