"""
OWLSSubfind data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from yt.utilities.exceptions import YTDomainOverflow
from yt.funcs import mylog

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.lib.geometry_utils import compute_morton

class IOHandlerOWLSSubfindHDF5(BaseIOHandler):
    _dataset_type = "subfind_hdf5"

    def __init__(self, ds):
        super(IOHandlerOWLSSubfindHDF5, self).__init__(ds)
        self.offset_fields = set([])

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda f: f.filename):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    coords = f[ptype]["CenterOfMass"].value.astype("float64")
                    coords = np.resize(coords, (pcount, 3))
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    yield ptype, (x, y, z)

    def _read_offset_particle_field(self, field, data_file, fh):
        field_data = np.empty(data_file.total_particles["FOF"], dtype="float64")
        fofindex = np.arange(data_file.total_particles["FOF"]) + data_file.index_start["FOF"]
        for offset_file in data_file.offset_files:
            if fh.filename == offset_file.filename:
                ofh = fh
            else:
                ofh = h5py.File(offset_file.filename, "r")
            subindex = np.arange(offset_file.total_offset) + offset_file.offset_start
            substart = max(fofindex[0] - subindex[0], 0)
            subend = min(fofindex[-1] - subindex[0], subindex.size - 1)
            fofstart = substart + subindex[0] - fofindex[0]
            fofend = subend + subindex[0] - fofindex[0]
            field_data[fofstart:fofend + 1] = ofh["SUBFIND"][field][substart:subend + 1]
        return field_data
                    
    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda f: f.filename):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    if pcount == 0: continue
                    coords = f[ptype]["CenterOfMass"].value.astype("float64")
                    coords = np.resize(coords, (pcount, 3))
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        if field in self.offset_fields:
                            field_data = \
                              self._read_offset_particle_field(field, data_file, f)
                        else:
                            if field == "particle_identifier":
                                field_data = \
                                  np.arange(data_file.total_particles[ptype]) + \
                                  data_file.index_start[ptype]
                            elif field in f[ptype]:
                                field_data = f[ptype][field].value.astype("float64")
                            else:
                                fname = field[:field.rfind("_")]
                                field_data = f[ptype][fname].value.astype("float64")
                                my_div = field_data.size / pcount
                                if my_div > 1:
                                    field_data = np.resize(field_data, (pcount, my_div))
                                    findex = int(field[field.rfind("_") + 1:])
                                    field_data = field_data[:, findex]
                        data = field_data[mask]
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = sum(data_file.total_particles.values())
        morton = np.empty(pcount, dtype='uint64')
        if pcount == 0: return morton
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        with h5py.File(data_file.filename, "r") as f:
            if not f.keys(): return None
            dx = np.finfo(f["FOF"]['CenterOfMass'].dtype).eps
            dx = 2.0*self.ds.quan(dx, "code_length")

            for ptype in data_file.ds.particle_types_raw:
                if data_file.total_particles[ptype] == 0: continue
                pos = f[ptype]["CenterOfMass"].value.astype("float64")
                pos = np.resize(pos, (data_file.total_particles[ptype], 3))
                pos = data_file.ds.arr(pos, "code_length")
                
                # These are 32 bit numbers, so we give a little lee-way.
                # Otherwise, for big sets of particles, we often will bump into the
                # domain edges.  This helps alleviate that.
                np.clip(pos, self.ds.domain_left_edge + dx,
                             self.ds.domain_right_edge - dx, pos)
                if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
                   np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                    raise YTDomainOverflow(pos.min(axis=0),
                                           pos.max(axis=0),
                                           self.ds.domain_left_edge,
                                           self.ds.domain_right_edge)
                regions.add_data_file(pos, data_file.file_id)
                morton[ind:ind+pos.shape[0]] = compute_morton(
                    pos[:,0], pos[:,1], pos[:,2],
                    data_file.ds.domain_left_edge,
                    data_file.ds.domain_right_edge)
                ind += pos.shape[0]
        return morton

    def _count_particles(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            pcount = {"FOF": f["FOF"].attrs["Number_of_groups"]}
            if "SUBFIND" in f:
                # We need this to figure out where the offset fields are stored.
                data_file.total_offset = f["SUBFIND"].attrs["Number_of_groups"]
                pcount["SUBFIND"] = f["FOF"].attrs["Number_of_subgroups"]
            else:
                data_file.total_offset = 0
                pcount["SUBFIND"] = 0
            return pcount

    def _identify_fields(self, data_file):
        fields = []
        pcount = data_file.total_particles
        if sum(pcount.values()) == 0: return fields, {}
        with h5py.File(data_file.filename, "r") as f:
            for ptype in self.ds.particle_types_raw:
                if data_file.total_particles[ptype] == 0: continue
                fields.append((ptype, "particle_identifier"))
                my_fields, my_offset_fields = \
                  subfind_field_list(f[ptype], ptype, data_file.total_particles)
                fields.extend(my_fields)
                self.offset_fields = self.offset_fields.union(set(my_offset_fields))
        return fields, {}

def subfind_field_list(fh, ptype, pcount):
    fields = []
    offset_fields = []
    for field in fh.keys():
        if "PartType" in field:
            # These are halo member particles
            continue
        elif isinstance(fh[field], h5py.Group):
            my_fields, my_offset_fields = \
              subfind_field_list(fh[field], ptype, pcount)
            fields.extend(my_fields)
            my_offset_fields.extend(offset_fields)
        else:
            if not fh[field].size % pcount[ptype]:
                my_div = fh[field].size / pcount[ptype]
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1:]
                if my_div > 1:
                    for i in range(int(my_div)):
                        fields.append((ptype, "%s_%d" % (fname, i)))
                else:
                    fields.append((ptype, fname))
            elif ptype == "SUBFIND" and \
              not fh[field].size % fh["/SUBFIND"].attrs["Number_of_groups"]:
                # These are actually FOF fields, but they were written after 
                # a load balancing step moved halos around and thus they do not
                # correspond to the halos stored in the FOF group.
                my_div = fh[field].size / fh["/SUBFIND"].attrs["Number_of_groups"]
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1:]
                if my_div > 1:
                    for i in range(int(my_div)):
                        fields.append(("FOF", "%s_%d" % (fname, i)))
                else:
                    fields.append(("FOF", fname))
                offset_fields.append(fname)
            else:
                mylog.warn("Cannot add field (%s, %s) with size %d." % \
                           (ptype, fh[field].name, fh[field].size))
                continue
    return fields, offset_fields
