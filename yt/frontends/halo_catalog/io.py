"""
HaloCatalog data-file handling function




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


class IOHandlerHaloCatalogHDF5(BaseIOHandler):
    _dataset_type = "halocatalog_hdf5"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "halos")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        pn = "particle_position_%s"
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                units = f[pn % "x"].attrs["units"]
                x, y, z = \
                  (self.ds.arr(f[pn % ax].value.astype("float64"), units)
                   for ax in "xyz")
                yield "halos", (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "halos")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        pn = "particle_position_%s"
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    units = f[pn % "x"].attrs["units"]
                    x, y, z = \
                      (self.ds.arr(f[pn % ax].value.astype("float64"), units)
                       for ax in "xyz")
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f[field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = data_file.header["num_halos"]
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        if pcount == 0: return None
        with h5py.File(data_file.filename, "r") as f:
            if not f.keys(): return None
            pos = np.empty((pcount, 3), dtype="float64")
            units = f["particle_position_x"].attrs["units"]
            dx = np.finfo(f['particle_position_x'].dtype).eps
            dx = 2.0 * self.ds.quan(dx, units).to("code_length")
            pos[:,0] = f["particle_position_x"].value
            pos[:,1] = f["particle_position_y"].value
            pos[:,2] = f["particle_position_z"].value
            pos = data_file.ds.arr(pos, units). to("code_length")
            dle = self.ds.domain_left_edge.to("code_length")
            dre = self.ds.domain_right_edge.to("code_length")
            # These are 32 bit numbers, so we give a little lee-way.
            # Otherwise, for big sets of particles, we often will bump into the
            # domain edges.  This helps alleviate that.
            np.clip(pos, dle + dx, dre - dx, pos)
            if np.any(pos.min(axis=0) < dle) or \
               np.any(pos.max(axis=0) > dre):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       dle, dre)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+pos.shape[0]] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2], dle, dre)
        return morton

    def _count_particles(self, data_file):
        return {'halos': data_file.header['num_halos']}

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            fields = [("halos", field) for field in f]
            units = dict([(("halos", field), 
                           f[field].attrs["units"]) for field in f])
        return fields, units
