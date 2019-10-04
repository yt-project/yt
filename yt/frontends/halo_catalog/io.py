import numpy as np

from yt.funcs import \
    mylog, \
    parse_h5_attr
from yt.units.yt_array import \
    uvstack
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.exceptions import YTDomainOverflow
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
        ptype = 'halos'
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        pn = "particle_position_%s"
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            with h5py.File(data_file.filename, "r") as f:
                units = parse_h5_attr(f[pn % "x"], "units")
                pos = data_file._get_particle_positions(ptype, f=f)
                x, y, z = (self.ds.arr(pos[:, i], units)
                           for i in range(3))
                yield "halos", (x, y, z)

    def _yield_coordinates(self, data_file):
        pn = "particle_position_%s"
        with h5py.File(data_file.filename, 'r') as f:
            units = parse_h5_attr(f[pn % "x"], "units")
            x, y, z = (self.ds.arr(f[pn % ax].value.astype("float64"), units)
                       for ax in 'xyz')
            pos = uvstack([x, y, z]).T
            pos.convert_to_units('code_length')
            yield 'halos', pos

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
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    units = parse_h5_attr(f[pn % "x"], "units")
                    pos = data_file._get_particle_positions(ptype, f=f)
                    x, y, z = (self.ds.arr(pos[:, i], units)
                               for i in range(3))
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f[field][si:ei][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = data_file.header["num_halos"]
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        if pcount == 0: return None
        ptype = 'halos'
        with h5py.File(data_file.filename, "r") as f:
            if not f.keys(): return None
            units = parse_h5_attr(f["particle_position_x"], "units")
            pos = data_file._get_particle_positions(ptype, f=f)
            pos = data_file.ds.arr(pos, units). to("code_length")
            dle = self.ds.domain_left_edge.to("code_length")
            dre = self.ds.domain_right_edge.to("code_length")
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
        si, ei = data_file.start, data_file.end
        nhalos = data_file.header['num_halos']
        if None not in (si, ei):
            nhalos = np.clip(nhalos - si, 0, ei - si)
        return {'halos': nhalos}

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            fields = [("halos", field) for field in f]
            units = dict([(("halos", field),
                           parse_h5_attr(f[field], "units"))
                          for field in f])
        return fields, units
