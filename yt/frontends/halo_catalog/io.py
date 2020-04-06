import numpy as np

from yt.funcs import mylog, parse_h5_attr
from yt.units.yt_array import uvstack
from yt.utilities.exceptions import YTDomainOverflow
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.lib.geometry_utils import compute_morton
from yt.utilities.on_demand_imports import _h5py as h5py


class IOHandlerHaloCatalogHDF5(BaseIOHandler):
    _dataset_type = "halocatalog_hdf5"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == "halos"
        ptype = "halos"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        pn = "particle_position_%s"
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            with h5py.File(data_file.filename, mode="r") as f:
                units = parse_h5_attr(f[pn % "x"], "units")
                pos = data_file._get_particle_positions(ptype, f=f)
                x, y, z = (self.ds.arr(pos[:, i], units) for i in range(3))
                yield "halos", (x, y, z)

    def _yield_coordinates(self, data_file):
        pn = "particle_position_%s"
        with h5py.File(data_file.filename, "r") as f:
            units = parse_h5_attr(f[pn % "x"], "units")
            x, y, z = (
                self.ds.arr(f[pn % ax].value.astype("float64"), units) for ax in "xyz"
            )
            pos = uvstack([x, y, z]).T
            pos.convert_to_units("code_length")
            yield "halos", pos

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == "halos"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        pn = "particle_position_%s"
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            with h5py.File(data_file.filename, mode="r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    units = parse_h5_attr(f[pn % "x"], "units")
                    pos = data_file._get_particle_positions(ptype, f=f)
                    x, y, z = (self.ds.arr(pos[:, i], units) for i in range(3))
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None:
                        continue
                    for field in field_list:
                        data = f[field][si:ei][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = data_file.header["num_halos"]
        morton = np.empty(pcount, dtype="uint64")
        mylog.debug(
            "Initializing index % 5i (% 7i particles)", data_file.file_id, pcount
        )
        ind = 0
        if pcount == 0:
            return None
        ptype = "halos"
        with h5py.File(data_file.filename, mode="r") as f:
            if not f.keys():
                return None
            units = parse_h5_attr(f["particle_position_x"], "units")
            pos = data_file._get_particle_positions(ptype, f=f)
            pos = data_file.ds.arr(pos, units).to("code_length")
            dle = self.ds.domain_left_edge.to("code_length")
            dre = self.ds.domain_right_edge.to("code_length")
            if np.any(pos.min(axis=0) < dle) or np.any(pos.max(axis=0) > dre):
                raise YTDomainOverflow(pos.min(axis=0), pos.max(axis=0), dle, dre)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind : ind + pos.shape[0]] = compute_morton(
                pos[:, 0], pos[:, 1], pos[:, 2], dle, dre
            )
        return morton

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        nhalos = data_file.header["num_halos"]
        if None not in (si, ei):
            nhalos = np.clip(nhalos - si, 0, ei - si)
        return {"halos": nhalos}

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            fields = [("halos", field) for field in f
                      if not isinstance(f[field], h5py.Group)]
            units = dict([(("halos", field),
                           parse_h5_attr(f[field], "units"))
                          for field in f])
        return fields, units

class IOHandlerHaloCatalogHaloHDF5(IOHandlerHaloCatalogHDF5):
    _dataset_type = "halo_catalog_halo_hdf5"

    def _read_particle_coords(self, chunks, ptf):
        pass

    def _read_particle_selection(self, dobj, fields):
        rv = {}
        ind = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)        # ON-DISK TO READ
        fsize = defaultdict(lambda: 0) # COUNT RV
        field_maps = defaultdict(list) # ptypes -> fields
        unions = self.ds.particle_unions
        # What we need is a mapping from particle types to return types
        for field in fields:
            ftype, fname = field
            fsize[field] = 0
            # We should add a check for p.fparticle_unions or something here
            if ftype in unions:
                for pt in unions[ftype]:
                    ptf[pt].append(fname)
                    field_maps[pt, fname].append(field)
            else:
                ptf[ftype].append(fname)
                field_maps[field].append(field)

        # Now we allocate
        psize = {dobj.ptype: dobj.particle_number}
        for field in fields:
            if field[0] in unions:
                for pt in unions[field[0]]:
                    fsize[field] += psize.get(pt, 0)
            else:
                fsize[field] += psize.get(field[0], 0)
        for field in fields:
            if field[1] in self._vector_fields:
                shape = (fsize[field], self._vector_fields[field[1]])
            elif field[1] in self._array_fields:
                shape = (fsize[field],)+self._array_fields[field[1]]
            elif field in self.ds.scalar_field_list:
                shape = (1,)
            else:
                shape = (fsize[field], )
            rv[field] = np.empty(shape, dtype="float64")
            ind[field] = 0
        # Now we read.
        for field_r, vals in self._read_particle_fields(dobj, ptf):
            # Note that we now need to check the mappings
            for field_f in field_maps[field_r]:
                my_ind = ind[field_f]
                rv[field_f][my_ind:my_ind + vals.shape[0],...] = vals
                ind[field_f] += vals.shape[0]
        # Now we need to truncate all our fields, since we allow for
        # over-estimating.
        for field_f in ind:
            rv[field_f] = rv[field_f][:ind[field_f]]
        return rv

    def _read_scalar_fields(self, dobj, scalar_fields):
        all_data = {}
        if not scalar_fields: return all_data
        pcount = 1
        with h5py.File(dobj.scalar_data_file.filename, "r") as f:
            for ptype, field_list in sorted(scalar_fields.items()):
                for field in field_list:
                    if field == "particle_identifier":
                        field_data = \
                          np.arange(dobj.scalar_data_file.total_particles[ptype]) + \
                          dobj.scalar_data_file.index_start[ptype]
                    elif field in f[ptype]:
                        field_data = f[ptype][field][()].astype("float64")
                    else:
                        fname = field[:field.rfind("_")]
                        field_data = f[ptype][fname][()].astype("float64")
                        my_div = field_data.size / pcount
                        if my_div > 1:
                            findex = int(field[field.rfind("_") + 1:])
                            field_data = field_data[:, findex]
                    data = np.array([field_data[dobj.scalar_index]])
                    all_data[(ptype, field)] = data
        return all_data

    def _read_member_fields(self, dobj, member_fields):
        all_data = defaultdict(lambda: np.empty(dobj.particle_number,
                                                dtype=np.float64))
        if not member_fields: return all_data
        field_start = 0
        for i, data_file in enumerate(dobj.field_data_files):
            start_index = dobj.field_data_start[i]
            end_index = dobj.field_data_end[i]
            pcount = end_index - start_index
            if pcount == 0: continue
            field_end = field_start + end_index - start_index
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(member_fields.items()):
                    for field in field_list:
                        field_data = all_data[(ptype, field)]
                        if field in f["IDs"]:
                            my_data = \
                              f["IDs"][field][start_index:end_index].astype("float64")
                        else:
                            fname = field[:field.rfind("_")]
                            my_data = \
                              f["IDs"][fname][start_index:end_index].astype("float64")
                            my_div = my_data.size / pcount
                            if my_div > 1:
                                findex = int(field[field.rfind("_") + 1:])
                                my_data = my_data[:, findex]
                        field_data[field_start:field_end] = my_data
            field_start = field_end
        return all_data

    def _read_particle_fields(self, dobj, ptf):
        # separate member particle fields from scalar fields
        scalar_fields = defaultdict(list)
        member_fields = defaultdict(list)
        for ptype, field_list in sorted(ptf.items()):
            for field in field_list:
                if (ptype, field) in self.ds.scalar_field_list:
                    scalar_fields[ptype].append(field)
                else:
                    member_fields[ptype].append(field)

        all_data = self._read_scalar_fields(dobj, scalar_fields)
        all_data.update(self._read_member_fields(dobj, member_fields))

        for field, field_data in all_data.items():
            yield field, field_data

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            scalar_fields = [("halos", field) for field in f
                        if not isinstance(f[field], h5py.Group)]
            units = dict([(("halos", field),
                           parse_h5_attr(f[field], "units"))
                          for field in f])
            if 'particles' in f:
                id_fields = [('halos', field) for field in f['particles']]
            else:
                id_fields = []

        return scalar_fields+id_fields, scalar_fields, id_fields, units
