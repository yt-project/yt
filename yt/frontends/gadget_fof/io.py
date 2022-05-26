from collections import defaultdict

import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


class IOHandlerGadgetFOFHDF5(BaseParticleIOHandler):
    _dataset_type = "gadget_fof_hdf5"

    def __init__(self, ds):
        super().__init__(ds)
        self.offset_fields = set()

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError(
            "IOHandlerGadgetFOFHDF5 _read_fluid_selection not implemented yet"
        )

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            with h5py.File(data_file.filename, mode="r") as f:
                for ptype in sorted(ptf):
                    coords = data_file._get_particle_positions(ptype, f=f)
                    if coords is None:
                        continue
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    yield ptype, (x, y, z), 0.0

    def _yield_coordinates(self, data_file):
        ptypes = self.ds.particle_types_raw
        with h5py.File(data_file.filename, mode="r") as f:
            for ptype in sorted(ptypes):
                pcount = data_file.total_particles[ptype]
                if pcount == 0:
                    continue
                coords = f[ptype][f"{ptype}Pos"][()].astype("float64")
                coords = np.resize(coords, (pcount, 3))
                yield ptype, coords

    def _read_offset_particle_field(self, field, data_file, fh):
        field_data = np.empty(data_file.total_particles["Group"], dtype="float64")
        fofindex = (
            np.arange(data_file.total_particles["Group"])
            + data_file.index_start["Group"]
        )
        for offset_file in data_file.offset_files:
            if fh.filename == offset_file.filename:
                ofh = fh
            else:
                ofh = h5py.File(offset_file.filename, mode="r")
            subindex = np.arange(offset_file.total_offset) + offset_file.offset_start
            substart = max(fofindex[0] - subindex[0], 0)
            subend = min(fofindex[-1] - subindex[0], subindex.size - 1)
            fofstart = substart + subindex[0] - fofindex[0]
            fofend = subend + subindex[0] - fofindex[0]
            field_data[fofstart : fofend + 1] = ofh["Subhalo"][field][
                substart : subend + 1
            ]
        return field_data

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            with h5py.File(data_file.filename, mode="r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    if pcount == 0:
                        continue
                    coords = data_file._get_particle_positions(ptype, f=f)
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None:
                        continue
                    for field in field_list:
                        if field in self.offset_fields:
                            field_data = self._read_offset_particle_field(
                                field, data_file, f
                            )
                        else:
                            if field == "particle_identifier":
                                field_data = (
                                    np.arange(data_file.total_particles[ptype])
                                    + data_file.index_start[ptype]
                                )
                            elif field in f[ptype]:
                                field_data = f[ptype][field][()].astype("float64")
                            else:
                                fname = field[: field.rfind("_")]
                                field_data = f[ptype][fname][()].astype("float64")
                                my_div = field_data.size / pcount
                                if my_div > 1:
                                    findex = int(field[field.rfind("_") + 1 :])
                                    field_data = field_data[:, findex]
                        data = field_data[si:ei][mask]
                        yield (ptype, field), data

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        pcount = {
            "Group": data_file.header["Ngroups_ThisFile"],
            "Subhalo": data_file.header["Nsubgroups_ThisFile"],
        }
        if None not in (si, ei):
            for ptype in pcount:
                pcount[ptype] = np.clip(pcount[ptype] - si, 0, ei - si)
        return pcount

    def _identify_fields(self, data_file):
        fields = []
        pcount = data_file.total_particles
        if sum(pcount.values()) == 0:
            return fields, {}
        with h5py.File(data_file.filename, mode="r") as f:
            for ptype in self.ds.particle_types_raw:
                if data_file.total_particles[ptype] == 0:
                    continue
                fields.append((ptype, "particle_identifier"))
                my_fields, my_offset_fields = subfind_field_list(
                    f[ptype], ptype, data_file.total_particles
                )
                fields.extend(my_fields)
                self.offset_fields = self.offset_fields.union(set(my_offset_fields))
        return fields, {}


class IOHandlerGadgetFOFHaloHDF5(IOHandlerGadgetFOFHDF5):
    _dataset_type = "gadget_fof_halo_hdf5"

    def _read_particle_coords(self, chunks, ptf):
        pass

    def _read_particle_selection(self, dobj, fields):
        rv = {}
        ind = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(list)  # ON-DISK TO READ
        fsize = defaultdict(lambda: 0)  # COUNT RV
        field_maps = defaultdict(list)  # ptypes -> fields
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
                shape = (fsize[field],) + self._array_fields[field[1]]
            elif field in self.ds.scalar_field_list:
                shape = (1,)
            else:
                shape = (fsize[field],)
            rv[field] = np.empty(shape, dtype="float64")
            ind[field] = 0
        # Now we read.
        for field_r, vals in self._read_particle_fields(dobj, ptf):
            # Note that we now need to check the mappings
            for field_f in field_maps[field_r]:
                my_ind = ind[field_f]
                rv[field_f][my_ind : my_ind + vals.shape[0], ...] = vals
                ind[field_f] += vals.shape[0]
        # Now we need to truncate all our fields, since we allow for
        # over-estimating.
        for field_f in ind:
            rv[field_f] = rv[field_f][: ind[field_f]]
        return rv

    def _read_scalar_fields(self, dobj, scalar_fields):
        all_data = {}
        if not scalar_fields:
            return all_data
        pcount = 1
        with h5py.File(dobj.scalar_data_file.filename, mode="r") as f:
            for ptype, field_list in sorted(scalar_fields.items()):
                for field in field_list:
                    if field == "particle_identifier":
                        field_data = (
                            np.arange(dobj.scalar_data_file.total_particles[ptype])
                            + dobj.scalar_data_file.index_start[ptype]
                        )
                    elif field in f[ptype]:
                        field_data = f[ptype][field][()].astype("float64")
                    else:
                        fname = field[: field.rfind("_")]
                        field_data = f[ptype][fname][()].astype("float64")
                        my_div = field_data.size / pcount
                        if my_div > 1:
                            findex = int(field[field.rfind("_") + 1 :])
                            field_data = field_data[:, findex]
                    data = np.array([field_data[dobj.scalar_index]])
                    all_data[(ptype, field)] = data
        return all_data

    def _read_member_fields(self, dobj, member_fields):
        all_data = defaultdict(lambda: np.empty(dobj.particle_number, dtype=np.float64))
        if not member_fields:
            return all_data
        field_start = 0
        for i, data_file in enumerate(dobj.field_data_files):
            start_index = dobj.field_data_start[i]
            end_index = dobj.field_data_end[i]
            pcount = end_index - start_index
            if pcount == 0:
                continue
            field_end = field_start + end_index - start_index
            with h5py.File(data_file.filename, mode="r") as f:
                for ptype, field_list in sorted(member_fields.items()):
                    for field in field_list:
                        field_data = all_data[(ptype, field)]
                        if field in f["IDs"]:
                            my_data = f["IDs"][field][start_index:end_index].astype(
                                "float64"
                            )
                        else:
                            fname = field[: field.rfind("_")]
                            my_data = f["IDs"][fname][start_index:end_index].astype(
                                "float64"
                            )
                            my_div = my_data.size / pcount
                            if my_div > 1:
                                findex = int(field[field.rfind("_") + 1 :])
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
        fields = []
        scalar_fields = []
        id_fields = {}
        with h5py.File(data_file.filename, mode="r") as f:
            for ptype in self.ds.particle_types_raw:
                fields.append((ptype, "particle_identifier"))
                scalar_fields.append((ptype, "particle_identifier"))
                my_fields, my_offset_fields = subfind_field_list(
                    f[ptype], ptype, data_file.total_particles
                )
                fields.extend(my_fields)
                scalar_fields.extend(my_fields)

                if "IDs" not in f:
                    continue
                id_fields = [(ptype, field) for field in f["IDs"]]
                fields.extend(id_fields)
        return fields, scalar_fields, id_fields, {}


def subfind_field_list(fh, ptype, pcount):
    fields = []
    offset_fields = []
    for field in fh.keys():
        if isinstance(fh[field], h5py.Group):
            my_fields, my_offset_fields = subfind_field_list(fh[field], ptype, pcount)
            fields.extend(my_fields)
            my_offset_fields.extend(offset_fields)
        else:
            if not fh[field].size % pcount[ptype]:
                my_div = fh[field].size / pcount[ptype]
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1 :]
                if my_div > 1:
                    for i in range(int(my_div)):
                        fields.append((ptype, "%s_%d" % (fname, i)))
                else:
                    fields.append((ptype, fname))
            elif (
                ptype == "Subhalo"
                and not fh[field].size % fh["/Subhalo"].attrs["Number_of_groups"]
            ):
                # These are actually Group fields, but they were written after
                # a load balancing step moved halos around and thus they do not
                # correspond to the halos stored in the Group group.
                my_div = fh[field].size / fh["/Subhalo"].attrs["Number_of_groups"]
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1 :]
                if my_div > 1:
                    for i in range(int(my_div)):
                        fields.append(("Group", "%s_%d" % (fname, i)))
                else:
                    fields.append(("Group", fname))
                offset_fields.append(fname)
            else:
                mylog.warning(
                    "Cannot add field (%s, %s) with size %d.",
                    ptype,
                    fh[field].name,
                    fh[field].size,
                )
                continue
    return fields, offset_fields
