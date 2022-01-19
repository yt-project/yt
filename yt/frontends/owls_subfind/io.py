import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py

_pos_names = ["CenterOfMass", "CentreOfMass"]


class IOHandlerOWLSSubfindHDF5(BaseParticleIOHandler):
    _dataset_type = "subfind_hdf5"
    _position_name = None

    def __init__(self, ds):
        super().__init__(ds)
        self.offset_fields = set()

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

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
                    pcount = data_file.total_particles[ptype]
                    coords = f[ptype][self._position_name][()].astype("float64")
                    coords = np.resize(coords, (pcount, 3))
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    yield ptype, (x, y, z), 0.0

    def _yield_coordinates(self, data_file):
        ptypes = self.ds.particle_types_raw
        with h5py.File(data_file.filename, mode="r") as f:
            for ptype in sorted(ptypes):
                pcount = data_file.total_particles[ptype]
                coords = f[ptype][self._position_name][()].astype("float64")
                coords = np.resize(coords, (pcount, 3))
                yield ptype, coords

    def _read_offset_particle_field(self, field, data_file, fh):
        field_data = np.empty(data_file.total_particles["FOF"], dtype="float64")
        fofindex = (
            np.arange(data_file.total_particles["FOF"]) + data_file.index_start["FOF"]
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
            field_data[fofstart : fofend + 1] = ofh["SUBFIND"][field][
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
            with h5py.File(data_file.filename, mode="r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = data_file.total_particles[ptype]
                    if pcount == 0:
                        continue
                    coords = f[ptype][self._position_name][()].astype("float64")
                    coords = np.resize(coords, (pcount, 3))
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
                                    field_data = np.resize(
                                        field_data, (int(pcount), int(my_div))
                                    )
                                    findex = int(field[field.rfind("_") + 1 :])
                                    field_data = field_data[:, findex]
                        data = field_data[mask]
                        yield (ptype, field), data

    def _count_particles(self, data_file):
        with h5py.File(data_file.filename, mode="r") as f:
            pcount = {"FOF": get_one_attr(f["FOF"], ["Number_of_groups", "Ngroups"])}
            if "SUBFIND" in f:
                # We need this to figure out where the offset fields are stored.
                data_file.total_offset = get_one_attr(
                    f["SUBFIND"], ["Number_of_groups", "Ngroups"]
                )
                pcount["SUBFIND"] = get_one_attr(
                    f["FOF"], ["Number_of_subgroups", "Nsubgroups"]
                )
            else:
                data_file.total_offset = 0
                pcount["SUBFIND"] = 0

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
                self._identify_position_name(f[ptype])
                fields.append((ptype, "particle_identifier"))
                my_fields, my_offset_fields = subfind_field_list(
                    f[ptype], ptype, data_file.total_particles
                )
                fields.extend(my_fields)
                self.offset_fields = self.offset_fields.union(set(my_offset_fields))
        return fields, {}

    def _identify_position_name(self, fh):
        if self._position_name is not None:
            return
        for pname in _pos_names:
            if pname in fh:
                self._position_name = pname
                return


def get_one_attr(fh, attrs, default=None, error=True):
    """
    Try getting from a list of attrs. Return the first one that exists.
    """
    for attr in attrs:
        if attr in fh.attrs:
            return fh.attrs[attr]
    if error:
        raise RuntimeError(
            f"Could not find any of these attributes: {attrs}. "
            f"Available attributes: {fh.attrs.keys()}"
        )
    return default


def subfind_field_list(fh, ptype, pcount):
    fields = []
    offset_fields = []
    for field in fh.keys():
        if "PartType" in field:
            # These are halo member particles
            continue
        elif isinstance(fh[field], h5py.Group):
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
                ptype == "SUBFIND"
                and not fh[field].size % fh["/SUBFIND"].attrs["Number_of_groups"]
            ):
                # These are actually FOF fields, but they were written after
                # a load balancing step moved halos around and thus they do not
                # correspond to the halos stored in the FOF group.
                my_div = fh[field].size / fh["/SUBFIND"].attrs["Number_of_groups"]
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1 :]
                if my_div > 1:
                    for i in range(int(my_div)):
                        fields.append(("FOF", "%s_%d" % (fname, i)))
                else:
                    fields.append(("FOF", fname))
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
