from collections import defaultdict

import numpy as np

from yt.frontends.gadget_fof.io import IOHandlerGadgetFOFHaloHDF5
from yt.funcs import parse_h5_attr
from yt.units.yt_array import uvstack  # type: ignore
from yt.utilities.io_handler import BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


class IOHandlerYTHaloCatalog(BaseParticleIOHandler):
    _dataset_type = "ythalocatalog"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set()
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
                yield "halos", (x, y, z), 0.0

    def _yield_coordinates(self, data_file):
        pn = "particle_position_%s"
        with h5py.File(data_file.filename, mode="r") as f:
            units = parse_h5_attr(f[pn % "x"], "units")
            x, y, z = (
                self.ds.arr(f[pn % ax][()].astype("float64"), units) for ax in "xyz"
            )
            pos = uvstack([x, y, z]).T
            pos.convert_to_units("code_length")
            yield "halos", pos

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set()
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

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        nhalos = data_file.header["num_halos"]
        if None not in (si, ei):
            nhalos = np.clip(nhalos - si, 0, ei - si)
        return {"halos": nhalos}

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, mode="r") as f:
            fields = [
                ("halos", field) for field in f if not isinstance(f[field], h5py.Group)
            ]
            units = {("halos", field): parse_h5_attr(f[field], "units") for field in f}
        return fields, units


class HaloDatasetIOHandler:
    """
    Base class for io handlers to load halo member particles.
    """

    def _read_particle_coords(self, chunks, ptf):
        pass

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

    # This will be refactored.
    _read_particle_selection = IOHandlerGadgetFOFHaloHDF5._read_particle_selection


# ignoring type in this mixing to circunvent this error from mypy
# Definition of "_read_particle_fields" in base class "HaloDatasetIOHandler"
# is incompatible with definition in base class "IOHandlerYTHaloCatalog"
#
# it may not be possible to refactor out of this situation without breaking downstream
class IOHandlerYTHalo(HaloDatasetIOHandler, IOHandlerYTHaloCatalog):  # type: ignore
    _dataset_type = "ythalo"

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, mode="r") as f:
            scalar_fields = [
                ("halos", field) for field in f if not isinstance(f[field], h5py.Group)
            ]
            units = {("halos", field): parse_h5_attr(f[field], "units") for field in f}
            if "particles" in f:
                id_fields = [("halos", field) for field in f["particles"]]
            else:
                id_fields = []

        return scalar_fields + id_fields, scalar_fields, id_fields, units

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
                        my_data = f["particles"][field][start_index:end_index].astype(
                            "float64"
                        )
                        field_data[field_start:field_end] = my_data
            field_start = field_end
        return all_data

    def _read_scalar_fields(self, dobj, scalar_fields):
        all_data = {}
        if not scalar_fields:
            return all_data
        with h5py.File(dobj.scalar_data_file.filename, mode="r") as f:
            for ptype, field_list in sorted(scalar_fields.items()):
                for field in field_list:
                    data = np.array([f[field][dobj.scalar_index]]).astype("float64")
                    all_data[(ptype, field)] = data
        return all_data
