import os

import numpy as np

from yt.utilities.io_handler import BaseParticleIOHandler

from .definitions import halo_dts


class IOHandlerRockstarBinary(BaseParticleIOHandler):
    _dataset_type = "rockstar_binary"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._halo_dt = halo_dts[self.ds.parameters["format_revision"]]

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
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            pcount = data_file.header["num_halos"]
            if pcount == 0:
                continue
            with open(data_file.filename, "rb") as f:
                pos = data_file._get_particle_positions(ptype, f=f)
                yield "halos", (pos[:, i] for i in range(3)), 0.0

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
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            pcount = data_file.header["num_halos"]
            if pcount == 0:
                continue
            with open(data_file.filename, "rb") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pos = data_file._get_particle_positions(ptype, f=f)
                    x, y, z = (pos[:, i] for i in range(3))
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    f.seek(data_file._position_offset, os.SEEK_SET)
                    halos = np.fromfile(f, dtype=self._halo_dt, count=pcount)
                    if mask is None:
                        continue
                    for field in field_list:
                        data = halos[field][si:ei][mask].astype("float64")
                        yield (ptype, field), data

    def _yield_coordinates(self, data_file):
        # Just does halos
        pcount = data_file.header["num_halos"]
        with open(data_file.filename, "rb") as f:
            f.seek(data_file._position_offset, os.SEEK_SET)
            halos = np.fromfile(f, dtype=self._halo_dt, count=pcount)
            pos = np.empty((halos.size, 3), dtype="float64")
            pos[:, 0] = halos["particle_position_x"]
            pos[:, 1] = halos["particle_position_y"]
            pos[:, 2] = halos["particle_position_z"]
            yield "halos", pos

    def _count_particles(self, data_file):
        nhalos = data_file.header["num_halos"]
        si, ei = data_file.start, data_file.end
        if None not in (si, ei):
            nhalos = np.clip(nhalos - si, 0, ei - si)
        return {"halos": nhalos}

    def _identify_fields(self, data_file):
        fields = [("halos", f) for f in self._halo_dt.fields if "padding" not in f]
        return fields, {}
