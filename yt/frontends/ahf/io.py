from operator import attrgetter

import numpy as np

from yt.utilities.io_handler import BaseParticleIOHandler


class IOHandlerAHFHalos(BaseParticleIOHandler):
    _particle_reader = False
    _dataset_type = "ahf"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z), hsml).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        for data_file in self._get_data_files(chunks, ptf):
            pos = data_file._get_particle_positions("halos")
            x, y, z = (pos[:, i] for i in range(3))
            yield "halos", (x, y, z), 0.0

    def _yield_coordinates(self, data_file):
        halos = data_file.read_data(usecols=["Xc", "Yc", "Zc"])
        x = halos["Xc"].astype("float64")
        y = halos["Yc"].astype("float64")
        z = halos["Zc"].astype("float64")
        yield "halos", np.asarray((x, y, z)).T

    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        for data_file in self._get_data_files(chunks, ptf):
            si, ei = data_file.start, data_file.end
            cols = []
            for field_list in ptf.values():
                cols.extend(field_list)
            cols = list(set(cols))
            halos = data_file.read_data(usecols=cols)
            pos = data_file._get_particle_positions("halos")
            x, y, z = (pos[:, i] for i in range(3))
            yield "halos", (x, y, z)
            mask = selector.select_points(x, y, z, 0.0)
            del x, y, z
            if mask is None:
                continue
            for ptype, field_list in sorted(ptf.items()):
                for field in field_list:
                    data = halos[field][si:ei][mask].astype("float64")
                    yield (ptype, field), data

    def _count_particles(self, data_file):
        halos = data_file.read_data(usecols=["ID"])
        nhalos = len(halos["ID"])
        si, ei = data_file.start, data_file.end
        if None not in (si, ei):
            nhalos = np.clip(nhalos - si, 0, ei - si)
        return {"halos": nhalos}

    def _identify_fields(self, data_file):
        fields = [("halos", f) for f in data_file.col_names]
        return fields, {}

    # Helper methods

    def _get_data_files(self, chunks, ptf):
        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == "halos"
        # Get data_files
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        data_files = sorted(data_files, key=attrgetter("filename"))
        yield from data_files
