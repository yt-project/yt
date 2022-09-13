import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseParticleIOHandler
from yt.utilities.on_demand_imports import _requests as requests


class IOHandlerHTTPStream(BaseParticleIOHandler):
    _dataset_type = "http_particle_stream"
    _vector_fields = {"Coordinates": 3, "Velocity": 3, "Velocities": 3}

    def __init__(self, ds):
        self._url = ds.base_url
        # This should eventually manage the IO and cache it
        self.total_bytes = 0
        super().__init__(ds)

    def _open_stream(self, data_file, field):
        # This does not actually stream yet!
        ftype, fname = field
        s = f"{self._url}/{data_file.file_id}/{ftype}/{fname}"
        mylog.info("Loading URL %s", s)
        resp = requests.get(s)
        if resp.status_code != 200:
            raise RuntimeError
        self.total_bytes += len(resp.content)
        return resp.content

    def _identify_fields(self, data_file):
        f = []
        for ftype, fname in self.ds.parameters["field_list"]:
            f.append((str(ftype), str(fname)))
        return f, {}

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            for ptype in ptf:
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0] / 3.0, 3)
                yield ptype, (c[:, 0], c[:, 1], c[:, 2]), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            for ptype, field_list in sorted(ptf.items()):
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0] / 3.0, 3)
                mask = selector.select_points(c[:, 0], c[:, 1], c[:, 2], 0.0)
                del c
                if mask is None:
                    continue
                for field in field_list:
                    s = self._open_stream(data_file, (ptype, field))
                    c = np.frombuffer(s, dtype="float64")
                    if field in self._vector_fields:
                        c.shape = (c.shape[0] / 3.0, 3)
                    data = c[mask, ...]
                    yield (ptype, field), data

    def _count_particles(self, data_file):
        return self.ds.parameters["particle_count"][data_file.file_id]
