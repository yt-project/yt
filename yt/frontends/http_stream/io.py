"""
HTTPStream data-file handling function




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import \
    get_requests, \
    mylog
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton

class IOHandlerHTTPStream(BaseIOHandler):
    _dataset_type = "http_particle_stream"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")

    def __init__(self, ds):
        if get_requests() is None:
            raise ImportError(
                "This functionality depends on the requests package")
        self._url = ds.base_url
        # This should eventually manage the IO and cache it
        self.total_bytes = 0
        super(IOHandlerHTTPStream, self).__init__(ds)

    def _open_stream(self, data_file, field):
        # This does not actually stream yet!
        ftype, fname = field
        s = "%s/%s/%s/%s" % (self._url,
            data_file.file_id, ftype, fname)
        mylog.info("Loading URL %s", s)
        requests = get_requests()
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
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            for ptype in ptf:
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0]/3.0, 3)
                yield ptype, (c[:,0], c[:,1], c[:,2])

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            for ptype, field_list in sorted(ptf.items()):
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0]/3.0, 3)
                mask = selector.select_points(
                            c[:,0], c[:,1], c[:,2], 0.0)
                del c
                if mask is None: continue
                for field in field_list:
                    s = self._open_stream(data_file, (ptype, field))
                    c = np.frombuffer(s, dtype="float64")
                    if field in self._vector_fields:
                        c.shape = (c.shape[0]/3.0, 3)
                    data = c[mask, ...]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        header = self.ds.parameters
        ptypes = header["particle_count"][data_file.file_id].keys()
        pcount = sum(header["particle_count"][data_file.file_id].values())
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        for ptype in ptypes:
            s = self._open_stream(data_file, (ptype, "Coordinates"))
            c = np.frombuffer(s, dtype="float64")
            c.shape = (c.shape[0]/3.0, 3)
            regions.add_data_file(c, data_file.file_id,
                                  data_file.ds.filter_bbox)
            morton[ind:ind+c.shape[0]] = compute_morton(
                c[:,0], c[:,1], c[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge,
                data_file.ds.filter_bbox)
            ind += c.shape[0]
        return morton

    def _count_particles(self, data_file):
        return self.ds.parameters["particle_count"][data_file.file_id]
