"""
SDF data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.funcs import mylog
from yt.utilities.exceptions import YTDomainOverflow
from yt.utilities.lib.geometry_utils import compute_morton

CHUNKSIZE = 32**3

class IOHandlerSDF(BaseIOHandler):
    _dataset_type = "sdf_particles"

    @property
    def _handle(self):
        return self.ds.sdf_container

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in sorted(data_files):
            yield "dark_matter", (
                self._handle['x'], self._handle['y'], self._handle['z'])

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in sorted(data_files):
            for ptype, field_list in sorted(ptf.items()):
                x = self._handle['x']
                y = self._handle['y']
                z = self._handle['z']
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None: continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.ds.parameters["particle_mass"]
                    else:
                        data = self._handle[field][mask]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        x, y, z = (self._handle[ax] for ax in 'xyz')
        pcount = x.size

        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        while ind < pcount:
            npart = min(CHUNKSIZE, pcount - ind)
            pos = np.empty((npart, 3), dtype=x.dtype)
            pos[:,0] = x[ind:ind+npart]
            pos[:,1] = y[ind:ind+npart]
            pos[:,2] = z[ind:ind+npart]
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+npart] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge)
            ind += CHUNKSIZE
        return morton

    def _identify_fields(self, data_file):
        fields = [("dark_matter", v) for v in self._handle.keys()]
        fields.append(("dark_matter", "mass"))
        return fields, {}

    def _count_particles(self, data_file):
        pcount = self._handle['x'].size
        if (pcount > 1e9):
            mylog.warn("About to load %i particles into memory. " % (pcount) +
                       "You may want to consider a midx-enabled load")
        return {'dark_matter': pcount}


class IOHandlerHTTPSDF(IOHandlerSDF):
    _dataset_type = "http_sdf_particles"

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in data_files:
            pcount = self._handle['x'].size
            yield "dark_matter", (
                self._handle['x'][:pcount], self._handle['y'][:pcount], self._handle['z'][:pcount])

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        assert(len(ptf) == 1)
        assert(ptf.keys()[0] == "dark_matter")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert(len(data_files) == 1)
        for data_file in data_files:
            pcount = self._handle['x'].size
            for ptype, field_list in sorted(ptf.items()):
                x = self._handle['x'][:pcount]
                y = self._handle['y'][:pcount]
                z = self._handle['z'][:pcount]
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None: continue
                for field in field_list:
                    if field == "mass":
                        if self.ds.field_info._mass_field is None:
                            pm = 1.0
                            if 'particle_mass' in self.ds.parameters:
                                pm = self.ds.parameters['particle_mass']
                            else:
                                raise RuntimeError
                            data = pm * np.ones(mask.sum(), dtype="float64")
                        else:
                            data = self._handle[self.ds.field_info._mass_field][:][mask]
                    else:
                        data = self._handle[field][:][mask]
                    yield (ptype, field), data

    def _count_particles(self, data_file):
        return {'dark_matter': self._handle['x'].http_array.shape}


class IOHandlerSIndexSDF(IOHandlerSDF):
    _dataset_type = "midx_sdf_particles"


    def _read_particle_coords(self, chunks, ptf):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        for dd in self.ds.midx.iter_bbox_data(
            dle, dre,
            ['x','y','z']):
            yield "dark_matter", (
                dd['x'], dd['y'], dd['z'])

    def _read_particle_fields(self, chunks, ptf, selector):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        required_fields = []
        for ptype, field_list in sorted(ptf.items()):
            for field in field_list:
                if field == "mass": continue
                required_fields.append(field)

        for dd in self.ds.midx.iter_bbox_data(
            dle, dre,
            required_fields):

            for ptype, field_list in sorted(ptf.items()):
                x = dd['x']
                y = dd['y']
                z = dd['z']
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None: continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.ds.parameters["particle_mass"]
                    else:
                        data = dd[field][mask]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        pcount = 0
        for dd in self.ds.midx.iter_bbox_data(
            dle, dre,
            ['x']):
            pcount += dd['x'].size

        morton = np.empty(pcount, dtype='uint64')
        ind = 0

        chunk_id = 0
        for dd in self.ds.midx.iter_bbox_data(
            dle, dre,
            ['x','y','z']):
            npart = dd['x'].size
            pos = np.empty((npart, 3), dtype=dd['x'].dtype)
            pos[:,0] = dd['x']
            pos[:,1] = dd['y']
            pos[:,2] = dd['z']
            if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
               np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       self.ds.domain_left_edge,
                                       self.ds.domain_right_edge)
            regions.add_data_file(pos, chunk_id)
            morton[ind:ind+npart] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge)
            ind += npart
        return morton

    def _count_particles(self, data_file):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        pcount_estimate = self.ds.midx.get_nparticles_bbox(dle, dre)
        if pcount_estimate > 1e9:
            mylog.warning("Filtering %i particles to find total."
                          % pcount_estimate + \
                          " You may want to reconsider your bounding box.")
        pcount = 0
        for dd in self.ds.midx.iter_bbox_data(
            dle, dre,
            ['x']):
            pcount += dd['x'].size
        return {'dark_matter': pcount}

    def _identify_fields(self, data_file):
        fields = [("dark_matter", v) for v in self._handle.keys()]
        fields.append(("dark_matter", "mass"))
        return fields, {}


class IOHandlerSIndexHTTPSDF(IOHandlerSIndexSDF):
    _dataset_type = "midx_http_sdf_particles"

