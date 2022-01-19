import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseParticleIOHandler


class IOHandlerSDF(BaseParticleIOHandler):
    _dataset_type = "sdf_particles"

    @property
    def _handle(self):
        return self.ds.sdf_container

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set()
        assert len(ptf) == 1
        assert ptf.keys()[0] == "dark_matter"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert len(data_files) == 1
        for _data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            yield "dark_matter", (
                self._handle["x"],
                self._handle["y"],
                self._handle["z"],
            ), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set()
        assert len(ptf) == 1
        assert ptf.keys()[0] == "dark_matter"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert len(data_files) == 1
        for _data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            for ptype, field_list in sorted(ptf.items()):
                x = self._handle["x"]
                y = self._handle["y"]
                z = self._handle["z"]
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None:
                    continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.ds.parameters["particle_mass"]
                    else:
                        data = self._handle[field][mask]
                    yield (ptype, field), data

    def _identify_fields(self, data_file):
        fields = [("dark_matter", v) for v in self._handle.keys()]
        fields.append(("dark_matter", "mass"))
        return fields, {}

    def _count_particles(self, data_file):
        pcount = self._handle["x"].size
        if pcount > 1e9:
            mylog.warning(
                "About to load %i particles into memory. "
                "You may want to consider a midx-enabled load",
                pcount,
            )
        return {"dark_matter": pcount}


class IOHandlerHTTPSDF(IOHandlerSDF):
    _dataset_type = "http_sdf_particles"

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set()
        assert len(ptf) == 1
        assert ptf.keys()[0] == "dark_matter"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert len(data_files) == 1
        for _data_file in data_files:
            pcount = self._handle["x"].size
            yield "dark_matter", (
                self._handle["x"][:pcount],
                self._handle["y"][:pcount],
                self._handle["z"][:pcount],
            ), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set()
        assert len(ptf) == 1
        assert ptf.keys()[0] == "dark_matter"
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        assert len(data_files) == 1
        for _data_file in data_files:
            pcount = self._handle["x"].size
            for ptype, field_list in sorted(ptf.items()):
                x = self._handle["x"][:pcount]
                y = self._handle["y"][:pcount]
                z = self._handle["z"][:pcount]
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None:
                    continue
                for field in field_list:
                    if field == "mass":
                        if self.ds.field_info._mass_field is None:
                            pm = 1.0
                            if "particle_mass" in self.ds.parameters:
                                pm = self.ds.parameters["particle_mass"]
                            else:
                                raise RuntimeError
                            data = pm * np.ones(mask.sum(), dtype="float64")
                        else:
                            data = self._handle[self.ds.field_info._mass_field][:][mask]
                    else:
                        data = self._handle[field][:][mask]
                    yield (ptype, field), data

    def _count_particles(self, data_file):
        return {"dark_matter": self._handle["x"].http_array.shape}


class IOHandlerSIndexSDF(IOHandlerSDF):
    _dataset_type = "midx_sdf_particles"

    def _read_particle_coords(self, chunks, ptf):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        for dd in self.ds.midx.iter_bbox_data(dle, dre, ["x", "y", "z"]):
            yield "dark_matter", (dd["x"], dd["y"], dd["z"]), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        required_fields = []
        for field_list in sorted(ptf.values()):
            for field in field_list:
                if field == "mass":
                    continue
                required_fields.append(field)

        for dd in self.ds.midx.iter_bbox_data(dle, dre, required_fields):

            for ptype, field_list in sorted(ptf.items()):
                x = dd["x"]
                y = dd["y"]
                z = dd["z"]
                mask = selector.select_points(x, y, z, 0.0)
                del x, y, z
                if mask is None:
                    continue
                for field in field_list:
                    if field == "mass":
                        data = np.ones(mask.sum(), dtype="float64")
                        data *= self.ds.parameters["particle_mass"]
                    else:
                        data = dd[field][mask]
                    yield (ptype, field), data

    def _count_particles(self, data_file):
        dle = self.ds.domain_left_edge.in_units("code_length").d
        dre = self.ds.domain_right_edge.in_units("code_length").d
        pcount_estimate = self.ds.midx.get_nparticles_bbox(dle, dre)
        if pcount_estimate > 1e9:
            mylog.warning(
                "Filtering %i particles to find total. "
                "You may want to reconsider your bounding box.",
                pcount_estimate,
            )
        pcount = 0
        for dd in self.ds.midx.iter_bbox_data(dle, dre, ["x"]):
            pcount += dd["x"].size
        return {"dark_matter": pcount}

    def _identify_fields(self, data_file):
        fields = [("dark_matter", v) for v in self._handle.keys()]
        fields.append(("dark_matter", "mass"))
        return fields, {}


class IOHandlerSIndexHTTPSDF(IOHandlerSIndexSDF):
    _dataset_type = "midx_http_sdf_particles"
