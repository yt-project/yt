import numpy as np

from yt.frontends.sph.io import IOHandlerSPH


class IOHandlerSwift(IOHandlerSPH):
    _dataset_type = "swift"

    def __init__(self, ds, *args, **kwargs):
        super().__init__(ds, *args, **kwargs)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    # NOTE: we refer to sub_files in the next sections, these sub_files may
    # actually be full data_files.
    # In the event data_files are too big, yt breaks them up into sub_files and
    # we sort of treat them as files in the chunking system
    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        # yt has the concept of sub_files, i.e, we break up big files into
        # virtual sub_files to deal with the chunking system
        cname = self.ds._particle_coordinates_name
        for sub_file in self._sorted_chunk_iterator(chunks):
            with sub_file.transaction() as f:
                # This does NOT double-read anymore
                for ptype in sorted(ptf):
                    if sub_file.total_particles[ptype] == 0:
                        continue
                    pos = sub_file._read_field(ptype, cname, handle=f)
                    pos = pos.astype("float64", copy=False)
                    if ptype == self.ds._sph_ptypes[0]:
                        hsml = self._get_smoothing_length(sub_file, handle=f)
                    else:
                        hsml = 0.0
                    yield ptype, (pos[:, 0], pos[:, 1], pos[:, 2]), hsml

    def _yield_coordinates(self, sub_file, needed_ptype=None):
        cname = self.ds._particle_coordinates_name
        with sub_file.transaction() as f:
            for ptype, _ in sub_file._nonzero_ptypes():
                if needed_ptype and ptype != needed_ptype:
                    continue
                pos = sub_file._read_field(ptype, cname, handle=f)
                pos = pos.astype("float64", copy=False)
                yield ptype, pos

    def _get_smoothing_length(self, sub_file, pdtype=None, pshape=None, handle=None):
        # We do not need the pdtype and the pshape, but some frontends do so we
        # accept them and then just ignore them
        ptype = self.ds._sph_ptypes[0]
        hsml = sub_file._read_field(ptype, "SmoothingLength", handle=handle)
        return hsml.astype("float64", copy=False)

    def _read_particle_data_file(self, sub_file, ptf, selector=None):
        # note: this frontend uses the variable name and terminology sub_file.
        # other frontends use data_file with the understanding that it may
        # actually be a sub_file.
        return_data = {}

        cname = self.ds._particle_coordinates_name
        with sub_file.transaction() as f:
            for ptype, field_list, _ in sub_file._nonzero_ptf(ptf):

                if selector:
                    coords = sub_file._read_field(ptype, cname, handle=f)
                    if ptype == "PartType0":
                        hsmls = self._get_smoothing_length(sub_file, handle=f)
                    else:
                        hsmls = 0.0
                    mask = selector.select_points(
                        coords[:, 0], coords[:, 1], coords[:, 2], hsmls
                    )
                    del coords
                    if mask is None:
                        continue

                for field in field_list:
                    data = sub_file._read_field(ptype, field, handle=f)
                    if selector:
                        data = data[mask, ...]
                    data.astype("float64", copy=False)
                    return_data[(ptype, field)] = data

        return return_data

    def _count_particles(self, data_file):
        with data_file.transaction() as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")

        # if this data_file was a sub_file, then we just extract the region
        # defined by the subfile
        si, ei = data_file.start, data_file.end
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
        npart = {f"PartType{i}": v for i, v in enumerate(pcount)}
        return npart

    def _identify_fields(self, data_file):

        with data_file.transaction() as f:

            fields = []

            for key in f.keys():
                if not key.startswith("PartType"):
                    continue

                g = f[key]
                if self.ds._particle_coordinates_name not in g:
                    continue

                ptype = str(key)
                for k in g.keys():
                    kk = k
                    if str(kk) == self.ds._particle_mass_name:
                        fields.append((ptype, "Mass"))
                        continue
                    if not hasattr(g[kk], "shape"):
                        continue
                    if len(g[kk].shape) > 1:
                        self._vector_fields[kk] = g[kk].shape[1]
                    fields.append((ptype, str(kk)))

        return fields, {}
