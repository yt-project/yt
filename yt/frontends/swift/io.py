import numpy as np

from yt.frontends.sph.io import IOHandlerSPH
from yt.utilities.on_demand_imports import _h5py as h5py


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
        chunks = list(chunks)
        sub_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                sub_files.update(obj.data_files)
        for sub_file in sorted(sub_files, key=lambda x: x.filename):
            si, ei = sub_file.start, sub_file.end
            f = h5py.File(sub_file.filename, mode="r")
            # This double-reads
            for ptype in sorted(ptf):
                if sub_file.total_particles[ptype] == 0:
                    continue
                pos = f[f"/{ptype}/Coordinates"][si:ei, :]
                pos = pos.astype("float64", copy=False)
                if ptype == self.ds._sph_ptypes[0]:
                    hsml = self._get_smoothing_length(sub_file)
                else:
                    hsml = 0.0
                yield ptype, (pos[:, 0], pos[:, 1], pos[:, 2]), hsml
            f.close()

    def _yield_coordinates(self, sub_file, needed_ptype=None):
        si, ei = sub_file.start, sub_file.end
        f = h5py.File(sub_file.filename, mode="r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        np.clip(pcount - si, 0, ei - si, out=pcount)
        pcount = pcount.sum()
        for key in f.keys():
            if (
                not key.startswith("PartType")
                or "Coordinates" not in f[key]
                or needed_ptype
                and key != needed_ptype
            ):
                continue
            pos = f[key]["Coordinates"][si:ei, ...]
            pos = pos.astype("float64", copy=False)
            yield key, pos
        f.close()

    def _get_smoothing_length(self, sub_file, pdtype=None, pshape=None):
        # We do not need the pdtype and the pshape, but some frontends do so we
        # accept them and then just ignore them
        ptype = self.ds._sph_ptypes[0]
        ind = int(ptype[-1])
        si, ei = sub_file.start, sub_file.end
        with h5py.File(sub_file.filename, mode="r") as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][ind].astype("int")
            pcount = np.clip(pcount - si, 0, ei - si)
            # we upscale to float64
            hsml = f[ptype]["SmoothingLength"][si:ei, ...]
            hsml = hsml.astype("float64", copy=False)
            return hsml

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        sub_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                sub_files.update(obj.data_files)

        for sub_file in sorted(sub_files, key=lambda x: x.filename):
            si, ei = sub_file.start, sub_file.end
            f = h5py.File(sub_file.filename, mode="r")
            for ptype, field_list in sorted(ptf.items()):
                if sub_file.total_particles[ptype] == 0:
                    continue
                g = f[f"/{ptype}"]
                # this should load as float64
                coords = g["Coordinates"][si:ei]
                if ptype == "PartType0":
                    hsmls = self._get_smoothing_length(sub_file)
                else:
                    hsmls = 0.0
                mask = selector.select_points(
                    coords[:, 0], coords[:, 1], coords[:, 2], hsmls
                )
                del coords
                if mask is None:
                    continue
                for field in field_list:
                    if field in ("Mass", "Masses"):
                        data = g[self.ds._particle_mass_name][si:ei][mask, ...]
                    else:
                        data = g[field][si:ei][mask, ...]

                    data.astype("float64", copy=False)
                    yield (ptype, field), data
            f.close()

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        f = h5py.File(data_file.filename, mode="r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        f.close()
        # if this data_file was a sub_file, then we just extract the region
        # defined by the subfile
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
        npart = {f"PartType{i}": v for i, v in enumerate(pcount)}
        return npart

    def _identify_fields(self, data_file):
        f = h5py.File(data_file.filename, mode="r")
        fields = []
        cname = self.ds._particle_coordinates_name  # Coordinates
        mname = self.ds._particle_mass_name  # Coordinates

        for key in f.keys():
            if not key.startswith("PartType"):
                continue

            g = f[key]
            if cname not in g:
                continue

            ptype = str(key)
            for k in g.keys():
                kk = k
                if str(kk) == mname:
                    fields.append((ptype, "Mass"))
                    continue
                if not hasattr(g[kk], "shape"):
                    continue
                if len(g[kk].shape) > 1:
                    self._vector_fields[kk] = g[kk].shape[1]
                fields.append((ptype, str(kk)))

        f.close()
        return fields, {}
