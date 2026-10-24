import os
from collections import defaultdict

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

from .definitions import FLUID_TYPES, PARTICLE_FIELD_NAMES

_FTYPE_TO_PREFIX = {ftype: prefix for prefix, ftype in FLUID_TYPES.items()}


def particle_field_list(descriptors):
    """Return the yt names of the fields of a particle file, given its
    (name, dtype, ncomp) descriptors, as (yt name, name, component)."""
    names = []
    for name, _, ncomp in descriptors:
        yname = PARTICLE_FIELD_NAMES.get(name, f"particle_{name}")
        if ncomp == 1:
            names.append((yname, name, 0))
        else:
            names.extend(
                (f"{yname}_{ax}", name, icomp) for icomp, ax in enumerate("xyz")
            )
    return names


class IOHandlerRAMSES2(BaseIOHandler):
    _dataset_type = "ramses2"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        tr = defaultdict(list)
        ftypes = {ft for ft, _ in fields}
        for chunk in chunks:
            for ft in ftypes:
                field_subs = [f for f in fields if f[0] == ft]
                field_list = [f for f in self.ds._fields_in_file if f[0] == ft]
                for subset in chunk.objs:
                    fn = subset.domain.fluid_fn(_FTYPE_TO_PREFIX[ft])
                    if not os.path.exists(fn):
                        continue
                    rv = subset.fill(fn, field_subs, selector, field_list)
                    for ftype, fname in field_subs:
                        d = rv.pop(fname)
                        if d.size == 0:
                            continue
                        mylog.debug(
                            "Filling %s with %s (%0.3e %0.3e) (%s zones)",
                            fname,
                            d.size,
                            d.min(),
                            d.max(),
                            d.size,
                        )
                        tr[ftype, fname].append(d)
        d = {}
        for field in fields:
            tmp = tr.pop(field, None)
            d[field] = np.concatenate(tmp) if tmp else np.empty(0, dtype="d")
        return d

    @staticmethod
    def _particle_count(fn):
        if not os.path.exists(fn):
            return 0
        return int(np.fromfile(fn, dtype="<i4", count=2)[1])

    def _particle_layout(self, fn, ptype):
        """Return the number of particles and the offset, dtype and number
        of components of each field of a particle file.

        Particle files contain a header (ndim, npart) followed by one
        contiguous array per field component, in the order of the header.
        """
        npart = self._particle_count(fn)
        if npart == 0:
            return 0, {}
        descriptors = self.ds._particle_descriptors[ptype]

        nbytes = os.path.getsize(fn) - 8
        # Old headers do not give the size of the ids: infer it from the file size
        n_id = sum(ncomp for _, dtype, ncomp in descriptors if dtype == "id")
        n_other = sum(
            ncomp * np.dtype(dtype).itemsize
            for _, dtype, ncomp in descriptors
            if dtype != "id"
        )
        id_size = (nbytes // npart - n_other) // max(n_id, 1)
        if nbytes != npart * (n_other + id_size * n_id) or (
            n_id > 0 and id_size not in (4, 8)
        ):
            raise RuntimeError(f"Could not infer the layout of particle file {fn}")

        layout = {}
        offset = 8
        for name, dtype, ncomp in descriptors:
            dtype = np.dtype(f"<i{id_size}" if dtype == "id" else dtype)
            layout[name] = (offset, dtype, ncomp)
            offset += ncomp * npart * dtype.itemsize
        return npart, layout

    def _read_particle_file(self, fn, ptype, fields):
        npart, layout = self._particle_layout(fn, ptype)
        lookup = {
            yname: (name, icomp)
            for yname, name, icomp in particle_field_list(
                self.ds._particle_descriptors[ptype]
            )
        }
        ds = self.ds
        tr = {}
        for yname in fields:
            name, icomp = lookup[yname]
            offset, dtype, _ = layout[name]
            data = np.fromfile(
                fn,
                dtype=dtype,
                count=npart,
                offset=offset + icomp * npart * dtype.itemsize,
            ).astype("float64")
            if name == "pos":
                # Positions are stored in [0, boxlen]. Wrap them periodically,
                # as float32 positions can round to exactly boxlen.
                data = np.mod(data / ds.parameters["boxlen"], 1.0)
            elif name == "angmom":
                # Lengths are in unit_l, whereas code_length is unit_l * boxlen
                data /= ds.parameters["boxlen"]
            elif name in ("birth_date", "merging_date") and ds.cosmological_simulation:
                # Times are stored as the proper time since today in units
                # of 1/H0 = unit_t / aexp**2
                aexp = ds.parameters["aexp"]
                data = ds.current_time.d + (data - ds.parameters["texp"]) / aexp**2
            tr[yname] = data
        return tr

    def _read_particle_coords(self, chunks, ptf):
        pn = "particle_position_%s"
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype in sorted(ptf):
                    fn = subset.domain.particle_fn(ptype)
                    if self._particle_count(fn) == 0:
                        continue
                    rv = self._read_particle_file(fn, ptype, [pn % ax for ax in "xyz"])
                    yield ptype, (rv[pn % "x"], rv[pn % "y"], rv[pn % "z"]), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = "particle_position_%s"
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    fn = subset.domain.particle_fn(ptype)
                    if self._particle_count(fn) == 0:
                        continue
                    pos_fields = [pn % ax for ax in "xyz"]
                    rv = self._read_particle_file(
                        fn, ptype, set(field_list) | set(pos_fields)
                    )
                    mask = selector.select_points(
                        rv[pn % "x"], rv[pn % "y"], rv[pn % "z"], 0.0
                    )
                    if mask is None:
                        continue
                    for fname in field_list:
                        yield (ptype, fname), rv[fname][mask]
