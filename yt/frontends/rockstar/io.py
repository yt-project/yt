import os

import numpy as np

from yt.utilities import fortran_utils as fpu
from yt.utilities.io_handler import BaseParticleIOHandler

from .definitions import halo_dts, header_dt


def _can_load_with_format(
    filename: str, header_fmt: tuple[str, int, str], halo_format: np.dtype
) -> bool:
    with open(filename, "rb") as f:
        header = fpu.read_cattrs(f, header_fmt, "=")
        Nhalos = header["num_halos"]
        Nparttot = header["num_particles"]
        halos = np.fromfile(f, dtype=halo_format, count=Nhalos)

        # Make sure all masses are > 0
        if np.any(halos["particle_mass"] <= 0):
            return False
        # Make sure number of particles sums to expected value
        if halos["num_p"].sum() != Nparttot:
            return False

    return True


class IOHandlerRockstarBinary(BaseParticleIOHandler):
    _dataset_type = "rockstar_binary"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._halo_dt = self.detect_rockstar_format(
            self.ds.filename,
            self.ds.parameters["format_revision"],
        )

    @staticmethod
    def detect_rockstar_format(
        filename: str,
        guess: int | str,
    ) -> bool:
        revisions = list(halo_dts.keys())
        if guess in revisions:
            revisions.pop(revisions.index(guess))
        revisions = [guess] + revisions
        for revision in revisions:
            if _can_load_with_format(filename, header_dt, halo_dts[revision]):
                return halo_dts[revision]
        raise RuntimeError(f"Could not detect Rockstar format for file {filename}")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.

        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == "halos"
        ptype = "halos"
        for data_file in self._sorted_chunk_iterator(chunks):
            pcount = data_file.header["num_halos"]
            if pcount == 0:
                continue
            with open(data_file.filename, "rb") as f:
                pos = data_file._get_particle_positions(ptype, f=f)
                yield "halos", (pos[:, i] for i in range(3)), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == "halos"
        for data_file in self._sorted_chunk_iterator(chunks):
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
