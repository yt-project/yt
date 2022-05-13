"""
AdaptaHOP data-file handling function




"""


from functools import partial
from operator import attrgetter
from typing import List, Tuple, Union

import numpy as np

from yt.utilities.cython_fortran_utils import FortranFile
from yt.utilities.io_handler import BaseParticleIOHandler

from .definitions import ATTR_T


class IOHandlerAdaptaHOPBinary(BaseParticleIOHandler):
    _dataset_type = "adaptahop_binary"

    _offsets = None  # Location of halos in the file
    _particle_positions = None  # Buffer of halo position

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError

    def _yield_coordinates(self, data_file):
        yield "halos", self._get_particle_positions()

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
        for data_file in sorted(data_files, key=attrgetter("filename")):
            pcount = (
                data_file.ds.parameters["nhalos"] + data_file.ds.parameters["nsubs"]
            )
            if pcount == 0:
                continue
            pos = self._get_particle_positions()
            yield ptype, [pos[:, i] for i in range(3)], 0.0

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

        def iterate_over_attributes(attr_list):
            for attr, *_ in attr_list:
                if isinstance(attr, tuple):
                    yield from attr
                else:
                    yield attr

        halo_attributes = self.ds._halo_attributes

        attr_pos = partial(_find_attr_position, halo_attributes=halo_attributes)

        for data_file in sorted(data_files, key=attrgetter("filename")):
            pcount = (
                data_file.ds.parameters["nhalos"] + data_file.ds.parameters["nsubs"]
            )
            if pcount == 0:
                continue
            ptype = "halos"
            field_list0 = sorted(ptf[ptype], key=attr_pos)
            field_list_pos = [f"raw_position_{k}" for k in "xyz"]
            field_list = sorted(set(field_list0 + field_list_pos), key=attr_pos)

            with FortranFile(self.ds.parameter_filename) as fpu:
                params = fpu.read_attrs(self.ds._header_attributes)

                todo = _todo_from_attributes(field_list, self.ds._halo_attributes)

                nhalos = params["nhalos"] + params["nsubs"]
                data = np.zeros((nhalos, len(field_list)))
                for ihalo in range(nhalos):
                    jj = 0
                    for it in todo:
                        if isinstance(it, int):
                            fpu.skip(it)
                        else:
                            tmp = fpu.read_attrs(it)
                            for key in iterate_over_attributes(it):
                                v = tmp[key]
                                if key not in field_list:
                                    continue
                                data[ihalo, jj] = v
                                jj += 1
            ipos = [field_list.index(k) for k in field_list_pos]
            w = self.ds.domain_width.to("code_length")[0].value / 2
            x, y, z = (data[:, i] + w for i in ipos)
            mask = selector.select_points(x, y, z, 0.0)
            del x, y, z

            if mask is None:
                continue
            for field in field_list0:
                i = field_list.index(field)
                yield (ptype, field), data[mask, i]

    def _count_particles(self, data_file):
        nhalos = data_file.ds.parameters["nhalos"] + data_file.ds.parameters["nsubs"]
        return {"halos": nhalos}

    def _identify_fields(self, data_file):
        fields = []

        for attr, _1, _2 in self.ds._halo_attributes:
            if isinstance(attr, str):
                fields.append(("halos", attr))
            else:
                for a in attr:
                    fields.append(("halos", a))
        return fields, {}

    # -----------------------------------------------------
    # Specific to AdaptaHOP
    def _get_particle_positions(self):
        """Read the particles and return them in code_units"""
        data = getattr(self, "_particle_positions", None)
        if data is not None:
            return data

        with FortranFile(self.ds.parameter_filename) as fpu:
            params = fpu.read_attrs(self.ds._header_attributes)

            todo = _todo_from_attributes(
                (
                    "particle_identifier",
                    "raw_position_x",
                    "raw_position_y",
                    "raw_position_z",
                ),
                self.ds._halo_attributes,
            )

            nhalos = params["nhalos"] + params["nsubs"]
            data = np.zeros((nhalos, 3))
            offset_map = np.zeros((nhalos, 2), dtype="int64")
            for ihalo in range(nhalos):
                ipos = fpu.tell()
                for it in todo:
                    if isinstance(it, int):
                        fpu.skip(it)
                    elif it[0][0] != "particle_identifier":
                        # Small optimisation here: we can read as vector
                        # dt = fpu.read_attrs(it)
                        # data[ihalo, 0] = dt['particle_position_x']
                        # data[ihalo, 1] = dt['particle_position_y']
                        # data[ihalo, 2] = dt['particle_position_z']
                        data[ihalo, :] = fpu.read_vector(it[0][-1])
                    else:
                        halo_id = fpu.read_int()
                        offset_map[ihalo, 0] = halo_id
                        offset_map[ihalo, 1] = ipos
        data = self.ds.arr(data, "code_length") + self.ds.domain_width / 2

        # Make sure halos are loaded in increasing halo_id order
        assert np.all(np.diff(offset_map[:, 0]) > 0)

        # Cache particle positions as one do not expect a large number of halos anyway
        self._particle_positions = data
        self._offsets = offset_map
        return data

    def members(self, ihalo):
        offset = self._offsets[ihalo, 1]
        todo = _todo_from_attributes(("particle_identities",), self.ds._halo_attributes)
        with FortranFile(self.ds.parameter_filename) as fpu:
            fpu.seek(offset)
            if isinstance(todo[0], int):
                fpu.skip(todo.pop(0))
            members = fpu.read_attrs(todo.pop(0))["particle_identities"]
        return members


def _todo_from_attributes(attributes: ATTR_T, halo_attributes: ATTR_T):
    # Helper function to generate a list of read-skip instructions given a list of
    # attributes. This is used to skip fields most of the fields when reading
    # the tree_brick files.
    iskip = 0
    todo: List[Union[int, List[Tuple[Union[Tuple[str, ...], str], int, str]]]] = []

    attributes = tuple(set(attributes))

    for i, (attrs, l, k) in enumerate(halo_attributes):
        attrs_list: Tuple[str, ...]
        if isinstance(attrs, tuple):
            if not all(isinstance(a, str) for a in attrs):
                raise TypeError(f"Expected a single str or a tuple of str, got {attrs}")
            attrs_list = attrs
        else:
            attrs_list = (attrs,)
        ok = False
        for attr in attrs_list:
            if attr in attributes:
                ok = True
                break

        if i == 0:
            if ok:
                state = "read"
                todo.append([])
            else:
                state = "skip"

        if ok:
            if state == "skip":
                # Switched from skip to read, store skip information and start
                # new read list
                todo.append(iskip)
                todo.append([])
                iskip = 0
            if not isinstance(todo[-1], list):
                raise TypeError
            todo[-1].append((attrs, l, k))
            state = "read"
        else:
            iskip += 1
            state = "skip"

    if state == "skip" and iskip > 0:
        todo.append(iskip)

    return todo


def _find_attr_position(key, halo_attributes: ATTR_T) -> int:
    j = 0
    for attrs, *_ in halo_attributes:
        if not isinstance(attrs, tuple):
            attrs = (attrs,)
        for a in attrs:
            if key == a:
                return j
            j += 1
    raise KeyError
