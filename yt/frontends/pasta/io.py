import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseIOHandler


class IOHandlerPasta(BaseIOHandler):
    _dataset_type = "pasta"

    def __init__(self, ds):
        super().__init__(ds)
        self.filename = ds.parameter_filename
        self._head = ds._head
        self._block_info = ds._block_info
        self._var_map = {name: i for i, name in enumerate(ds._var_names)}
        self._particle_cache = {}

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        rv = {}
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")
        if size == 0:
            return rv

        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )

        ind = dict.fromkeys(fields, 0)
        for chunk in chunks:
            data = self._read_chunk_data(chunk, fields)
            for g in chunk.objs:
                for field in fields:
                    ind[field] += g.select(
                        selector, data[g.id][field], rv[field], ind[field]
                    )
        return rv

    def _read_chunk_data(self, chunk, fields):
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        if len(fields) == 0:
            return rv

        field_sources = getattr(self.ds, "_field_sources", {})
        files_to_fields = {}
        for field in fields:
            ftype, fname = field
            if fname in field_sources:
                src = field_sources[fname]
                src_fn = src["filename"]
                if src_fn not in files_to_fields:
                    files_to_fields[src_fn] = []
                files_to_fields[src_fn].append((field, src))
            elif fname in self._var_map:
                if self.filename not in files_to_fields:
                    files_to_fields[self.filename] = []
                files_to_fields[self.filename].append(
                    (
                        field,
                        {
                            "filename": self.filename,
                            "head": self._head,
                            "block_info": self._block_info,
                            "v_idx": self._var_map[fname],
                            "is_face": False,
                            "is_edge": False,
                            "data_start": self._head["data_start"],
                            "coord_bytes": self.ds._coord_bytes,
                            "var_bytes": self.ds._var_bytes,
                            "block_bytes": self.ds._block_bytes,
                            "dfmt": self._head["byte_order"]
                            + ("f8" if self._head["data_size"] == 8 else "f4"),
                            "dsz": self._head["data_size"],
                        },
                    )
                )
            elif (
                fname.startswith("field_bf_")
                and self._head.get("num_face_field", 0) > 0
            ):
                if self.filename not in files_to_fields:
                    files_to_fields[self.filename] = []
                files_to_fields[self.filename].append(
                    (
                        field,
                        {
                            "filename": self.filename,
                            "head": self._head,
                            "block_info": self._block_info,
                            "v_idx": int(fname.split("_")[-1]) - 1,
                            "is_face": True,
                            "is_edge": False,
                            "data_start": self._head["data_start"],
                            "coord_bytes": self.ds._coord_bytes,
                            "var_bytes": self.ds._var_bytes,
                            "block_bytes": self.ds._block_bytes,
                            "dfmt": self._head["byte_order"]
                            + ("f8" if self._head["data_size"] == 8 else "f4"),
                            "dsz": self._head["data_size"],
                        },
                    )
                )

        info = self._block_info
        nc1, nc2, nc3 = info["ncells1"], info["ncells2"], info["ncells3"]
        tot_cells = nc1 * nc2 * nc3

        for src_fn, flist in files_to_fields.items():
            first_meta = flist[0][1]
            data_start = first_meta["data_start"]
            block_bytes = first_meta["block_bytes"]
            coord_bytes = first_meta["coord_bytes"]
            var_bytes = first_meta["var_bytes"]
            dfmt = first_meta["dfmt"]
            dsz = first_meta["dsz"]

            with open(src_fn, "rb") as f:
                for g in chunk.objs:
                    block_offset = data_start + g.id * block_bytes
                    for field, meta in flist:
                        if not meta["is_face"] and not meta["is_edge"]:
                            v_idx = meta["v_idx"]
                            var_offset = block_offset + coord_bytes + v_idx * var_bytes
                            f.seek(var_offset)
                            arr = np.fromfile(f, dtype=dfmt, count=tot_cells)
                            # Pasta shape on disk is (nc3, nc2, nc1); transpose to yt (nc1, nc2, nc3)
                            arr = arr.reshape((nc3, nc2, nc1)).T
                            rv[g.id][field] = arr.astype("float64")
                        elif meta["is_face"]:
                            arr = self._read_face_field(
                                f,
                                block_offset,
                                field[1],
                                meta["block_info"],
                                dsz,
                                dfmt,
                                coord_bytes,
                                var_bytes,
                            )
                            rv[g.id][field] = arr.astype("float64")
                        elif meta["is_edge"]:
                            arr = self._read_edge_field(
                                f,
                                block_offset,
                                field[1],
                                meta["block_info"],
                                dsz,
                                dfmt,
                                coord_bytes,
                                var_bytes,
                            )
                            rv[g.id][field] = arr.astype("float64")

        return rv

    def _read_face_field(
        self,
        f,
        block_offset,
        fname,
        info,
        dsz,
        dfmt,
        coord_bytes=None,
        var_bytes=None,
    ):
        if coord_bytes is None:
            coord_bytes = self.ds._coord_bytes
        if var_bytes is None:
            var_bytes = self.ds._var_bytes
        nc1, nc2, nc3 = info["ncells1"], info["ncells2"], info["ncells3"]
        face_start = block_offset + coord_bytes + info["nvar"] * var_bytes
        shapes = [
            ("field_bf_1", (nc3, nc2, nc1 + 1)),
            ("field_bf_2", (nc3, nc2 + 1, nc1)),
            ("field_bf_3", (nc3 + 1, nc2, nc1)),
        ]
        curr_offset = face_start
        for name, shape in shapes:
            count = shape[0] * shape[1] * shape[2]
            if name == fname:
                f.seek(curr_offset)
                raw = np.fromfile(f, dtype=dfmt, count=count).reshape(shape)
                # Average staggered component to cell centers
                for axis, ncell in ((0, nc3), (1, nc2), (2, nc1)):
                    if raw.shape[axis] == ncell + 1:
                        raw = 0.5 * (
                            np.take(raw, np.arange(ncell), axis=axis)
                            + np.take(raw, np.arange(1, ncell + 1), axis=axis)
                        )
                return raw.T
            curr_offset += count * dsz
        raise KeyError(f"Face field {fname} not found")

    def _read_edge_field(
        self,
        f,
        block_offset,
        fname,
        info,
        dsz,
        dfmt,
        coord_bytes=None,
        var_bytes=None,
    ):
        if coord_bytes is None:
            coord_bytes = self.ds._coord_bytes
        if var_bytes is None:
            var_bytes = self.ds._var_bytes
        nc1, nc2, nc3 = info["ncells1"], info["ncells2"], info["ncells3"]
        face_count = 0
        if info.get("num_face_field", 0) > 0:
            face_count = (
                nc3 * nc2 * (nc1 + 1) + nc3 * (nc2 + 1) * nc1 + (nc3 + 1) * nc2 * nc1
            )
        face_bytes = face_count * dsz
        edge_start = block_offset + coord_bytes + info["nvar"] * var_bytes + face_bytes
        shapes = [
            ("field_emf_1", (nc3 + 1, nc2 + 1, nc1)),
            ("field_emf_2", (nc3 + 1, nc2, nc1 + 1)),
            ("field_emf_3", (nc3, nc2 + 1, nc1 + 1)),
        ]
        curr_offset = edge_start
        for name, shape in shapes:
            count = shape[0] * shape[1] * shape[2]
            if name == fname:
                f.seek(curr_offset)
                raw = np.fromfile(f, dtype=dfmt, count=count).reshape(shape)
                for axis, ncell in ((0, nc3), (1, nc2), (2, nc1)):
                    if raw.shape[axis] == ncell + 1:
                        raw = 0.5 * (
                            np.take(raw, np.arange(ncell), axis=axis)
                            + np.take(raw, np.arange(1, ncell + 1), axis=axis)
                        )
                return raw.T
            curr_offset += count * dsz
        raise KeyError(f"Edge field {fname} not found")

    def _get_particle_prop(self, name):
        if name not in self._particle_cache:
            pinfo = getattr(self.ds, "_particle_info", None)
            if pinfo is None or name not in pinfo["prop_offsets"]:
                raise KeyError(
                    f"Particle property {name} not found in {getattr(self.ds, 'particle_filename', None)}"
                )
            offset = pinfo["prop_offsets"][name]
            dtype = pinfo["prop_dtypes"][name]
            npar = pinfo["npar"]
            with open(self.ds.particle_filename, "rb") as f:
                f.seek(offset)
                data = np.fromfile(f, dtype=dtype, count=npar)
            self._particle_cache[name] = data
        return self._particle_cache[name]

    def _get_particle_coords_raw(self):
        x1 = self._get_particle_prop("x1")
        x2 = self._get_particle_prop("x2")
        x3 = self._get_particle_prop("x3")
        return x1, x2, x3

    def _get_particle_field_data(self, field, start, end):
        pinfo = self.ds._particle_info
        geom = self.ds.geometry

        # Check raw property names
        if field in pinfo["prop_offsets"]:
            return self._get_particle_prop(field)[start:end]

        if field == "particle_index":
            return self._get_particle_prop("pid")[start:end]
        if field == "particle_mass":
            return self._get_particle_prop("mass")[start:end]

        # Cartesian positions
        if field in (
            "particle_position_x",
            "particle_position_y",
            "particle_position_z",
        ):
            x1, x2, x3 = self._get_particle_coords_raw()
            if geom == "cartesian":
                if field == "particle_position_x":
                    return x1[start:end]
                elif field == "particle_position_y":
                    return x2[start:end]
                else:
                    return x3[start:end]
            elif geom == "cylindrical":
                r = x1[start:end]
                phi = x2[start:end]
                z = x3[start:end]
                if field == "particle_position_x":
                    return r * np.cos(phi)
                elif field == "particle_position_y":
                    return r * np.sin(phi)
                else:
                    return z
            elif geom == "spherical":
                r = x1[start:end]
                th = x2[start:end]
                ph = x3[start:end]
                if field == "particle_position_x":
                    return r * np.sin(th) * np.cos(ph)
                elif field == "particle_position_y":
                    return r * np.sin(th) * np.sin(ph)
                else:
                    return r * np.cos(th)

        # Curvilinear positions
        if field == "particle_position_r":
            return self._get_particle_prop("x1")[start:end]
        if field == "particle_position_theta":
            return self._get_particle_prop("x2")[start:end]
        if field == "particle_position_phi":
            return self._get_particle_prop("x3")[start:end]

        # Cartesian velocities
        if field in (
            "particle_velocity_x",
            "particle_velocity_y",
            "particle_velocity_z",
        ):
            v1 = self._get_particle_prop("v1")[start:end]
            v2 = self._get_particle_prop("v2")[start:end]
            v3 = self._get_particle_prop("v3")[start:end]
            if geom == "cartesian":
                if field == "particle_velocity_x":
                    return v1
                elif field == "particle_velocity_y":
                    return v2
                else:
                    return v3
            elif geom == "cylindrical":
                phi = self._get_particle_prop("x2")[start:end]
                if field == "particle_velocity_x":
                    return v1 * np.cos(phi) - v2 * np.sin(phi)
                elif field == "particle_velocity_y":
                    return v1 * np.sin(phi) + v2 * np.cos(phi)
                else:
                    return v3
            elif geom == "spherical":
                th = self._get_particle_prop("x2")[start:end]
                ph = self._get_particle_prop("x3")[start:end]
                if field == "particle_velocity_x":
                    return (
                        v1 * np.sin(th) * np.cos(ph)
                        + v2 * np.cos(th) * np.cos(ph)
                        - v3 * np.sin(ph)
                    )
                elif field == "particle_velocity_y":
                    return (
                        v1 * np.sin(th) * np.sin(ph)
                        + v2 * np.cos(th) * np.sin(ph)
                        + v3 * np.cos(ph)
                    )
                else:
                    return v1 * np.cos(th) - v2 * np.sin(th)

        # Curvilinear velocities
        if field == "particle_velocity_r":
            return self._get_particle_prop("v1")[start:end]
        if field == "particle_velocity_theta":
            return self._get_particle_prop("v2")[start:end]
        if field == "particle_velocity_phi":
            return self._get_particle_prop("v3")[start:end]

        raise KeyError(f"Unknown particle field {field}")

    def _read_particle_coords(self, chunks, ptf):
        yield from (
            (ptype, xyz, 0.0)
            for ptype, xyz in self._read_particle_fields(chunks, ptf, None)
        )

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        pinfo = getattr(self.ds, "_particle_info", None)
        if pinfo is None or pinfo["npar"] == 0:
            return

        p_ind = self.ds.index.particle_indices
        c0, c1, c2 = self._get_particle_coords_raw()

        for chunk in chunks:
            for g in chunk.objs:
                count = g.NumberOfParticles
                if count == 0:
                    continue
                start = p_ind[g.id]
                end = p_ind[g.id + 1]

                gx = c0[start:end]
                gy = c1[start:end]
                gz = c2[start:end]

                for ptype, field_list in sorted(ptf.items()):
                    if selector is None:
                        yield ptype, (gx, gy, gz)
                        continue

                    mask = selector.select_points(gx, gy, gz, 0.0)
                    if mask is None:
                        continue

                    for field in field_list:
                        data = self._get_particle_field_data(field, start, end)
                        yield (ptype, field), data[mask]
