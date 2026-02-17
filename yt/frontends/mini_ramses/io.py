import os
import struct
from collections import defaultdict

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


class IOHandlerMiniRAMSES(BaseIOHandler):
    _dataset_type = "mini_ramses"
    _particle_reader = True

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Read fluid (hydro and gravity) data from mini-ramses output files
        # Use the selector to determine which cells to read
        tr = defaultdict(list)

        # Separate fields by type
        ftypes = {f[0] for f in fields}
        
        for chunk in chunks:
            for ft in ftypes:
                # Get all fields of this type
                field_subs = [f for f in fields if f[0] == ft]
                
                for subset in chunk.objs:
                    domain = subset.domain
                    ds = subset.ds
                    
                    # Determine file and field list based on field type
                    if ft == "mini-ramses":
                        fn = domain.hydro_fn
                        field_list = ds._fields_in_file
                    elif ft == "gravity":
                        fn = domain.grav_fn
                        field_list = ds._gravity_fields_in_file
                    else:
                        continue
                    
                    if not os.path.exists(fn):
                        continue
                    
                    # Use subset.fill to read with selection
                    rv = subset.fill(fn, field_subs, selector, field_list)
                    
                    for ftype, fname in field_subs:
                        d = rv.get(fname, np.empty(0, dtype="float64"))
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

        # Concatenate results
        rv = {}
        for field in fields:
            tmp = tr.pop(field, None)
            rv[field] = np.concatenate(tmp) if tmp else np.empty(0, dtype="float64")

        return rv

    def _read_particle_coords(self, chunks, ptf):
        for chunk in chunks:
            for subset in chunk.objs:
                domain = subset.domain
                ds = subset.ds
                ndim = ds.dimensionality
                boxlen = float(ds.domain_right_edge[0])

                for ptype in sorted(ptf):
                    prefix = self._ptype_to_prefix(ptype)
                    fn = os.path.join(
                        ds.root_folder,
                        f"{prefix}.{domain.domain_id:05d}",
                    )
                    if not os.path.exists(fn):
                        continue

                    pdata = self._read_particle_file(fn, ndim, prefix)
                    if pdata is None:
                        continue

                    npart = pdata["npart"]
                    pos = np.zeros((npart, 3), dtype="float64")
                    for i, ax in enumerate("xyz"[:ndim]):
                        pos[:, i] = pdata[f"pos_{ax}"]

                    # Normalize positions to code units
                    yield ptype, (
                        pos[:, 0] / boxlen,
                        pos[:, 1] / boxlen,
                        pos[:, 2] / boxlen,
                    ), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        for chunk in chunks:
            for subset in chunk.objs:
                domain = subset.domain
                ds = subset.ds
                ndim = ds.dimensionality
                boxlen = float(ds.domain_right_edge[0])

                for ptype, field_list in sorted(ptf.items()):
                    prefix = self._ptype_to_prefix(ptype)
                    fn = os.path.join(
                        ds.root_folder,
                        f"{prefix}.{domain.domain_id:05d}",
                    )
                    if not os.path.exists(fn):
                        continue

                    pdata = self._read_particle_file(fn, ndim, prefix)
                    if pdata is None:
                        continue

                    npart = pdata["npart"]
                    pos = np.zeros((npart, 3), dtype="float64")
                    for i, ax in enumerate("xyz"[:ndim]):
                        pos[:, i] = pdata[f"pos_{ax}"]

                    mask = selector.select_points(
                        pos[:, 0] / boxlen,
                        pos[:, 1] / boxlen,
                        pos[:, 2] / boxlen,
                        0.0,
                    )

                    if mask is None:
                        continue

                    for fname in field_list:
                        data = self._get_particle_field_data(
                            pdata, fname, ndim, boxlen
                        )
                        if data is not None:
                            yield (ptype, fname), data[mask]

    def _get_particle_field_data(self, pdata, fname, ndim, boxlen):
        """Map a yt particle field name to data from the particle file."""
        field_map = {
            "particle_position_x": "pos_x",
            "particle_position_y": "pos_y",
            "particle_position_z": "pos_z",
            "particle_velocity_x": "vel_x",
            "particle_velocity_y": "vel_y",
            "particle_velocity_z": "vel_z",
            "particle_mass": "mass",
            "particle_refinement_level": "level",
            "particle_identity": "birth_id",
            "particle_metallicity": "metallicity",
            "particle_birth_time": "birth_date",
        }

        if fname in field_map:
            key = field_map[fname]
            if key in pdata:
                return pdata[key].astype("float64")
        return None

    @staticmethod
    def _ptype_to_prefix(ptype):
        """Map yt particle type name to mini-ramses file prefix."""
        ptype_map = {
            "io": "part",
            "star": "star",
            "sink": "sink",
            "tree": "tree",
            "trac": "trac",
        }
        return ptype_map.get(ptype, "part")

    @staticmethod
    def _read_particle_file(fn, ndim, prefix):
        """Read a mini-ramses particle file (stream binary format).

        Mini-ramses output particle files use stream I/O with:
        - Header: ndim (int32) + npart (int32) = 8 bytes
        - Then contiguous arrays of float32 for positions, velocities, mass
        - Then int32 for level, int32 for birth_id
        """
        try:
            with open(fn, "rb") as f:
                file_ndim = struct.unpack("i", f.read(4))[0]
                npart = struct.unpack("i", f.read(4))[0]

                if npart == 0:
                    return None

                pdata = {"npart": npart}

                # Positions (float32)
                for ax in "xyz"[:ndim]:
                    pdata[f"pos_{ax}"] = np.frombuffer(
                        f.read(4 * npart), dtype="<f4"
                    ).copy()

                # Velocities (float32)
                for ax in "xyz"[:ndim]:
                    pdata[f"vel_{ax}"] = np.frombuffer(
                        f.read(4 * npart), dtype="<f4"
                    ).copy()

                # Mass (float32)
                pdata["mass"] = np.frombuffer(
                    f.read(4 * npart), dtype="<f4"
                ).copy()

                # Read remaining optional fields based on what's available
                remaining = os.path.getsize(fn) - f.tell()

                # Check for header to determine optional fields
                header_fn = os.path.join(
                    os.path.dirname(fn),
                    f"{prefix}_header.txt",
                )
                opt_fields = _detect_optional_particle_fields(
                    header_fn, prefix
                )

                for opt_name, opt_dtype, opt_size_per in opt_fields:
                    if remaining >= opt_size_per * npart:
                        pdata[opt_name] = np.frombuffer(
                            f.read(opt_size_per * npart), dtype=opt_dtype
                        ).copy()
                        remaining -= opt_size_per * npart

                # If no header was found, try to read level and birth_id
                if not opt_fields:
                    if remaining >= 4 * npart:
                        pdata["level"] = np.frombuffer(
                            f.read(4 * npart), dtype="<i4"
                        ).copy()
                        remaining -= 4 * npart
                    if remaining >= 4 * npart:
                        pdata["birth_id"] = np.frombuffer(
                            f.read(4 * npart), dtype="<i4"
                        ).copy()
                        remaining -= 4 * npart

                return pdata

        except (OSError, struct.error):
            return None

    def _count_particles(self, data_file):
        return {}

    def _identify_fields(self, data_file):
        return [], {}


def _detect_optional_particle_fields(header_fn, prefix):
    """Detect optional particle fields from header file."""
    fields = []
    if not os.path.exists(header_fn):
        return fields

    with open(header_fn) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "Particle fields" in line and i + 1 < len(lines):
                field_names = lines[i + 1].strip().split()
                # Skip pos, vel, mass (already read)
                for fname in field_names:
                    if fname in ("pos", "vel", "mass"):
                        continue
                    if fname == "potential":
                        fields.append(("potential", "<f4", 4))
                    elif fname == "metallicity":
                        fields.append(("metallicity", "<f4", 4))
                    elif fname == "accel":
                        fields.append(("accel_x", "<f4", 4))
                        fields.append(("accel_y", "<f4", 4))
                        fields.append(("accel_z", "<f4", 4))
                    elif fname == "angmom":
                        fields.append(("angmom_x", "<f4", 4))
                        fields.append(("angmom_y", "<f4", 4))
                        fields.append(("angmom_z", "<f4", 4))
                    elif fname == "birth_date":
                        fields.append(("birth_date", "<f4", 4))
                    elif fname == "merging_date":
                        fields.append(("merging_date", "<f4", 4))
                    elif fname == "level":
                        fields.append(("level", "<i4", 4))
                    elif fname == "birth_id":
                        fields.append(("birth_id", "<i4", 4))
                    elif fname == "merging_id":
                        fields.append(("merging_id", "<i8", 8))
                    elif fname == "tracking_id":
                        fields.append(("tracking_id", "<i8", 8))
                break
    return fields
