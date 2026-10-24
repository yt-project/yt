import os
from functools import cached_property
from pathlib import Path

import numpy as np

from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import RAMSESOctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.cosmology import Cosmology

from .definitions import (
    FLUID_TYPES,
    LEGACY_PARTICLE_FIELDS,
    OUTPUT_DIR_RE,
    PARTICLE_PREFIXES,
    RAMSES2_FILE_RE,
)
from .fields import RAMSES2FieldInfo
from .io import particle_field_list

# Offsets (in units of the oct width) of the 8 cells of an oct, in the
# RAMSES ordering (x varies fastest).
_CELL_OFFSETS = np.array(
    [[(i >> idim) & 1 for idim in range(3)] for i in range(8)], dtype="int64"
)


def _encode_keys(ckey):
    """Pack an (N, 3) array of integer oct coordinates into int64 keys.

    Each coordinate uses 21 bits, which is enough for octs up to level 22.
    """
    ckey = ckey.astype("int64")
    return ckey[:, 0] | (ckey[:, 1] << 21) | (ckey[:, 2] << 42)


def _decode_keys(keys):
    """Inverse of _encode_keys."""
    mask = (1 << 21) - 1
    return np.stack([keys & mask, (keys >> 21) & mask, (keys >> 42) & mask], axis=1)


class RAMSES2FileSanitizer:
    """Find the output folder and info file from any file or folder
    that belongs to a RAMSES II (or mini-ramses) output."""

    def __init__(self, filename):
        self.original_filename = filename
        self.root_folder = None
        self.info_fname = None

        path = Path(filename)
        for candidate in (path, path.resolve()):
            if candidate.is_dir():
                folder = candidate
            elif candidate.is_file() and RAMSES2_FILE_RE.fullmatch(candidate.name):
                folder = candidate.parent
            else:
                continue
            if self._check_folder(folder):
                self.root_folder = folder
                self.info_fname = folder / "info.txt"
                break

    @staticmethod
    def _check_folder(folder):
        if not OUTPUT_DIR_RE.fullmatch(folder.name):
            return False
        info_fname = folder / "info.txt"
        if not info_fname.is_file():
            return False
        # Unlike RAMSES, the first entry of the info file is the number of
        # files the output is split into.
        try:
            with open(info_fname) as f:
                return f.readline().strip().startswith("nfile")
        except (OSError, UnicodeDecodeError):
            return False

    @property
    def is_valid(self):
        return self.root_folder is not None

    def validate(self):
        if not self.is_valid:
            raise ValueError(
                f"Could not determine RAMSES II output from {self.original_filename}"
            )


class RAMSES2DomainFile:
    """A single file of a RAMSES II output.

    RAMSES II writes one file per group of CPUs. Within each file,
    octs are stored level by level, which makes each file a contiguous
    chunk of the Hilbert curve at every level. Each file gets its own
    octree containing its octs, plus placeholder octs (attached to a
    dummy domain) needed to connect them to the root mesh or to mark
    cells refined by octs stored in another file.
    """

    def __init__(self, ds, domain_id):
        self.ds = ds
        self.domain_id = domain_id
        self.amr_fn = os.path.join(ds.root_folder, f"amr.{domain_id:05d}")

    def __repr__(self):
        return f"RAMSES2DomainFile: {self.domain_id}"

    def fluid_fn(self, prefix):
        return os.path.join(self.ds.root_folder, f"{prefix}.{self.domain_id:05d}")

    def particle_fn(self, prefix):
        return os.path.join(self.ds.root_folder, f"{prefix}.{self.domain_id:05d}")

    @cached_property
    def level_count(self):
        """Number of octs at each level (starting at levelmin) in this file."""
        nlevels = self.ds.parameters["levelmax"] - self.ds.parameters["levelmin"] + 1
        if not os.path.exists(self.amr_fn):
            return np.zeros(nlevels, dtype="int64")
        header = np.fromfile(self.amr_fn, dtype="<i4", count=3 + nlevels)
        ndim, levelmin, levelmax = header[:3]
        if (levelmin, levelmax) != (
            self.ds.parameters["levelmin"],
            self.ds.parameters["levelmax"],
        ):
            raise RuntimeError(
                f"Inconsistent levels in {self.amr_fn}: found "
                f"levelmin={levelmin}, levelmax={levelmax}."
            )
        return header[3:].astype("int64")

    @cached_property
    def level_offsets(self):
        """Index of the first oct of each level in this file."""
        return np.concatenate([[0], np.cumsum(self.level_count)[:-1]])

    @cached_property
    def max_level(self):
        """Deepest level (relative to levelmin) with octs in this file."""
        nonzero = np.flatnonzero(self.level_count)
        return int(nonzero[-1]) if nonzero.size else 0

    def fluid_header_size(self):
        # ndim, nvar, levelmin, levelmax, noct(levelmin:levelmax)
        return 4 * (4 + len(self.level_count))

    def _read_amr(self):
        """Return the oct coordinates and refinement maps at each level."""
        level_count = self.level_count
        noct = int(level_count.sum())
        header_size = 4 * (3 + len(level_count))
        if noct == 0:
            empty = [np.empty((0, 3), dtype="int32")] * len(level_count)
            return empty, [np.empty(0, dtype="int32")] * len(level_count)
        record_size = (os.path.getsize(self.amr_fn) - header_size) // max(noct, 1)

        # Recent versions store the refinement map of each oct as a bitmask
        # (one int32), older ones as one 4-byte logical per cell.
        nref = record_size // 4 - 3
        if nref not in (1, 8) or (
            record_size * noct + header_size != os.path.getsize(self.amr_fn)
        ):
            raise RuntimeError(f"Could not infer the oct record size of {self.amr_fn}")
        dt = np.dtype([("ckey", "<i4", (3,)), ("refined", "<i4", (nref,))])
        data = np.fromfile(self.amr_fn, dtype=dt, count=noct, offset=header_size)

        refined = data["refined"]
        if nref == 1:
            refined = refined[:, 0]
        else:
            refined = ((refined != 0) << np.arange(8)).sum(axis=1)

        ckeys = np.split(data["ckey"], np.cumsum(level_count)[:-1])
        refined = np.split(refined, np.cumsum(level_count)[:-1])
        return ckeys, refined

    @cached_property
    def oct_handler(self):
        ds = self.ds
        nlevels = len(self.level_count)

        if os.path.exists(self.amr_fn):
            ckeys, refined = self._read_amr()
        else:
            ckeys = [np.empty((0, 3), dtype="int32")] * nlevels
            refined = [np.empty(0, dtype="int32")] * nlevels
        own = [_encode_keys(ckey) for ckey in ckeys]

        # Placeholder octs. These are
        # 1. the children of refined cells stored in another file
        #    (so that the refined cells are not counted as leaves), and
        # 2. the parents of octs in this file stored in another file.
        other = [np.empty(0, dtype="int64") for _ in range(nlevels)]
        for ilevel in range(nlevels - 1):
            iparent, icell = np.nonzero((refined[ilevel][:, None] >> np.arange(8)) & 1)
            children = 2 * ckeys[ilevel][iparent].astype("int64") + _CELL_OFFSETS[icell]
            other[ilevel + 1] = np.setdiff1d(_encode_keys(children), own[ilevel + 1])
        for ilevel in range(nlevels - 1, 0, -1):
            keys = np.concatenate([own[ilevel], other[ilevel]])
            parents = _encode_keys(_decode_keys(keys) // 2)
            parents = np.setdiff1d(parents, own[ilevel - 1])
            other[ilevel - 1] = np.union1d(other[ilevel - 1], parents)

        # Placeholders all go into a dummy domain, after the last file
        dummy_domain = ds.nfile + 1
        domain_counts = np.zeros(dummy_domain, dtype="int64")
        domain_counts[self.domain_id - 1] = sum(o.size for o in own)
        domain_counts[dummy_domain - 1] = sum(o.size for o in other)
        root_nodes = own[0].size + other[0].size

        oct_handler = RAMSESOctreeContainer(
            ds.domain_dimensions // 2,
            ds.domain_left_edge,
            ds.domain_right_edge,
        )
        oct_handler.allocate_domains(domain_counts, root_nodes)
        mylog.debug(
            "Reading domain AMR % 4i (%s octs, %s placeholders)",
            self.domain_id,
            domain_counts[self.domain_id - 1],
            domain_counts[dummy_domain - 1],
        )

        # Octs are added level by level, so that their parents always exist.
        # Own octs are added in on-disk order, as the file index of an oct is
        # its position within its level.
        for ilevel in range(nlevels):
            dx = 1.0 / 2 ** (ilevel + ds.parameters["levelmin"] - 1)
            for domain, keys in (
                (self.domain_id, own[ilevel]),
                (dummy_domain, other[ilevel]),
            ):
                if keys.size == 0:
                    continue
                pos = (_decode_keys(keys) + 0.5) * dx
                oct_handler.add(domain, ilevel, np.asfortranarray(pos, dtype="float64"))
        oct_handler.finalize()
        return oct_handler

    def included(self, selector):
        if getattr(selector, "domain_id", None) is not None:
            return selector.domain_id == self.domain_id
        if not os.path.exists(self.amr_fn):
            # Without AMR data (e.g. output_amr=.false.), we cannot tell
            return True
        domain_ids = self.oct_handler.domain_identify(selector)
        return self.domain_id in domain_ids


class RAMSES2DomainSubset(OctreeSubset):
    _domain_offset = 1
    _block_order = "F"

    def __init__(self, base_region, domain, ds, num_zones=2, num_ghost_zones=0):
        if num_ghost_zones > 0:
            raise NotImplementedError(
                "Ghost zones are not supported yet by the RAMSES II frontend."
            )
        super().__init__(base_region, domain, ds, num_zones, num_ghost_zones)

    @property
    def oct_handler(self):
        return self.domain.oct_handler

    def fill(self, fn, fields, selector, field_list):
        """Read the selected cells of a fluid file.

        Each oct is stored on disk as a (nvar, 2**ndim) block of float32,
        so a cell is located by its oct index (level offset + position
        within the level), the variable index and its index within the oct.
        """
        oct_handler = self.oct_handler
        all_fields = [f for _, f in field_list]
        cell_count = selector.count_oct_cells(oct_handler, self.domain_id)
        data = {f: np.zeros(cell_count, dtype="float64") for _, f in fields}
        if cell_count == 0:
            return data

        level_inds, cell_inds, file_inds = oct_handler.file_index_octs(
            selector, self.domain_id, cell_count
        )
        oct_inds = self.domain.level_offsets[level_inds] + file_inds

        nvar = len(all_fields)
        noct = int(self.domain.level_count.sum())
        values = np.memmap(
            fn,
            dtype="<f4",
            mode="r",
            offset=self.domain.fluid_header_size(),
            shape=(noct, nvar, 8),
        )
        for _, fname in fields:
            ivar = all_fields.index(fname)
            data[fname][:] = values[oct_inds, ivar, cell_inds]
        return data

    def retrieve_ghost_zones(self, ngz, fields, smoothed=False):
        raise NotImplementedError(
            "Ghost zones are not supported yet by the RAMSES II frontend."
        )


class RAMSES2Index(OctreeIndex):
    def __init__(self, ds, dataset_type="ramses2"):
        self.fluid_field_list = ds._fields_in_file
        self.dataset_type = dataset_type
        super().__init__(ds, dataset_type)

    def _initialize_oct_handler(self):
        self.domains = [
            RAMSES2DomainFile(self.dataset, i + 1) for i in range(self.dataset.nfile)
        ]
        self.num_grids = sum(int(dom.level_count.sum()) for dom in self.domains)

    @cached_property
    def max_level(self):
        return max(dom.max_level for dom in self.domains)

    def _detect_output_fields(self):
        ds = self.dataset
        self.fluid_field_list = list(ds._fields_in_file)
        self.particle_field_list = [
            (ptype, fname)
            for ptype, fields in ds._particle_fields_in_file.items()
            for fname in fields
        ]
        self.field_list = self.particle_field_list + self.fluid_field_list

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            domains = [dom for dom in self.domains if dom.included(dobj.selector)]
            if len(domains) >= 1:
                mylog.debug("Identified %s intersecting domains", len(domains))
            base_region = getattr(dobj, "base_region", dobj)
            dobj._chunk_info = [
                RAMSES2DomainSubset(base_region, dom, self.dataset) for dom in domains
            ]
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            g = og.retrieve_ghost_zones(ngz, []) if ngz > 0 else og
            yield YTDataChunk(dobj, "spatial", [g], None)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None, cache=cache)

    def _get_particle_type_counts(self):
        return {
            ptype: sum(
                self.io._particle_count(dom.particle_fn(ptype)) for dom in self.domains
            )
            for ptype in self.dataset.particle_types_raw
        }


class RAMSES2Dataset(Dataset):
    _index_class = RAMSES2Index
    _field_info_class = RAMSES2FieldInfo
    _dataset_type = "ramses2"

    gamma = 1.4
    nfile = 1

    def __init__(
        self,
        filename,
        dataset_type="ramses2",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        sanitizer = RAMSES2FileSanitizer(filename)
        sanitizer.validate()
        self.root_folder = str(sanitizer.root_folder)

        self.fluid_types += tuple(FLUID_TYPES.values())
        self._fields_in_file = []
        self._particle_fields_in_file = {}
        self._particle_descriptors = {}
        self.storage_filename = storage_filename

        super().__init__(
            str(sanitizer.info_fname),
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

    @property
    def basename(self):
        return os.path.basename(self.root_folder)

    def __str__(self):
        return self.basename

    def _set_code_unit_attributes(self):
        boxlen = self.parameters["boxlen"]
        length_unit = self.parameters["unit_l"]
        density_unit = self.parameters["unit_d"]
        time_unit = self.parameters["unit_t"]

        mass_unit = density_unit * length_unit**3
        magnetic_unit = np.sqrt(4 * np.pi * mass_unit / (time_unit**2 * length_unit))
        pressure_unit = density_unit * (length_unit / time_unit) ** 2

        setdefaultattr(self, "density_unit", self.quan(density_unit, "g/cm**3"))
        setdefaultattr(self, "magnetic_unit", self.quan(magnetic_unit, "gauss"))
        setdefaultattr(self, "pressure_unit", self.quan(pressure_unit, "dyne/cm**2"))
        setdefaultattr(self, "time_unit", self.quan(time_unit, "s"))
        setdefaultattr(self, "mass_unit", self.quan(mass_unit, "g"))
        setdefaultattr(
            self, "velocity_unit", self.quan(length_unit, "cm") / self.time_unit
        )
        # Only the length unit gets scaled by boxlen, so that the
        # domain spans [0, 1] in code units, as in the RAMSES frontend.
        setdefaultattr(self, "length_unit", self.quan(length_unit * boxlen, "cm"))

    def _parse_parameter_file(self):
        params = {}
        with open(self.parameter_filename) as f:
            for line in f:
                if "=" not in line:
                    continue
                key, val = (s.strip() for s in line.split("=", 1))
                try:
                    params[key] = int(val)
                except ValueError:
                    params[key] = float(val)
        self.parameters.update(params)

        self.unique_identifier = int(os.stat(self.parameter_filename).st_ctime)
        self.nfile = params["nfile"]
        self.dimensionality = params["ndim"]
        if self.dimensionality != 3:
            raise NotImplementedError(
                "The RAMSES II frontend only supports 3D outputs, "
                f"got ndim={self.dimensionality}."
            )
        self.gamma = params["gamma"]
        self.refine_by = 2
        self.min_level = params["levelmin"] - 1
        self.max_level = params["levelmax"] - params["levelmin"]
        if params["levelmax"] > 22:
            raise NotImplementedError(
                "The RAMSES II frontend supports up to 22 levels of refinement, "
                f"got levelmax={params['levelmax']}."
            )
        self.domain_dimensions = np.ones(3, dtype="int64") * 2 ** params["levelmin"]
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype="float64")
        self._periodicity = (True, True, True)

        # Same criterion as in the RAMSES frontend
        is_cosmological = not (
            params["time"] >= 0 and params["H0"] == 1 and params["aexp"] == 1
        )
        if is_cosmological:
            self.cosmological_simulation = True
            self.current_redshift = 1.0 / params["aexp"] - 1.0
            self.omega_lambda = params["omega_l"]
            self.omega_matter = params["omega_m"]
            self.omega_radiation = 0.0
            self.hubble_constant = params["H0"] / 100.0
            co = Cosmology(
                hubble_constant=self.hubble_constant,
                omega_matter=self.omega_matter,
                omega_lambda=self.omega_lambda,
            )
            self.current_time = (
                co.t_from_z(self.current_redshift).to("s").d / (params["unit_t"])
            )
        else:
            self.cosmological_simulation = False
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.current_time = params["time"]

        self._detect_fluid_fields()
        self._detect_particle_fields()

    def _detect_fluid_fields(self):
        for prefix, ftype in FLUID_TYPES.items():
            header_fn = os.path.join(self.root_folder, f"{prefix}_header.txt")
            if not os.path.exists(header_fn):
                continue
            with open(header_fn) as f:
                for line in f:
                    if line.startswith("variable #"):
                        name = line.split(":", 1)[1].strip()
                        self._fields_in_file.append((ftype, name))

    def _detect_particle_fields(self):
        """Read the fields of each particle type from <prefix>_header.txt.

        Each particle type is named after the prefix of its files (part,
        star, sink, ...). Recent headers list one field per line as
        "name dtype ncomp"; older ones list the field names on one line.
        """
        for ptype in PARTICLE_PREFIXES:
            header_fn = os.path.join(self.root_folder, f"{ptype}_header.txt")
            if not os.path.exists(header_fn):
                continue
            with open(header_fn) as f:
                lines = [line.strip() for line in f]
            lines = [
                line for line in lines[lines.index("Particle fields") + 1 :] if line
            ]
            tokens = [line.split() for line in lines]
            # A legacy header is a single line of names, which can also have
            # three entries: check that the last token is a component count
            if all(len(t) == 3 and t[2].isdigit() for t in tokens):
                descriptors = [
                    (name, dtype, int(ncomp)) for name, dtype, ncomp in tokens
                ]
            else:
                names = lines[0].split()
                unknown = [name for name in names if name not in LEGACY_PARTICLE_FIELDS]
                if unknown:
                    mylog.warning(
                        "Skipping %s particles: unknown fields %s in %s",
                        ptype,
                        unknown,
                        header_fn,
                    )
                    continue
                descriptors = [(name, *LEGACY_PARTICLE_FIELDS[name]) for name in names]
            self._particle_descriptors[ptype] = descriptors
            self._particle_fields_in_file[ptype] = [
                yname for yname, _, _ in particle_field_list(descriptors)
            ]

        self.particle_types = self.particle_types_raw = tuple(
            self._particle_fields_in_file
        )

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        return RAMSES2FileSanitizer(filename).is_valid
