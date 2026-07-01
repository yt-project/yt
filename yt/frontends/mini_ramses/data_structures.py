import os
import struct
from pathlib import Path

import numpy as np

from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import RAMSESOctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex

from .definitions import (
    MINI_RAMSES_FILE_RE,
    OUTPUT_DIR_RE,
)
from .fields import MiniRAMSESFieldInfo


class MiniRAMSESFileSanitizer:
    """Handle the different files that can be passed and associated
    safely to a mini-ramses output."""

    root_folder = None
    info_fname = None

    def __init__(self, filename):
        paths_to_try = (Path(filename), Path(filename).resolve())

        self.original_filename = filename
        self.is_valid = False

        for path in paths_to_try:
            if self._test_path(path):
                break

    def _test_path(self, path):
        """Test if path leads to a valid mini-ramses output."""
        # If it's a file, check if it's a known mini-ramses file
        if path.is_file():
            if MINI_RAMSES_FILE_RE.match(path.name):
                return self._test_parent_folder(path.parent)
            # Could be the info.txt file directly
            if path.name == "info.txt":
                return self._test_parent_folder(path.parent)
            return False
        # If it's a directory, check if it's an output directory
        if path.is_dir():
            return self._test_parent_folder(path)
        return False

    def _test_parent_folder(self, folder):
        """Check if folder is a valid mini-ramses output directory."""
        match = OUTPUT_DIR_RE.match(folder.name)
        if not match:
            return False

        info_path = folder / "info.txt"
        if not info_path.exists():
            return False

        # Verify this is mini-ramses (not RAMSES) by checking for
        # mini-ramses-specific file patterns (e.g., amr.NNNNN not amr_NNNNN.outNNNNN)
        # and info.txt format (has "nfile" as first line)
        if not self._is_mini_ramses_info(info_path):
            return False

        self.root_folder = folder
        self.info_fname = info_path
        self.is_valid = True
        return True

    @staticmethod
    def _is_mini_ramses_info(info_path):
        """Check if this is a mini-ramses info file (has nfile as first field)."""
        try:
            with open(info_path) as f:
                first_line = f.readline().strip()
                return first_line.startswith("nfile")
        except (OSError, UnicodeDecodeError):
            return False

    def validate(self):
        """Raise an error if the file is not valid."""
        if not self.is_valid:
            raise ValueError(
                f"Cannot identify {self.original_filename} as a "
                "mini-ramses output"
            )


class MiniRAMSESDomainFile:
    """Manage a single domain (CPU file) of a mini-ramses output."""

    _last_mask = None
    _last_selector_id = None

    def __init__(self, ds, domain_id):
        self.ds = ds
        self.domain_id = domain_id
        self._level_count = None
        self._level_offsets = None
        self._octree = None

        iout = int(ds.basename.split("_")[1])
        basename = ds.root_folder

        # File paths for this domain
        icpu_str = f"{domain_id:05d}"
        self.amr_fn = os.path.join(basename, f"amr.{icpu_str}")
        self.hydro_fn = os.path.join(basename, f"hydro.{icpu_str}")
        self.grav_fn = os.path.join(basename, f"grav.{icpu_str}")

    @property
    def level_count(self):
        if self._level_count is not None:
            return self._level_count
        self._read_amr_header()
        return self._level_count

    @property
    def level_offsets(self):
        """Return cumulative oct count up to (but not including) each level."""
        if self._level_offsets is not None:
            return self._level_offsets
        self._read_amr_header()
        return self._level_offsets

    def _read_amr_header(self):
        """Read the AMR header from the stream-based binary file."""
        if not os.path.exists(self.amr_fn):
            self._level_count = np.zeros(
                self.ds.max_level - self.ds.min_level + 1, dtype="int64"
            )
            self._level_offsets = np.zeros_like(self._level_count)
            return

        with open(self.amr_fn, "rb") as f:
            ndim = struct.unpack("i", f.read(4))[0]
            levelmin = struct.unpack("i", f.read(4))[0]
            nlevelmax = struct.unpack("i", f.read(4))[0]

            noct = np.zeros(nlevelmax, dtype="int32")
            for ilevel in range(levelmin - 1, nlevelmax):
                noct[ilevel] = struct.unpack("i", f.read(4))[0]

        nlevels = self.ds.max_level - self.ds.min_level + 1
        self._level_count = np.zeros(nlevels, dtype="int64")
        for ilevel in range(levelmin - 1, nlevelmax):
            idx = ilevel - (self.ds.min_level - 1)
            if 0 <= idx < nlevels:
                self._level_count[idx] = noct[ilevel]

        # Compute cumulative offsets: offset[i] = sum of octs at levels < i
        self._level_offsets = np.zeros(nlevels, dtype="int64")
        cumsum = 0
        for i in range(nlevels):
            self._level_offsets[i] = cumsum
            cumsum += self._level_count[i]

    def _read_hydro_header(self):
        """Read the hydro header to get nvar."""
        if not os.path.exists(self.hydro_fn):
            return 0
        with open(self.hydro_fn, "rb") as f:
            ndim = struct.unpack("i", f.read(4))[0]
            nvar = struct.unpack("i", f.read(4))[0]
        return nvar

    def _read_grav_header(self):
        """Read the gravity header to get nvar."""
        if not os.path.exists(self.grav_fn):
            return 0
        with open(self.grav_fn, "rb") as f:
            ndim = struct.unpack("i", f.read(4))[0]
            nvar = struct.unpack("i", f.read(4))[0]
        return nvar


class MiniRAMSESDomainSubset(OctreeSubset):
    _domain_offset = 1
    _block_order = "F"

    _base_domain = None

    def __init__(
        self, base_region, domain, ds, over_refine_factor=1, num_ghost_zones=0
    ):
        super().__init__(
            base_region, domain, ds, over_refine_factor, num_ghost_zones
        )
        self._base_domain = domain

    @property
    def oct_handler(self):
        return self.ds.index.oct_handler

    def fill(self, fn, fields, selector, field_list):
        """Read and fill field data from a mini-ramses binary file.
        
        Parameters
        ----------
        fn : str
            Filename to read from
        fields : list
            List of (ftype, fname) tuples to read
        selector : SelectorObject
            The selector determining which cells to include
        field_list : list
            List of all fields in the file (for determining indices)
            
        Returns
        -------
        data : dict
            Dictionary of {fname: array} for each requested field
        """
        ndim = self.ds.dimensionality
        twotondim = 2**ndim
        oct_handler = self.oct_handler
        domain = self._base_domain
        
        # Count selected cells
        cell_count = selector.count_oct_cells(oct_handler, self.domain_id)
        
        # Get field names only (without ftype)
        field_names = [f for ft, f in fields]
        
        # Initialize output data
        data = {}
        for fname in field_names:
            data[fname] = np.zeros(cell_count, dtype="float64")
        
        if cell_count == 0:
            return data
        
        # Get indices of selected cells
        # level_inds: which level (0-indexed relative to min_level)
        # cell_inds: which cell within the oct (0-7)
        # file_inds: which oct within that level
        level_inds, cell_inds, file_inds = oct_handler.file_index_octs(
            selector, self.domain_id, cell_count
        )
        
        # Read all data from file into a flat array
        all_data = self._read_all_cell_data(fn, field_list)
        
        if all_data is None:
            return data
        
        # Get level offsets to compute global oct index
        level_offsets = domain.level_offsets
        
        # Compute global oct index: offset for this level + oct index within level
        # file_inds gives oct index within each level
        # level_offsets[level] gives cumulative count of octs at levels < level
        global_oct_inds = np.zeros(cell_count, dtype="int64")
        for i in range(cell_count):
            level = level_inds[i]
            if level < len(level_offsets):
                global_oct_inds[i] = level_offsets[level] + file_inds[i]
            else:
                global_oct_inds[i] = file_inds[i]
        
        # Cell index in flat array: oct_idx * 8 + cell_idx
        indices = global_oct_inds * twotondim + cell_inds
        
        # Select the requested cells
        for ftype, fname in fields:
            if fname in all_data:
                data[fname] = all_data[fname][indices]
        
        return data
    
    def _read_all_cell_data(self, fn, field_list):
        """Read all cell data from file into arrays indexed by oct.
        
        Returns dict with field name -> flat array where index = oct_idx * 8 + cell_idx
        """
        if not os.path.exists(fn):
            return None
            
        ndim = self.ds.dimensionality
        twotondim = 2**ndim
        
        with open(fn, "rb") as f:
            # Read header
            f_ndim = struct.unpack("i", f.read(4))[0]
            f_nvar = struct.unpack("i", f.read(4))[0]
            levelmin = struct.unpack("i", f.read(4))[0]
            nlevelmax = struct.unpack("i", f.read(4))[0]
            
            noct = np.zeros(nlevelmax, dtype="int32")
            for ilevel in range(levelmin - 1, nlevelmax):
                noct[ilevel] = struct.unpack("i", f.read(4))[0]
            
            total_octs = noct.sum()
            total_cells = total_octs * twotondim
            
            # Build field index map
            all_field_names = [f for ft, f in field_list]
            field_indices = {}
            for ftype, fname in field_list:
                if fname in all_field_names:
                    try:
                        fidx = all_field_names.index(fname)
                        if fidx < f_nvar:
                            field_indices[fname] = fidx
                    except ValueError:
                        pass
            
            # Pre-allocate arrays
            raw_data = {fname: np.zeros(total_cells, dtype="float64") 
                       for fname in field_indices}
            
            # Read all data
            # Layout is: for each oct, nvar blocks of 8 cells each (field-major within oct)
            cell_offset = 0
            for ilevel in range(levelmin - 1, nlevelmax):
                ncache = noct[ilevel]
                if ncache == 0:
                    continue
                    
                for igrid in range(ncache):
                    # Each oct: nvar fields * twotondim cells * 4 bytes
                    # Layout is (nvar, twotondim) - field-major
                    oct_data = np.frombuffer(
                        f.read(4 * twotondim * f_nvar),
                        dtype="<f4"
                    ).reshape(f_nvar, twotondim)
                    
                    # Store each field's data for this oct's cells
                    for fname, fidx in field_indices.items():
                        for icell in range(twotondim):
                            raw_data[fname][cell_offset + icell] = oct_data[fidx, icell]
                    
                    cell_offset += twotondim
            
        return raw_data

    def retrieve_ghost_zones(self, ngz, fields, smoothed=False):
        new_subset = MiniRAMSESDomainSubset(
            self.base_region,
            self._base_domain,
            self.ds,
            self._num_ghost_zones,
            num_ghost_zones=ngz,
        )
        return new_subset


class MiniRAMSESIndex(OctreeIndex):
    def __init__(self, ds, dataset_type="mini_ramses"):
        self.fluid_field_list = ds._fields_in_file
        self.over_refine_factor = 1
        self.dataset_type = dataset_type
        super().__init__(ds, dataset_type)

    def _initialize_oct_handler(self):
        """Build the octree from mini-ramses AMR data."""
        ds = self.dataset
        ndomains = ds.nfile

        self.domains = [
            MiniRAMSESDomainFile(ds, i + 1) for i in range(ndomains)
        ]

        # domain_counts is an array with oct count per domain
        domain_counts = np.array([d.level_count.sum() for d in self.domains])
        # Root nodes are octs at the minimum level (index 0 in level_count)
        root_nodes = sum(int(d.level_count[0]) for d in self.domains)
        mylog.info("Allocating %s octs", domain_counts.sum())

        self.oct_handler = RAMSESOctreeContainer(
            ds.domain_dimensions / 2,
            ds.domain_left_edge,
            ds.domain_right_edge,
        )
        self.oct_handler.allocate_domains(domain_counts, root_nodes)

        # Read AMR data from each domain and add to octree
        for dom in self.domains:
            self._read_amr_data(dom)

        self.oct_handler.finalize()

    def _read_amr_data(self, domain):
        """Read AMR data from a mini-ramses domain file and populate octree."""
        if not os.path.exists(domain.amr_fn):
            return

        ds = self.dataset
        ndim = ds.dimensionality
        twotondim = 2**ndim

        with open(domain.amr_fn, "rb") as f:
            ndim_file = struct.unpack("i", f.read(4))[0]
            levelmin = struct.unpack("i", f.read(4))[0]
            nlevelmax = struct.unpack("i", f.read(4))[0]

            noct = np.zeros(nlevelmax, dtype="int32")
            for ilevel in range(levelmin - 1, nlevelmax):
                noct[ilevel] = struct.unpack("i", f.read(4))[0]

                # Read grid data: for each level, for each grid:
            # ckey (ndim ints) + refined (twotondim ints)
            for ilevel in range(levelmin - 1, nlevelmax):
                ncache = noct[ilevel]
                if ncache == 0:
                    continue

                # ilevel is 0-indexed (0 = level 1)
                # For RAMSES octree, we pass level relative to levelmin
                # ilevel=6 means level 7, and if levelmin=7, this is the root level (0)
                yt_level = ilevel - (levelmin - 1)

                # Read all grids at this level in batch
                pos_all = np.zeros((ncache, 3), dtype="float64", order="F")

                for igrid in range(ncache):
                    # Read Cartesian key (ndim int32 values)
                    ckey = np.array(
                        struct.unpack(f"{ndim}i", f.read(4 * ndim)),
                        dtype="int64",
                    )

                    # Read refinement map (twotondim int32 values)
                    refined = np.array(
                        struct.unpack(f"{twotondim}i", f.read(4 * twotondim)),
                        dtype="int32",
                    )

                    # Compute position in code units from Cartesian key
                    # ckey gives grid indices at this level
                    # For levelmin=7, root octs form a 2^(levelmin-1) = 64^3 grid
                    # Number of octs per dimension at level L is 2^(L-1)
                    # Position in normalized [0,1] = (ckey + 0.5) / 2^(L-1)
                    # Then scale by boxlen to get code units
                    boxlen = float(ds.domain_right_edge[0] - ds.domain_left_edge[0])
                    n_octs_per_dim = 2 ** ilevel  # For level 7 (ilevel=6): 64
                    dx = boxlen / n_octs_per_dim
                    for idim in range(ndim):
                        pos_all[igrid, idim] = (ckey[idim] + 0.5) * dx

                # Add all octs at this level to octree in batch
                self.oct_handler.add(
                    domain.domain_id,
                    yt_level,
                    pos_all,
                )

    def _detect_output_fields(self):
        ds = self.dataset
        # Combine hydro and gravity fields
        self.fluid_field_list = ds._fields_in_file + ds._gravity_fields_in_file

        # Particle fields - check for particle files
        self.particle_field_list = []
        ptype_map = {
            "part": "io",
            "star": "star",
            "sink": "sink",
        }

        for prefix, ptype in ptype_map.items():
            fn = os.path.join(ds.root_folder, f"{prefix}.00001")
            if os.path.exists(fn):
                # Read header to determine what fields exist
                header_fn = os.path.join(
                    ds.root_folder, f"{prefix}_header.txt"
                )
                pfields = self._detect_particle_fields(
                    fn, header_fn, prefix, ptype
                )
                for f in pfields:
                    self.particle_field_list.append(f)

        self.field_list = self.particle_field_list + self.fluid_field_list

    def _detect_particle_fields(self, fn, header_fn, prefix, ptype):
        """Detect available particle fields from header file."""
        fields = []
        ndim = self.dataset.dimensionality

        if os.path.exists(header_fn):
            with open(header_fn) as f:
                lines = f.readlines()
                # Look for "Particle fields" line
                for i, line in enumerate(lines):
                    if "Particle fields" in line:
                        if i + 1 < len(lines):
                            field_names = lines[i + 1].strip().split()
                            for fname in field_names:
                                if fname == "pos":
                                    for ax in "xyz"[:ndim]:
                                        fields.append(
                                            (ptype, f"particle_position_{ax}")
                                        )
                                elif fname == "vel":
                                    for ax in "xyz"[:ndim]:
                                        fields.append(
                                            (ptype, f"particle_velocity_{ax}")
                                        )
                                elif fname == "mass":
                                    fields.append((ptype, "particle_mass"))
                                elif fname == "level":
                                    fields.append(
                                        (ptype, "particle_refinement_level")
                                    )
                                elif fname == "birth_id":
                                    fields.append(
                                        (ptype, "particle_identity")
                                    )
                                elif fname == "metallicity":
                                    fields.append(
                                        (ptype, "particle_metallicity")
                                    )
                                elif fname == "birth_date":
                                    fields.append(
                                        (ptype, "particle_birth_time")
                                    )
                                elif fname == "potential":
                                    fields.append(
                                        (ptype, "particle_potential")
                                    )
                                elif fname == "accel":
                                    for ax in "xyz"[:ndim]:
                                        fields.append(
                                            (
                                                ptype,
                                                f"particle_acceleration_{ax}",
                                            )
                                        )
                                elif fname == "angmom":
                                    for ax in "xyz"[:ndim]:
                                        fields.append(
                                            (
                                                ptype,
                                                f"particle_angular_momentum_{ax}",
                                            )
                                        )
                        break
        else:
            # Default particle fields based on reading the binary header
            for ax in "xyz"[:ndim]:
                fields.append((ptype, f"particle_position_{ax}"))
            for ax in "xyz"[:ndim]:
                fields.append((ptype, f"particle_velocity_{ax}"))
            fields.append((ptype, "particle_mass"))
            fields.append((ptype, "particle_refinement_level"))
            fields.append((ptype, "particle_identity"))

        return fields

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            domains = [dom for dom in self.domains]
            base_region = getattr(dobj, "base_region", dobj)
            chunks = [
                MiniRAMSESDomainSubset(
                    base_region,
                    d,
                    self.dataset,
                    self.over_refine_factor,
                )
                for d in domains
            ]
            dobj._chunk_info = chunks
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [])
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g], None)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None)


class MiniRAMSESDataset(Dataset):
    _index_class = MiniRAMSESIndex
    _field_info_class = MiniRAMSESFieldInfo
    _dataset_type = "mini_ramses"

    # Attributes set during parsing
    gamma = 1.4
    nfile = 1

    def __init__(
        self,
        filename,
        dataset_type="mini_ramses",
        fields=None,
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        self.fluid_types += ("mini-ramses", "gravity")
        self._fields_in_file = []
        self._gravity_fields_in_file = []
        self.storage_filename = storage_filename
        self.default_species_fields = default_species_fields

        # Sanitize input filename
        sanitizer = MiniRAMSESFileSanitizer(filename)
        sanitizer.validate()
        self.root_folder = str(sanitizer.root_folder)

        super().__init__(
            str(sanitizer.info_fname),
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    @property
    def basename(self):
        return os.path.basename(self.root_folder)

    def _set_code_unit_attributes(self):
        """Set code unit attributes from info file parameters."""
        setdefaultattr(self, "length_unit", self.quan(self.unit_l, "cm"))
        setdefaultattr(self, "mass_unit", self.quan(
            self.unit_d * self.unit_l**3, "g"
        ))
        setdefaultattr(self, "time_unit", self.quan(self.unit_t, "s"))
        setdefaultattr(
            self, "velocity_unit",
            self.quan(self.unit_l / self.unit_t, "cm/s"),
        )
        setdefaultattr(
            self,
            "magnetic_unit",
            self.quan(
                np.sqrt(4.0 * np.pi * self.unit_d)
                * self.unit_l
                / self.unit_t,
                "gauss",
            ),
        )

    def _parse_parameter_file(self):
        """Parse the mini-ramses info.txt file."""
        info_path = os.path.join(self.root_folder, "info.txt")

        # Read all key=value pairs from info.txt
        params = {}
        with open(info_path) as f:
            for line in f:
                line = line.strip()
                if "=" in line:
                    key, val = line.split("=", 1)
                    key = key.strip()
                    val = val.strip()
                    try:
                        if "." in val or "E" in val.upper():
                            params[key] = float(val)
                        else:
                            params[key] = int(val)
                    except ValueError:
                        params[key] = val

        # Extract parameters
        # mini-ramses info.txt has: nfile, ncpu, ndim, levelmin, levelmax,
        # ngridmax, nstep_coarse, boxlen, time, texp, aexp, H0,
        # omega_m, omega_l, omega_k, omega_b, gamma, unit_l, unit_d, unit_t
        self.nfile = params.get("nfile", 1)
        self.dimensionality = params.get("ndim", 3)
        self.min_level = params.get("levelmin", 1)
        self.max_level = params.get("levelmax", 1)
        self.domain_dimensions = np.ones(3, dtype="int64") * 2 ** (
            self.min_level
        )

        self.gamma = params.get("gamma", 1.4)
        boxlen = params.get("boxlen", 1.0)

        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype="float64") * boxlen

        # Physical parameters
        self.unit_l = params.get("unit_l", 1.0)
        self.unit_d = params.get("unit_d", 1.0)
        self.unit_t = params.get("unit_t", 1.0)

        # Cosmology parameters
        self.omega_matter = params.get("omega_m", 0.0)
        self.omega_lambda = params.get("omega_l", 0.0)
        self.omega_radiation = 0.0
        self.hubble_constant = params.get("H0", 0.0) / 100.0  # Convert to H100
        current_time = params.get("time", 0.0)
        self.current_time = current_time

        aexp = params.get("aexp", 1.0)
        h0 = params.get("H0", 0.0)

        # Determine if this is a cosmological simulation
        # Same logic as RAMSES: if time >= 0 and H0 == 1 and aexp == 1
        # then it's NOT cosmological (all three must be true)
        is_cosmological = not (
            current_time >= 0 and h0 == 1.0 and aexp == 1.0
        )
        
        if is_cosmological:
            self.cosmological_simulation = True
            self.current_redshift = 1.0 / aexp - 1.0
        else:
            self.cosmological_simulation = False
            self.current_redshift = 0.0

        # Store all parameters
        self.parameters = params
        self.parameters["boxlen"] = boxlen
        self.parameters["time"] = current_time
        self.parameters["gamma"] = self.gamma

        self.unique_identifier = int(
            os.stat(self.parameter_filename).st_ctime
        )

        # Periodicity
        self._periodicity = (True, True, True)

        # Refine by 2 in each dimension
        self.refine_by = 2

        # Detect fluid fields (hydro and gravity)
        self._detect_fluid_fields()
        self._detect_gravity_fields()

        # Detect particle types
        self._detect_particle_types()

    def _detect_fluid_fields(self):
        """Detect available fluid fields from hydro header or first hydro file."""
        # Try reading hydro_header.txt first
        header_fn = os.path.join(self.root_folder, "hydro_header.txt")
        if os.path.exists(header_fn):
            self._fields_in_file = self._parse_hydro_header(header_fn)
            return

        # Try reading the first hydro file header
        hydro_fn = os.path.join(self.root_folder, "hydro.00001")
        if os.path.exists(hydro_fn):
            with open(hydro_fn, "rb") as f:
                ndim = struct.unpack("i", f.read(4))[0]
                nvar = struct.unpack("i", f.read(4))[0]
            self._fields_in_file = self._default_hydro_fields(nvar)
            return

        self._fields_in_file = []

    def _parse_hydro_header(self, header_fn):
        """Parse the mini-ramses hydro_header.txt file."""
        fields = []
        field_map = {
            "density": ("mini-ramses", "Density"),
            "velocity_x": ("mini-ramses", "x-velocity"),
            "velocity_y": ("mini-ramses", "y-velocity"),
            "velocity_z": ("mini-ramses", "z-velocity"),
            "thermal_pressure": ("mini-ramses", "Pressure"),
            "magnetic_field_x": ("mini-ramses", "B_x"),
            "magnetic_field_y": ("mini-ramses", "B_y"),
            "magnetic_field_z": ("mini-ramses", "B_z"),
            "metal_mass_fraction": ("mini-ramses", "Metallicity"),
        }
        with open(header_fn) as f:
            for line in f:
                line = line.strip()
                if line.startswith("nvar"):
                    continue
                if line.startswith("variable"):
                    # Parse "variable # N: name"
                    parts = line.split(":")
                    if len(parts) >= 2:
                        fname = parts[1].strip()
                        if fname in field_map:
                            fields.append(field_map[fname])
                        else:
                            fields.append(("mini-ramses", fname))
        return fields

    def _default_hydro_fields(self, nvar):
        """Return default field names based on nvar count."""
        ndim = self.dimensionality
        fields = [("mini-ramses", "Density")]
        for ax in "xyz"[:ndim]:
            fields.append(("mini-ramses", f"{ax}-velocity"))
        fields.append(("mini-ramses", "Pressure"))
        # MHD fields
        if nvar > ndim + 2:
            remaining = nvar - (ndim + 2)
            if remaining >= 3:
                for ax in "xyz":
                    fields.append(("mini-ramses", f"B_{ax}"))
                remaining -= 3
            for i in range(remaining):
                fields.append(("mini-ramses", f"scalar_{i + 1}"))
        return fields

    def _detect_gravity_fields(self):
        """Detect available gravity fields from grav header or first grav file."""
        # Try reading grav_header.txt first
        header_fn = os.path.join(self.root_folder, "grav_header.txt")
        if os.path.exists(header_fn):
            self._gravity_fields_in_file = self._parse_grav_header(header_fn)
            return

        # Try reading the first grav file header
        grav_fn = os.path.join(self.root_folder, "grav.00001")
        if os.path.exists(grav_fn):
            with open(grav_fn, "rb") as f:
                ndim = struct.unpack("i", f.read(4))[0]
                nvar = struct.unpack("i", f.read(4))[0]
            self._gravity_fields_in_file = self._default_gravity_fields(nvar)
            return

        self._gravity_fields_in_file = []

    def _parse_grav_header(self, header_fn):
        """Parse the mini-ramses grav_header.txt file."""
        fields = []
        field_map = {
            "potential": ("gravity", "Potential"),
            "accel_x": ("gravity", "x-acceleration"),
            "accel_y": ("gravity", "y-acceleration"),
            "accel_z": ("gravity", "z-acceleration"),
        }
        with open(header_fn) as f:
            for line in f:
                line = line.strip()
                if line.startswith("nvar"):
                    continue
                if line.startswith("variable"):
                    # Parse "variable # N: name"
                    parts = line.split(":")
                    if len(parts) >= 2:
                        fname = parts[1].strip()
                        if fname in field_map:
                            fields.append(field_map[fname])
                        else:
                            fields.append(("gravity", fname))
        return fields

    def _default_gravity_fields(self, nvar):
        """Return default gravity field names based on nvar count."""
        fields = [("gravity", "Potential")]
        for ax in "xyz"[:min(nvar - 1, 3)]:
            fields.append(("gravity", f"{ax}-acceleration"))
        return fields

    def _detect_particle_types(self):
        """Detect available particle types from output files."""
        particle_types = []
        ptype_map = {
            "part": "io",
            "star": "star",
            "sink": "sink",
        }
        for prefix, ptype in ptype_map.items():
            fn = os.path.join(self.root_folder, f"{prefix}.00001")
            if os.path.exists(fn):
                particle_types.append(ptype)

        if particle_types:
            self.particle_types = tuple(particle_types)
            self.particle_types_raw = tuple(particle_types)
        else:
            self.particle_types = ()
            self.particle_types_raw = ()

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if cls._missing_load_requirements():
            return False
        return MiniRAMSESFileSanitizer(filename).is_valid

    def __str__(self):
        return self.basename
