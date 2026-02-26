"""
Data structures for Dyablo frontend.
"""

from pathlib import Path

import h5py
import numpy as np

from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.static_output import Dataset
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import OctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.file_handler import HDF5FileHandler
from yt.utilities.logger import ytLogger as mylog

from .definitions import HYDRO_FILE_PATTERN, PARTICLE_FILE_PATTERN
from .fields import DyabloFieldInfo
from .io import DyabloIOHandler


class DyabloOctreeIndex(OctreeIndex):
    """Octree Index for Dyablo data with leaf-only blocks."""

    domain_id = 1  # Dyablo uses a single domain

    def __init__(self, ds, dataset_type):
        self.dataset_type = dataset_type
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = Path(self.index_filename).parent
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _initialize_oct_handler(self):
        """Initialize octree handler from Dyablo connectivity/coordinates."""
        mylog.debug("Initializing Dyablo Octree Handler")

        # Get connectivity and coordinates from dataset
        connectivity = self.ds._connectivity
        coordinates = self.ds._coordinates
        n_cells = self.ds._n_cells
        block_size = self.ds.block_size  # (N, M, L)

        # Store for later use by IO handler
        self.connectivity = connectivity
        self.coordinates = coordinates
        self.n_cells = n_cells
        self.block_size = block_size

        # Calculate number of blocks
        N, M, L = block_size
        cells_per_block = N * M * L
        n_blocks = n_cells // cells_per_block

        if n_cells % cells_per_block != 0:
            mylog.warning(
                f"Number of cells ({n_cells}) is not divisible by "
                f"cells per block ({cells_per_block}). "
                "Some cells will be ignored."
            )
        self._n_blocks = n_blocks
        self._n_cells = n_cells
        self._cells_per_block = cells_per_block

        # OctreeContainer needs the number of root octs (= n_blocks),
        # which is domain_dimensions / block_size.
        n_oct_dims = self.ds.domain_dimensions // np.array(block_size, dtype=np.int32)
        self.oct_handler = OctreeContainer(
            n_oct_dims,
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
            num_zones=block_size,
        )

        # Compute block centers and levels
        # Each block contains cells_per_block cells in C-order
        block_centers = np.zeros((n_blocks, 3), dtype=np.float64)
        block_levels = np.zeros(n_blocks, dtype=np.uint64)

        # First pass: find the maximum block width (coarsest level)
        max_block_width = 0.0
        for block_id in range(n_blocks):
            # Get cell indices for this block
            cell_start = block_id * cells_per_block
            cell_end = cell_start + cells_per_block

            # Get corners for all cells in this block
            block_connectivity = connectivity[cell_start:cell_end]
            block_corners = coordinates[block_connectivity, :]

            # Compute block bounds from all cell corners
            # (cells_per_block * 8, 3)
            all_corners = block_corners.reshape(-1, 3)
            block_min = np.min(all_corners, axis=0)
            block_max = np.max(all_corners, axis=0)
            block_width = block_max - block_min

            # Track maximum block width
            chunk_max_width = np.max(block_width)
            max_block_width = max(max_block_width, chunk_max_width)

        # Second pass: compute centers and levels
        for block_id in range(n_blocks):
            # Get cell indices for this block
            cell_start = block_id * cells_per_block
            cell_end = cell_start + cells_per_block

            # Get corners for all cells in this block
            block_connectivity = connectivity[cell_start:cell_end]
            block_corners = coordinates[
                block_connectivity, :
            ]  # (cells_per_block, 8, 3)

            # Compute block center: mean of all cell centers
            # (cells_per_block, 3)
            cell_centers = np.mean(block_corners, axis=1)
            block_centers[block_id] = np.mean(cell_centers, axis=0)

            # Compute block size and level
            all_corners = block_corners.reshape(-1, 3)
            block_min = np.min(all_corners, axis=0)
            block_max = np.max(all_corners, axis=0)
            block_width = block_max - block_min

            # Level is determined by block size relative to coarsest block
            max_dim_width = np.max(block_width)
            refinement_ratio = max_block_width / max_dim_width
            level = np.round(np.log2(refinement_ratio)).astype(np.uint64)
            block_levels[block_id] = level

        # Count number of blocks - note that a block at level l
        # requires its parent blocks at levels < l to be present
        # so we need to add them as well
        levels, counts = np.unique(block_levels, return_counts=True)
        n_blocks_with_parents = 0
        for lvl, count in zip(levels, counts, strict=True):
            for i in range(lvl + 1):
                n_blocks_with_parents += count / (8 ** (lvl - i))

        n_blocks_with_parents = int(np.ceil(n_blocks_with_parents))
        # Allocate space for blocks
        self.oct_handler.allocate_domains([n_blocks_with_parents])

        # Add all blocks as octs at their appropriate levels
        # Each block gets file_ind from 0 to n_blocks-1
        self.oct_handler.add(
            1,  # domain 1
            -1,  # use levels array
            block_centers,  # positions
            levels=block_levels,
        )

        self.oct_handler.finalize()

    def _detect_output_fields(self):
        """Detect available fields in the dataset."""
        dsl = set()

        # Open the hydro file and look for available fields
        try:
            with h5py.File(self.dataset.parameter_filename, "r") as f:
                for key in f.keys():
                    # Skip metadata/coordinate datasets
                    skip_keys = (
                        "metadata",
                        "scalar_data",
                        "connectivity",
                        "coordinates",
                        "units",
                    )
                    if key in skip_keys:
                        continue

                    obj = f[key]
                    if isinstance(obj, h5py.Dataset):
                        dsl.add(("dyablo", key))
        except Exception as e:
            mylog.warning(f"Error detecting output fields: {e}")

        self.fluid_field_list = list(dsl)

        # Detect particle fields from particle files
        particle_field_list = []
        particle_types = set()
        particle_files = self.dataset._particle_filename or []
        for pfile in particle_files:
            pmatch = PARTICLE_FILE_PATTERN.match(Path(pfile).name)
            if not pmatch:
                continue
            ptype = pmatch.group(2)
            try:
                with h5py.File(pfile, "r") as f:
                    particle_field_list.append((ptype, "particle_position_x"))
                    particle_field_list.append((ptype, "particle_position_y"))
                    particle_field_list.append((ptype, "particle_position_z"))
                    for key in (_ for _ in f.keys() if isinstance(f[_], h5py.Dataset)):
                        if key == "coordinates":
                            continue
                        yt_key = f"particle_{key}"
                        particle_field_list.append((ptype, yt_key))
                        particle_types.add(ptype)
            except Exception as e:
                mylog.warning(f"Error reading particle file {pfile}: {e}")

        # Ensure all detected particles have been registered
        for ptype in (
            _ for _ in particle_types if _ not in self.dataset.particle_types
        ):
            raise RuntimeError(
                f"Detected particle type '{ptype}' from file {pfile} "
                "but it was not registered in dataset.particle_types. "
                "Please ensure the particle file is named correctly and "
                "that _parse_parameter_file is correctly detecting particle types."
            )

        self.particle_field_list = particle_field_list
        self.field_list = self.fluid_field_list + self.particle_field_list

    @property
    def max_level(self):
        """Return the maximum refinement level."""
        return self.ds.max_level

    def _identify_base_chunk(self, dobj):
        """Identify the chunks that make up a data object."""
        if getattr(dobj, "_chunk_info", None) is None:
            base_region = getattr(dobj, "base_region", dobj)
            # For Dyablo, we have a single domain containing all octs
            subset = [
                DyabloOctreeSubset(
                    base_region,
                    self,
                    self.dataset,
                )
            ]
            dobj._chunk_info = subset
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        """Return all chunks."""
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        """Yield spatial chunks."""
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g], None)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        """Yield IO chunks."""
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None, cache=cache)


class DyabloOctreeSubset(OctreeSubset):
    """Octree subset for Dyablo data with NxMxL blocks.

    Unlike standard octrees where each oct contains 2x2x2 cells,
    Dyablo uses blocks as octs, where each block contains NxMxL cells.
    This class overrides coordinate methods to properly handle this
    non-standard structure.
    """

    _domain_offset = 1
    _block_order = "C"

    def __init__(self, base_region, domain, ds):
        super().__init__(base_region, domain, ds, num_zones=ds.block_size)

        self._current_particle_type = "all"
        self._current_fluid_type = self.ds.default_fluid_type

    @property
    def oct_handler(self):
        return self.domain.oct_handler

    def fill(self, fields, selector):
        """
        Fill field data for selected cells.

        This method queries the octree for which blocks are selected,
        reads all cells from those blocks, then applies the selector
        at the cell level.

        Parameters
        ----------
        fields : list
            List of (ftype, fname) tuples to read
        selector : SelectorObject
            Selector object for filtering cells

        Returns
        -------
        dict
            Dictionary mapping (ftype, fname) to selected data arrays
        """
        oct_handler = self.oct_handler
        cell_count = selector.count_oct_cells(oct_handler)

        # Initialize data container
        data = {field: np.zeros(cell_count, "float64") for field in fields}

        # Early exit if no cells are selected
        if cell_count == 0:
            return data

        _level_inds, cell_inds, file_inds = oct_handler.file_index_octs(
            selector, 1, cell_count
        )

        fname = self.ds._hydro_filename
        cells_per_block = self.index._cells_per_block
        block_size = self.ds.block_size  # (N, M, L)

        # The oct visitor returns cell_inds in C-order (ix slowest, iz fastest):
        #   cell_ind = ix*M*L + iy*L + iz
        # The HDF5 file stores cells in Fortran-order (ix fastest, iz slowest):
        #   file_cell_ind = ix + N*iy + N*M*iz
        # Remap by decoding C-order then re-encoding as F-order.
        ix, iy, iz = np.unravel_index(cell_inds, block_size, order="C")
        file_cell_inds = np.ravel_multi_index((ix, iy, iz), block_size, order="F")

        # Flat HDF5 index: block_offset + cell_within_block
        indices = file_inds * cells_per_block + file_cell_inds

        with h5py.File(self.ds._hydro_filename, "r") as f:
            for ftype, fname in fields:
                field_path = f"/{fname}"
                if field_path in f:
                    # Apply cell-level mask
                    data[ftype, fname] = f[field_path][:][indices].astype("float64")
                else:
                    raise KeyError(
                        f"Field {field_path} not found in HDF5 file {self.ds._hydro_filename}"
                    )

        return data


class DyabloDataset(Dataset):
    """Dataset class for Dyablo code output files."""

    _index_class = DyabloOctreeIndex
    _field_info_class = DyabloFieldInfo
    _file_class = HDF5FileHandler
    _io_class = DyabloIOHandler

    def __init__(
        self,
        filename,
        dataset_type="dyablo",
        storage_filename=None,
        units_override=None,
        unit_system="code",
        # Metadata overrides for when they're not in the file
        domain_left_edge=None,
        domain_right_edge=None,
        domain_dimensions=None,
        periodicity=None,
        block_size=None,
        max_level=None,
    ):
        """
        Initialize a Dyablo dataset.

        Parameters
        ----------
        filename : str
            Path to the Dyablo HDF5 hydro output file
        dataset_type : str, optional
            Type of dataset, default is "dyablo"
        storage_filename : str, optional
            Name for the storage file (if using h5 storage backend)
        units_override : dict, optional
            Dictionary of unit overrides
        unit_system : str, optional
            Unit system to use, default is "code"
        domain_left_edge : array-like, optional
            Override for domain left edge
        domain_right_edge : array-like, optional
            Override for domain right edge
        domain_dimensions : array-like, optional
            Override for domain dimensions
        periodicity : tuple, optional
            Periodicity for each dimension (default: (False, False, False))
        block_size : tuple, optional
            Block size (N, M, L) - inferred from data if not provided
        max_level : int, optional
            Maximum refinement level - inferred from data if not provided
        """
        self.domain_left_edge_override = domain_left_edge
        self.domain_right_edge_override = domain_right_edge
        self.domain_dimensions_override = domain_dimensions
        self.periodicity_override = periodicity
        self.block_size_override = block_size
        self.max_level_override = max_level

        super().__init__(
            filename,
            dataset_type=dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.storage_filename = storage_filename

    @staticmethod
    def _infer_block_structure(connectivity, coordinates, n_cells):
        """
        Infer block size (N, M, L) and number of blocks from spatial cell
        layout.

        Assumes cells are C-ordered within blocks: cells advance in
        x (fastest), then y, then z (slowest).

        Parameters
        ----------
        connectivity : ndarray
            Connectivity array (n_cells, 8)
        coordinates : ndarray
            Coordinates array (n_vertices, 3)
        n_cells : int
            Total number of cells

        Returns
        -------
        tuple
            (N, M, L) block size
        tuple
            (N_blocks, M_blocks, L_blocks) number of blocks in each dimension
        """
        # Get first cell position and size
        corners_0 = coordinates[connectivity[0, :], :]
        pos_0 = np.mean(corners_0, axis=0)
        cell_size = np.max(corners_0, axis=0) - np.min(corners_0, axis=0)

        NML = np.ones(3, dtype=int)

        jump_size = 1
        for idim in range(3):
            spacing = [0, 0, 0]
            spacing[idim] = cell_size[idim]

            i = jump_size
            ii = 1
            prev_pos = pos_0
            # Keep on iterating so long as the position change
            # is one dx along current dimension
            while i < n_cells:
                corners_i = coordinates[connectivity[i, :], :]
                pos = np.mean(corners_i, axis=0)
                # Changed in another dimension, so we're done
                if not np.isclose(pos - prev_pos, spacing).all():
                    NML[idim] = ii
                    jump_size *= ii
                    break
                prev_pos = pos
                ii += 1
                i += jump_size
        # Now we have the number of cells/block.
        # Compute number of blocks in each dimension
        domain_size = np.max(coordinates, axis=0) - np.min(coordinates, axis=0)

        # Find smallest largest cell size
        cells_per_block = np.prod(NML)
        top_left_cells = coordinates[connectivity[::cells_per_block, :], :]
        cell_sizes = np.ptp(top_left_cells, axis=1)

        largest_cell = np.max(cell_sizes, axis=0)

        n_blocks = np.ceil(domain_size / (largest_cell * NML)).astype(int)
        return NML, n_blocks

    def _parse_parameter_file(self):
        """Parse metadata from the HDF5 file."""
        # Open the hydro file
        with h5py.File(self.parameter_filename, "r") as f:
            # Try to read metadata
            metadata = {}
            if "/metadata" in f:
                metadata = dict(f["/metadata"].attrs)
            if "/scalar_data" in f:
                scalar_attrs = dict(f["/scalar_data"].attrs)
                metadata.update(scalar_attrs)

            coordinates = f["/coordinates"][:]
            connectivity = f["/connectivity"][:]

        # Also store simulation name into metadata
        sim_name = HYDRO_FILE_PATTERN.match(Path(self.parameter_filename).name).group(1)
        metadata["name"] = sim_name

        self.parameters.update(metadata)

        self.root_folder = Path(self.parameter_filename).parent

        n_cells = connectivity.shape[0] if len(connectivity.shape) > 0 else 1

        # Infer domain bounds from coordinate extrema
        if self.domain_left_edge_override is not None:
            domain_left = self.domain_left_edge_override
        else:
            domain_left = np.min(coordinates, axis=0)

        if self.domain_right_edge_override is not None:
            domain_right = self.domain_right_edge_override
        else:
            domain_right = np.max(coordinates, axis=0)

        # Get periodicity
        if "periodicity" in metadata:
            periodicity = tuple(bool(metadata["periodicity"][i]) for i in range(3))
        elif self.periodicity_override is not None:
            periodicity = self.periodicity_override
        else:
            periodicity = (False, False, False)

        # Get block size - try metadata first, then infer from connectivity
        if "block_size" in metadata:
            block_size = tuple(int(x) for x in metadata["block_size"][:3])
        elif self.block_size_override is not None:
            block_size = self.block_size_override
        else:
            # Infer from spatial layout
            block_size, n_blocks = self._infer_block_structure(
                connectivity, coordinates, n_cells
            )
            mylog.info(
                "Inferred block size: %s, number of blocks: %s",
                block_size,
                n_blocks,
            )

        # Infer min level
        min_level = np.max(np.log2(block_size).astype(int))

        # Get max refinement level
        if "max_level" in metadata:
            max_level = int(metadata["max_level"])
        elif self.max_level_override is not None:
            max_level = self.max_level_override
        else:
            max_level = min_level  # Default: no refinement

        current_time = float(metadata.get("time", 0.0))

        # Set dataset parameters
        self.domain_left_edge = np.array(domain_left, dtype=np.float64)
        self.domain_right_edge = np.array(domain_right, dtype=np.float64)
        # domain_dimensions must be total cells at the coarsest level so that
        # the QuadTree projection (which uses cell icoords = pos*nz + ind) does
        # not overflow. OctreeContainer receives n_blocks (= domain_dimensions //
        # block_size) as its oct_domain_dimensions.
        self.domain_dimensions = np.array(n_blocks, dtype=np.int32) * np.array(
            block_size, dtype=np.int32
        )

        self._periodicity = periodicity
        self.refine_by = 2  # Dyablo uses factor-of-2 refinement
        self.current_time = current_time
        self.current_redshift = 0
        self.cosmological_simulation = False
        self.min_level = min_level
        self.max_level = max_level
        self.block_size = block_size
        # If we only have one block in a given dimension,
        # decrease dimensionality
        self.dimensionality = 3 - (n_blocks == 1).sum()

        # Store connectivity and coordinates for later use
        self._connectivity = connectivity
        self._coordinates = coordinates
        self._n_cells = n_cells

        # Store file information
        self._hydro_filename = self.parameter_filename
        self._particle_filename = self._find_particle_file()

        # Set particle_types from detected particle files so the 'all' union
        # is built correctly in create_field_info (before it runs).
        ptypes = []
        for pfile in self._particle_filename or []:
            m = PARTICLE_FILE_PATTERN.match(Path(pfile).name)
            if m:
                ptype = m.group(2)
                if ptype not in ptypes:
                    ptypes.append(ptype)
        if ptypes:
            self.particle_types = self.particle_types_raw = tuple(ptypes)

    def _find_particle_file(self):
        """Find and validate particle file matching the hydro file."""
        hydro_match = HYDRO_FILE_PATTERN.match(Path(self.parameter_filename).name)
        if not hydro_match:
            return None

        name, iteration = hydro_match.groups()
        directory = Path(self.parameter_filename).parent

        # Look for particle files with matching iteration
        pattern = directory.glob(f"{name}_particles_*_iter{iteration}.h5")
        particle_files = sorted(pattern)

        if particle_files:
            return particle_files  # Return list of particle files
        return None

    @staticmethod
    def _is_valid(filename, *args, **kwargs) -> bool:
        """Check if file is a valid Dyablo output."""
        if not Path(filename).exists():
            return False

        if not Path(filename).suffix == ".h5":
            return False

        # Check if it looks like a Dyablo hydro file
        if not HYDRO_FILE_PATTERN.match(Path(filename).name):
            return False

        # Verify it's an HDF5 file with expected Dyablo datasets
        try:
            with h5py.File(filename, "r") as f:
                # Check for key Dyablo datasets
                required_fields = ["/rho", "/coordinates"]
                if any(field not in f for field in required_fields):
                    return False
            return True
        except (OSError, Exception):
            return False

    def _set_code_unit_attributes(self):
        """Set code unit attributes."""
        # Try to read units from file, otherwise use defaults
        try:
            with h5py.File(self._hydro_filename, "r") as f:
                if "/units" in f:
                    units_data = f["/units"][:]
                    # Assume format: [length_unit, mass_unit, time_unit, ...]
                    if len(units_data) >= 3:
                        length_unit = float(units_data[0])
                        mass_unit = float(units_data[1])
                        time_unit = float(units_data[2])
                    else:
                        raise ValueError("Not enough units data")
                else:
                    raise ValueError("No units found in file")
        except Exception:
            # Fall back to defaults
            mylog.warning(
                "Could not read units from HDF5 file. Using default code units."
            )
            length_unit = 1.0
            mass_unit = 1.0
            time_unit = 1.0

        # Set units
        self.length_unit = self.quan(length_unit, "cm")
        self.mass_unit = self.quan(mass_unit, "g")
        self.time_unit = self.quan(time_unit, "s")
        self.velocity_unit = self.length_unit / self.time_unit

        # Register dyablo as a fluid type
        self.fluid_types += ("dyablo",)
