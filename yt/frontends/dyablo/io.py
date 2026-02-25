"""
IO handler for Dyablo frontend.
"""

from collections.abc import Iterable

import h5py
import numpy as np

from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.selection_routines import SelectorObject
from yt.utilities.io_handler import BaseIOHandler


class DyabloIOHandler(BaseIOHandler):
    """IO handler for reading Dyablo data."""

    _dataset_type = "dyablo"

    def _read_fluid_selection(
        self,
        chunks: Iterable[YTDataChunk],
        selector: SelectorObject,
        fields: Iterable[tuple],
        size: int,
    ):
        """
        Read fluid data from blocks selected by selector.

        Parameters
        ----------
        chunks : iterable
            Iterable of chunk information (geometry chunks)
        selector : yt.geometry.selection_routines.SelectorObject
            SelectorObject object to filter cells
        fields : iterable
            Iterable of (ftype, fname) tuples to read
        size : int
            Expected size of data

        Returns
        -------
        dict
            Dictionary mapping (ftype, fname) to arrays of selected data
        """
        from collections import defaultdict

        rv = defaultdict(list)
        fields = list(fields)

        # Iterate through chunks and subsets (similar to RAMSES)
        for chunk in chunks:
            for subset in chunk.objs:
                # Use the fill method which queries the octree
                data = subset.fill(fields, selector)

                for field in fields:
                    rv[field].append(data[field])

        # Concatenate all data for each field
        return {
            field: np.concatenate(rv[field]) if rv[field] else np.array([])
            for field in fields
        }

    def _read_chunk_data(
        self,
        chunk: YTDataChunk,
        fields: Iterable[tuple],
    ):
        """
        Read a full chunk of data and cache it (optional optimization).

        Parameters
        ----------
        chunk : yt.geometry.geometry_handler.YTDataChunk
            Chunk to read
        fields : iterable
            Fields to read

        Returns
        -------
        dict
            Dictionary mapping fields to data arrays
        """
        # This is an optional optimization for caching
        # For now, we'll just pass
        pass

    def _read_particle_coords(
        self,
        chunks: Iterable[YTDataChunk],
        ptf: dict[str, list[str]],
    ):
        """
        Read particle coordinates.

        Parameters
        ----------
        chunks : iterable
            Iterable of chunk information
        ptf : dict
            Dictionary mapping particle types to lists of fields

        Yields
        ------
        tuple
            (ptype, (x, y, z)) for each particle type in ptf
        """
        filename = self.ds._particle_filename
        if filename is None or not filename:
            return

        # Handle single file or list of files
        if isinstance(filename, str):
            filename = [filename]

        for ptype in ptf:
            coords_list = []
            for fname in filename:
                try:
                    with h5py.File(fname, "r") as f:
                        if "/coordinates" in f:
                            data = f["/coordinates"][:]
                            # Extract x, y, z coordinates
                            if data.ndim == 2 and data.shape[1] >= 3:
                                coords_list.append(data)
                except Exception:
                    # File might not exist or be unreadable
                    continue

            if coords_list:
                coords = np.concatenate(coords_list)
                x = coords[:, 0]
                y = coords[:, 1]
                z = coords[:, 2]
                yield (ptype, (x, y, z))

    def _read_particle_fields(
        self,
        chunks: Iterable[YTDataChunk],
        ptf: dict[str, list[str]],
        selector: SelectorObject,
    ):
        """
        Read particle field data.

        Parameters
        ----------
        chunks : iterable
            Iterable of chunk information
        ptf : dict
            Dictionary mapping particle types to lists of fields
        selector : yt.geometry.selection_routines.SelectorObject
            SelectorObject object for filtering particles

        Yields
        ------
        tuple
            ((ptype, field), data_array) for each field in ptf
        """
        filename = self.ds._particle_filename
        if filename is None or not filename:
            return

        # Handle single file or list of files
        if isinstance(filename, str):
            filename = [filename]

        for ptype in ptf:
            for field_name in ptf[ptype]:
                data_list = []
                for fname in filename:
                    try:
                        with h5py.File(fname, "r") as f:
                            # Try different field name patterns
                            coords_path = "/coordinates"
                            mass_path = "/mass"
                            stripped = field_name.replace("particle_", "")
                            stripped_path = f"/{stripped}"

                            field_candidates = [
                                f"/{field_name}",
                                stripped_path,
                                coords_path if "position" in field_name else None,
                                mass_path if field_name == "mass" else None,
                            ]

                            for field_path in field_candidates:
                                if field_path and field_path in f:
                                    field_data = f[field_path][:]
                                    # Handle position fields
                                    if "position" in field_name and field_data.ndim > 1:
                                        axis_map = {"x": 0, "y": 1, "z": 2}
                                        axis = axis_map.get(
                                            field_name.split("_")[-1], 0
                                        )
                                        if axis < field_data.shape[1]:
                                            field_data = field_data[:, axis]
                                    data_list.append(field_data.ravel())
                                    break
                    except Exception:
                        continue

                if data_list:
                    field_data = np.concatenate(data_list)
                    yield ((ptype, field_name), field_data)
