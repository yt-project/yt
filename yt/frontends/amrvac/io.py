"""
AMRVAC-specific IO functions



"""
import os

import numpy as np
from more_itertools import always_iterable

from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.on_demand_imports import _f90nml as f90nml


def read_amrvac_namelist(parfiles):
    """Read one or more parfiles, and return a unified f90nml.Namelist object.

    This function replicates the patching logic of MPI-AMRVAC where redundant parameters
    only retain last-in-line values, with the exception of `&filelist:base_filename`,
    which is accumulated. When passed a single file, this function acts as a mere
    wrapper of f90nml.read().

    Parameters
    ----------
    parfiles : str, os.Pathlike, byte, or an iterable returning those types
        A file path, or a list of file paths to MPI-AMRVAC configuration parfiles.

    Returns
    -------
    unified_namelist : f90nml.Namelist
        A single namelist object. The class inherits from ordereddict.

    """
    parfiles = (os.path.expanduser(pf) for pf in always_iterable(parfiles))

    # first merge the namelists
    namelists = [f90nml.read(parfile) for parfile in parfiles]
    unified_namelist = f90nml.Namelist()
    for nml in namelists:
        unified_namelist.patch(nml)

    # accumulate `&filelist:base_filename`
    base_filename = "".join(
        nml.get("filelist", {}).get("base_filename", "") for nml in namelists
    )
    unified_namelist["filelist"]["base_filename"] = base_filename

    return unified_namelist


class AMRVACIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "amrvac"

    def __init__(self, ds):
        BaseIOHandler.__init__(self, ds)
        self.ds = ds
        self.datfile = ds.parameter_filename
        header = self.ds.parameters
        self.block_shape = np.append(header["block_nx"], header["nw"])

    def _read_particle_coords(self, chunks, ptf):
        """Not implemented yet."""
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        """Not implemented yet."""
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        raise NotImplementedError

    def _read_data(self, fid, grid, field):
        """Retrieve field data from a grid.

        Parameters
        ----------
        fid: file descriptor (open binary file with read access)

        grid : yt.frontends.amrvac.data_structures.AMRVACGrid
            The grid from which data is to be read.
        field : str
            A field name.

        Returns
        -------
        data : np.ndarray
            A 3D array of float64 type representing grid data.

        """
        ileaf = grid.id
        offset = grid._index.block_offsets[ileaf]
        field_idx = self.ds.parameters["w_names"].index(field)

        field_shape = self.block_shape[:-1]
        count = np.prod(field_shape)
        byte_size_field = count * 8  # size of a double

        fid.seek(offset + byte_size_field * field_idx)
        data = np.fromfile(fid, "=f8", count=count)
        data.shape = field_shape[::-1]
        data = data.T
        # Always convert data to 3D, as grid.ActiveDimensions is always 3D
        while len(data.shape) < 3:
            data = data[..., np.newaxis]
        return data

    def _read_fluid_selection(self, chunks, selector, fields, size):
        """Retrieve field(s) data in a selected region of space.

        Parameters
        ----------
        chunks : generator
            A generator for multiple chunks, each of which contains a list of grids.

        selector : yt.geometry.selection_routines.SelectorObject
            A spatial region selector.

        fields : list
            A list of tuples (ftype, fname).

        size : np.int64
            The cumulative number of objs contained in all chunks.

        Returns
        -------
        data_dict : dict
            keys are the (ftype, fname) tuples, values are arrays that have been masked
            using whatever selector method is appropriate. Arrays have dtype float64.
        """

        # @Notes from Niels:
        # The chunks list has YTDataChunk objects containing the different grids.
        # The list of grids can be obtained by doing eg.
        # grids_list = chunks[0].objs or chunks[1].objs etc.
        # Every element in "grids_list" is then an AMRVACGrid object,
        # and has hence all attributes of a grid :
        # (Level, ActiveDimensions, LeftEdge, etc.)

        chunks = list(chunks)
        data_dict = {}  # <- return variable

        if isinstance(selector, GridSelector):
            if not len(chunks) == len(chunks[0].objs) == 1:
                raise RuntimeError
            grid = chunks[0].objs[0]

            with open(self.datfile, "rb") as fh:
                for ftype, fname in fields:
                    data_dict[ftype, fname] = self._read_data(fh, grid, fname)
        else:
            if size is None:
                size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)

            for field in fields:
                data_dict[field] = np.empty(size, dtype="float64")

            # nb_grids = sum(len(chunk.objs) for chunk in chunks)
            with open(self.datfile, "rb") as fh:
                ind = 0
                for chunk in chunks:
                    for grid in chunk.objs:
                        nd = 0
                        for field in fields:
                            ftype, fname = field
                            data = self._read_data(fh, grid, fname)
                            nd = grid.select(selector, data, data_dict[field], ind)
                        ind += nd

        return data_dict

    def _read_chunk_data(self, chunk, fields):
        """Not implemented yet."""
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        # it should be used by _read_fluid_selection instead of _read_data
        raise NotImplementedError
