"""
AMRVAC-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.geometry.selection_routines import \
    GridSelector
import numpy as np
from .datreader import get_block_info, get_single_block_data


class AMRVACIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = 'amrvac'

    def __init__(self, ds):
        BaseIOHandler.__init__(self, ds)
        self.ds = ds
        self.datfile = ds.parameter_filename

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        raise NotImplementedError

    def _read_data(self, grid, field):
        ileaf = grid.id
        offset = grid._index.block_offsets[ileaf]
        block_shape = grid._index.block_shapes[ileaf]
        with open(self.datfile, "rb") as istream:
            block = get_single_block_data(istream, offset, block_shape)
        dim = self.ds.dimensionality
        field_idx = self.ds.parameters['w_names'].index(field)

        # Always convert into 3D array, as grid.ActiveDimensions is always 3D
        if dim == 1:
            data = block[:, field_idx]
            data = data[:, np.newaxis, np.newaxis]
        elif dim == 2:
            data = block[:, :, field_idx]
            data = data[:, :, np.newaxis]
        else:
            data = block[:, :, :, field_idx]
        return data

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # This needs to allocate a set of arrays inside a dictionary, where the
        # keys are the (ftype, fname) tuples and the values are arrays that
        # have been masked using whatever selector method is appropriate.  The
        # dict gets returned at the end and it should be flat, with selected
        # data.  Note that if you're reading grid data, you might need to
        # special-case a grid selector object.
        # Also note that "chunks" is a generator for multiple chunks, each of
        # which contains a list of grids. The returned numpy arrays should be
        # in 64-bit float and contiguous along the z direction. Therefore, for
        # a C-like input array with the dimension [x][y][z] or a
        # Fortran-like input array with the dimension (z,y,x), a matrix
        # transpose is required (e.g., using np_array.transpose() or
        # np_array.swapaxes(0,2)).
        data_dict = {}
        chunks = list(chunks)
        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            grid = chunks[0].objs[0]
            for ftype, fname in fields:
                data_dict[ftype, fname] = self._read_data(grid, fname)
            return data_dict

        '''
        @Notes from Niels:
        The chunks list has YTDataChunk objects containing the different grids.
        The list of grids can be obtained by doing eg. grids_list = chunks[0].objs or chunks[1].objs etc.
        Every element in "grids_list" is then an AMRVACGrid object, and has hence all attributes of a grid
            (Level, ActiveDimensions, LeftEdge, etc.)
        '''
        
        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))

        for field in fields:
            data_dict[field] = np.empty(size, dtype='float64')

        nb_grids = sum(len(chunk.objs) for chunk in chunks)

        ind = 0
        for chunk in chunks:
            for grid in chunk.objs:
                nd = 0
                for field in fields:
                    ftype, fname = field
                    data = self._read_data(grid, fname)
                    nd = grid.select(selector, data, data_dict[field], ind)
                ind += nd

        return data_dict


    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        # it should be used by _read_fluid_selection instead of _read_data
        raise NotImplementedError
