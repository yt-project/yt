"""
TIGER-specific IO functions


Authors:
 * Matthew Turk 


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

class IOHandlerTiger(BaseIOHandler):
    _data_style = "tiger"
    _offset = 36

    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)
        self._memmaps = {}

    def _read_data(self, grid, field):
        fn = grid.pf.basename + grid.hierarchy.file_mapping[field]
        LD = np.array(grid.left_dims, dtype='int64')
        SS = np.array(grid.ActiveDimensions, dtype='int64')
        RS = np.array(grid.pf.root_size, dtype='int64')
        data = au.read_tiger_section(fn, LD, SS, RS).astype("float64")
        return data
