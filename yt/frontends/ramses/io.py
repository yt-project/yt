"""
RAMSES-specific IO


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

from collections import defaultdict
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerRAMSES(BaseIOHandler):
    _data_style = "ramses"

    def __init__(self, ramses_tree, *args, **kwargs):
        self.ramses_tree = ramses_tree
        BaseIOHandler.__init__(self, *args, **kwargs)

    def _read_data(self, grid, field):
        tr = np.zeros(grid.ActiveDimensions, dtype='float64')
        filled = np.zeros(grid.ActiveDimensions, dtype='int32')
        to_fill = grid.ActiveDimensions.prod()
        grids = [grid]
        l_delta = 0
        varindex = self.ramses_tree.field_ind[field]
        while to_fill > 0 and len(grids) > 0:
            next_grids = []
            for g in grids:
                to_fill -= self.ramses_tree.read_grid(varindex, field,
                        grid.get_global_startindex(), grid.ActiveDimensions,
                        tr, filled, g.Level, 2**l_delta, g.locations)
                next_grids += g.Parent
            grids = next_grids
            l_delta += 1
        return tr

    def preload(self, grids, sets):
        if len(grids) == 0: return
        domain_keys = defaultdict(list)
        pf_field_list = grids[0].pf.h.field_list
        sets = [dset for dset in list(sets) if dset in pf_field_list]
        exc = self._read_exception
        for g in grids:
            domain_keys[g.domain].append(g)
        for domain, grids in domain_keys.items():
            mylog.debug("Starting read of domain %s (%s)", domain, sets)
            for field in sets:
                for g in grids:
                    self.queue[g.id][field] = self._read_data(g, field)
                print "Clearing", field, domain
                self.ramses_tree.clear_tree(field, domain - 1)
        mylog.debug("Finished read of %s", sets)

    def modify(self, data): return data
