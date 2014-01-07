"""
Halo filter object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .halo_callbacks import HaloCallback

class HaloFilter(HaloCallback):
    def __init__(self, function, args, kwargs):
        HaloCallback.__init__(self, function, args, kwargs)

    def __call__(self, halo_catalog, halo):
        return self.function(halo_catalog, halo, *self.args, **self.kwargs)
